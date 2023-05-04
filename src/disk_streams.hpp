// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//    http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef SRC_CPP_DISK_STREAMS_HPP_
#define SRC_CPP_DISK_STREAMS_HPP_
#include "bitfield.hpp"
#include "bitfield_index.hpp"
#include "disk.hpp"
#include "pos_constants.hpp"
#include "util.hpp"
#include "stream_buffer.hpp"


struct IBlockWriter{
	virtual void Write( StreamBuffer & block ) = 0;
	virtual void Flush() {};
	virtual void Close() = 0;
	virtual ~IBlockWriter() = default;
};

struct IBlockReader{
	virtual uint32_t Read( StreamBuffer & block ) = 0;
	virtual bool atEnd() const = 0;
	virtual void Close() = 0;
	virtual ~IBlockReader() = default;
};

struct IBlockWriterReader : public IBlockWriter, IBlockReader {
	virtual void Remove() = 0;
};

// ============== File ========================
struct BlockedFileStream : public IBlockWriterReader {
	const int32_t read_align_size;
	BlockedFileStream( const std::string &fileName, uint32_t read_align_size = 1 )
		: read_align_size(read_align_size), file_name(fileName)	{
		assert( read_align_size > 0 );
	}

	void Write( StreamBuffer & block ) override {
		assert( read_position == 0 ); // do not read and write same time
		if( !disk ){
			assert( write_position == 0 );
			disk.reset( new FileDisk(file_name) );
		}
		disk->Write( write_position, block.get(), block.used() );
		write_position += block.used();
	}

	uint32_t Read( StreamBuffer & block ) override {
		uint32_t to_read = std::min( BUF_SIZE/read_align_size*read_align_size
																 , write_position - read_position );
		block.ensureSize( to_read ).setUsed( to_read );
		if( to_read > 0 ){
			assert( read_position < write_position );

			disk->Read( read_position, block.get(), to_read );
			read_position += to_read;
			assert( read_position <= write_position );
		}

		return to_read;
	}

	bool atEnd() const override{ return read_position >= write_position; };
	void Flush() override { if( disk ) disk->Flush(); };
	void Close() override { if( disk ) disk->Close(); };
	void Remove() override { if(disk) disk->Remove( true ); }
	~BlockedFileStream() { Remove(); }
private:
	std::unique_ptr<FileDisk> disk;
	fs::path file_name;
	uint64_t write_position = 0;
	uint64_t read_position = 0;
};


// ============== Buffer ========================
struct BlockBufferedWriter : public IBlockWriter{
	BlockBufferedWriter( IBlockWriter* disk, uint16_t entry_size = 1 )
		: disk(disk),write_size(BUF_SIZE/entry_size*entry_size ){}

	inline void Write( StreamBuffer & block ) override {
		assert( buffer.used() < write_size );
		for( uint32_t processed = 0; processed < block.used(); ){
			uint32_t to_process = std::min( write_size - buffer.used(), block.used() - processed );
			buffer.add( block.get() + processed, to_process );
			processed += to_process;

			assert( buffer.used() <= write_size );

			if( buffer.used() == write_size ){
				disk->Write( buffer );
				buffer.ensureSize( BUF_SIZE );
			}
		}
	};

	void Close() override {
		if( disk ){
			if( buffer.used() > 0 )
				disk->Write( buffer );

			disk->Close();
			disk.reset();
			buffer.reset(); // clear ram
		}
	};
	~BlockBufferedWriter(){ Close(); }

	StreamBuffer buffer;
private:
	std::unique_ptr<IBlockWriter> disk;
	const uint32_t write_size;
};

// ============== Cache ========================


/* CachedFileStream class stores in cache as many data as possible before wrtting it to file
	Warning! CachedFileStream class is not thread safe i.e. IO from more than one thread leads to unpredictable.
	Warning! Writing after start to read can lead to unperdictable.
*/
struct BlockCachedFile: IBlockWriterReader, ICacheConsumer {

	// @max_buf_size - is a size of buffer supposed to be full. if such side provided it is not copied but replaced
	BlockCachedFile( const fs::path &fileName, MemoryManager &memory_manager, uint32_t read_align = 1 )
		: memory_manager(memory_manager), cache( BUF_SIZE/read_align*read_align), read_align(read_align)
	{
		consumer_idx = memory_manager.registerConsumer( (ICacheConsumer*)this );
		file_name = fileName;
	}

	void Write( StreamBuffer & block ) override{
		assert( read_position == 0 ); // cannot read and write same time
		uint32_t block_position = 0;
		uint32_t block_used = block.used();

		if( consumer_idx >= 0 // we are caching
				&& cache_sync.try_lock() // and we can lock the cache
				&& consumer_idx >= 0 ) // and after lock we still caching :)
		{
			assert( (block_used%read_align) == 0 ); // check new data alligned

			// check we can add to current cache
			block_position = std::min( cache.bufSize() - cache.buf_size, block_used );
			if( block_position > 0 ){ // we really can fill last buffer of current cache
				memcpy( cache.buffer() + cache.bufSize(), block.get(), block_position );
				write_position += block_position;
				if( block_position == block_used ) return; // all data stored in the cache
			}

			auto buf = memory_manager.consumerRequest(); // request next buffer for cache
			if( buf != nullptr ){ // the buffer was provided
				if( block_position == 0 && block.size() == BUF_SIZE ){ // simplest case - just replace the buffer
					cache.add( block.release(buf), write_position, block_used );
					write_position += block_used;
					cache_sync.unlock();
					return;
				} else { // need to copy from current block to cache
					do{
						uint32_t to_copy = std::min( (uint32_t)BUF_SIZE, block_used - block_position );
						memcpy( buf, block.get() + block_position, to_copy );
						cache.add( buf, write_position, to_copy ); // add next buffer to cache
						write_position += to_copy;
						block_position += to_copy;
						buf = block_position >= block_used ? nullptr : memory_manager.consumerRequest(); // request next buffer
					}while( buf != nullptr ); // contuinue up to no next buffer or no need in next buffer
					cache_sync.unlock();
					if( block_position < block_used ){ // not all block was stored in cache
						write_position += DiskWrite( write_position, block.get()+block_position, block_used - block_position );
					}
					return;
				}
			}
			cache_sync.unlock();
		}

		write_position += DiskWrite( write_position, block.get() + block_position, block_used - block_position );
	}

	// warning this is not thread safe function
	uint32_t Read( StreamBuffer & block ) override{

		// first read we need to lock
		if( read_position == 0 ) cache_sync.lock();

		if( consumer_idx >= 0  ){
			// in order prevent locks on read we unregister started to read stream
			memory_manager.unregisterConsumer( this, consumer_idx );
			consumer_idx = -1;
		}

		uint32_t to_read = cache.buf_size;

		if( !cache.isEmpty() ) { // we have cached buffers
			assert( read_position <= cache.startPosition() );

				if( read_position == cache.startPosition() ){
				// simplest case
				auto buf_size = cache.bufSize();

				if( block.size() == BUF_SIZE && block.get() != nullptr ){// return previous buffer to memory manager
					memory_manager.consumerRelease( block.release( cache.buffer() ), buf_size );
					block.setUsed( buf_size );
				}else{
					block.reset( cache.buffer(), BUF_SIZE, buf_size );
					memory_manager.consumerRelease( nullptr, buf_size );
				}
				cache.moveNext();
				if( read_position == 0 ) cache_sync.unlock(); // do we need unlock?
				read_position += buf_size;
				return buf_size;
			}

			to_read = std::min( (uint64_t)to_read, cache.startPosition() - read_position ); // read up to cache start
		} else
			to_read = std::min( (uint64_t)to_read, write_position - read_position );

		// TODO return previous buffer to memory manager.
		block.ensureSize( to_read ).setUsed( to_read );
		if( to_read > 0 )
			disk->Read( read_position, block.get(), to_read );

		if( read_position == 0 ) cache_sync.unlock(); // do we need unlock?

		read_position += to_read;
		return to_read;
	};


	bool atEnd() const override { return read_position >= write_position; };

	cache_free_status FreeCache( int64_t size_to_free ) override{

		if( read_position > 0 ) // no cache operations when read is started.
			return cache.isEmpty() ? FULL_CLEAN : NO_CLEAN;


		// we try to sync to prevent dead locks because it can be
		// in read mode or currently in writting than we do not clean
		if( cache_sync.try_lock() ){
			if( !cache.isEmpty() ){
				for( ; size_to_free > 0 && !cache.isEmpty(); cache.moveNext(), size_to_free -= BUF_SIZE ){
					DiskWrite( cache.startPosition(), cache.buffer(), cache.bufSize() );
					memory_manager.consumerRelease( cache.buffer(), 0 );
				}

				if( cache.isEmpty() ){
					consumer_idx = -1; // if cache is empty the consumer is unregistered.
					cache.reset(); // reset start counter.
				}

				if( wasClosed ){
					std::lock_guard<std::mutex> lk( file_sync );
					disk->Close();
				}
			}
			cache_sync.unlock();
			return cache.isEmpty()? FULL_CLEAN : PARTIAL_CLEAN;
		}

		return cache.isEmpty()?FULL_CLEAN : NO_CLEAN;
	}

	void Close() override{ if( disk ) { std::lock_guard<std::mutex> lk( file_sync ); disk->Close(); wasClosed = true; } }

	void Remove() override { if( disk ){ std::lock_guard<std::mutex> lk( file_sync ); disk->Remove( true ); disk.reset(); } }

	~BlockCachedFile(){
		if( consumer_idx >= 0 ){
			memory_manager.unregisterConsumer( this, consumer_idx );
			for( ; !cache.isEmpty(); cache.moveNext() ){
				memory_manager.consumerRelease( cache.buffer(), cache.bufSize() );
			}
		}
		Remove(); // remove not used anymore file
	}
private:
	struct CacheStorage {
		const uint32_t buf_size;
		uint32_t last_buf_size;
		std::vector<uint8_t*> bufs;
		std::vector<uint64_t> positions;
		uint32_t start_idx = 0;

		CacheStorage( uint32_t buf_size ) : buf_size(buf_size), last_buf_size(buf_size){}

		inline bool isEmpty() const { return start_idx >= bufs.size(); }
		inline uint64_t startPosition() const { return positions[start_idx]; }
		// current buffer
		inline uint8_t* buffer() const { return bufs[start_idx]; }
		// the size of current buffer
		inline uint32_t bufSize() const { return (start_idx+1) == bufs.size() ? last_buf_size : buf_size; }
		inline uint64_t size() const { return (bufs.size() - start_idx)*buf_size - buf_size + last_buf_size; }
		inline void addToLastBuffer( uint32_t size ) { last_buf_size += size;
																								 assert( last_buf_size <= buf_size ); }
		inline void reset() { start_idx = 0; bufs.clear(); positions.clear(); }

		inline void add( uint8_t* buf, uint64_t position, uint32_t size ){
			assert( last_buf_size == buf_size ); // last only buffer can be other size

			last_buf_size = size;
			bufs.push_back( buf );
			positions.push_back( position );
		}

		void moveNext(){ start_idx++; }

		~CacheStorage(){
			for( uint32_t i = start_idx; i < bufs.size(); i++ )
				delete [] bufs[i]; // this mem is released from memory manager in descrutor of Stream
		}
	};

	fs::path file_name;
	MemoryManager &memory_manager;
	int32_t consumer_idx;
	std::unique_ptr<FileDisk> disk;
	bool wasClosed = false;
	uint64_t write_position = 0;
	uint64_t read_position = 0;
	std::mutex cache_sync;
	std::mutex file_sync;
	CacheStorage cache;
	const uint32_t read_align;


	inline uint32_t DiskWrite( uint64_t pos, uint8_t * buf, uint32_t buf_size ){
		std::lock_guard<std::mutex> lk( file_sync );
		if( !disk ) {
			disk.reset( new FileDisk(file_name) );
			file_name.clear();
		}

		disk->Write( pos, buf, buf_size );
		return buf_size;
	}
};  // end of CachedFileStream


// ============== Sequenc Compacter ========================
struct BlockSequenceCompacterWriter : public IBlockWriter{
	BlockSequenceCompacterWriter( IBlockWriter *disk,
													 const uint16_t &entry_size, const uint8_t &begin_bits )
		: disk( disk ) , entry_size(entry_size), begin_bits_( begin_bits&7 )
		, begin_bytes( begin_bits>>3 ), mask( 0xff00 >> begin_bits_ )
	{
	}

	void Write( StreamBuffer & block ) override {
		CompactBuffer( block.get(), block.used() );
	}

	void Close() override {
		if(disk){
			if( buffer.used() > 0 )
				disk->Write( buffer );
			disk->Close();
			disk.reset();
		}
		buffer.reset();
	}

	~BlockSequenceCompacterWriter(){Close();}
private:
	std::unique_ptr<IBlockWriter> disk;
	uint64_t last_write_value = 0;
	StreamBuffer buffer;
	const uint16_t entry_size;
	const uint8_t begin_bits_;
	const uint8_t begin_bytes;
	const uint8_t mask;

	void CompactBuffer( uint8_t *buf, uint32_t buf_size ){
		assert( (buf_size%entry_size) == 0 );

		auto buffer_idx = buffer.used();
		auto cur_buf = buffer.get();
#define setNext( v ) { assert( buffer_idx < buffer.size() ); cur_buf[buffer_idx++] = v; }

		for( uint32_t i = 0; i < buf_size; ){
//			uint64_t next_val = Util::ExtractNum32( buf + i + begin_bytes, begin_bits_ );
			uint64_t next_val = Util::ExtractNum64( buf + i + begin_bytes, begin_bits_, 32 );
			uint64_t diff = next_val - last_write_value;
			if( diff < 128 ){
				// to save one byte of diff
				setNext( diff );
			} else if( diff < 16512 /* 2^14+128 */ ){
				// to save 2 bytes b10xx xxxx xxxx xxxx: x = diff - 128;
				diff -= 128;
				assert( (diff >> 8) < 64 );
				setNext( 128 + (diff>>8) );
				setNext( diff );
			} else if( diff < 2113664 /* 2^21 + 128 + 2^14 */ ){
				// to save 3 bytes ;
				diff -= 16512;
				setNext( 192 + (diff>>16) );
				setNext( diff>>8 );
				setNext( diff );
			} else if( diff < 270549120 /* 2^28 + 2^21 + 2^14 + 128 */){
				// to save 4 bytes
				diff = (diff-2113664) + (((uint32_t)0xe)<<28);
				((uint32_t*)(cur_buf + buffer_idx))[0] = bswap_32( diff );
				buffer_idx += 4;
			}	else {
				// save 255 & full value
				setNext( 255 );
				((uint32_t*)(cur_buf + buffer_idx))[0] = bswap_32( next_val );
				buffer_idx += 4;
			}

			// Process cuted entry
			for( uint32_t j = 0; j < begin_bytes; j++ )
				setNext( buf[i++] );

			// save one more bytes just in case need to do it i.e. there is a next byte
			if( begin_bits_ || (begin_bytes + 4) < entry_size ){
				setNext( (buf[i]&mask) | (buf[i+4]&~mask) );
				i++;
			}
			i += 4; // move saved bytes forward
			for( uint32_t j = begin_bytes + 5; j < entry_size; j++ )
				setNext( buf[i++] );

			// set result
			last_write_value = next_val;

			if( (buffer_idx + 1 + entry_size ) > BUF_SIZE ){ // next value could not fit buffer
				memset( cur_buf + buffer_idx, 0, BUF_SIZE - buffer_idx ); // remove valgring Error :)
				disk->Write( buffer.setUsed( BUF_SIZE ) );
				cur_buf = buffer.setUsed( buffer_idx = 0 ).get(); // the buffer can be changed after writing.
				last_write_value = 0;
			}
		}

		buffer.setUsed( buffer_idx );
	}

	friend struct BlockSequenceCompacterReader;

};

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
struct BlockSequenceCompacterReader : public IBlockReader {
	BlockSequenceCompacterReader( IBlockReader *disk, const uint16_t &entry_size, uint8_t begin_bits )
		: disk( disk ), entry_size(entry_size), begin_bits_( begin_bits&7 ), begin_bytes( begin_bits>>3 )
		, mask( 0xff00 >> begin_bits_ )
	{	}

	BlockSequenceCompacterReader(  IBlockReader *disk, BlockSequenceCompacterWriter &writer )
		: disk(disk), entry_size(writer.entry_size), begin_bits_( writer.begin_bits_ ), begin_bytes( writer.begin_bytes )
		, mask( writer.mask )
	{
		read_buffer.swap( writer.buffer );
		writer.last_write_value = 0; // it is not necceesary because writer shouldn't be writted after this.
	}

	uint32_t Read( StreamBuffer & block ) override{
		block.setUsed( 0 );
		if( read_buffer.used() == 0 // in case we get buffer from outside
				&& disk->Read( read_buffer ) == 0 )
			return 0; // end of input stream

		return GrowBuffer( read_buffer, block );
	}

	inline bool atEnd() const override { return read_buffer.used() == 0 && disk->atEnd(); }
	void Close() override { if(disk){ disk->Close(); disk.reset(); }read_buffer.reset(); }

	~BlockSequenceCompacterReader(){ if( disk ) disk->Close(); }
private:
	std::unique_ptr<IBlockReader> disk;
	StreamBuffer read_buffer;
	const uint16_t entry_size;
	const uint8_t begin_bits_;
	const uint8_t begin_bytes;
	const uint8_t mask;

	// Returns growed buffer size
	uint32_t GrowBuffer( StreamBuffer &from, StreamBuffer &to ){

		to.ensureSize( from.used()/(entry_size-3)*entry_size ); // if all entries are compacted this is the max buffer size

		uint32_t j = 0;
		auto src = from.get();
		auto dst = to.get();

		for( uint64_t i = 0, last_value = 0, buf_size = from.used();
				 ( i + (buf_size < BUF_SIZE? 0 : entry_size ) + 1 ) <= buf_size;
				 j += entry_size )
			i += RecreateEntry( last_value, src + i, dst + j );

		assert( j < to.size() );

		from.setUsed( 0 ); // clear used for atEnd func
		to.setUsed( j  );
		return j;
	}



	// Recornstruct the entry, update l
	inline uint32_t RecreateEntry( uint64_t &last_value, const uint8_t * src_buf, uint8_t * buf  ) const {
		uint32_t ret = 1, value_to_insert;
		uint8_t compress_code = src_buf[0];

		if( compress_code < 128 )
			value_to_insert = last_value + compress_code;
		else if( compress_code < 192 ){
			value_to_insert = last_value + 128 + ((((uint32_t)(compress_code - 128))<<8) | src_buf[1]);
			ret = 2;
		} else if( compress_code < 224 ){
			value_to_insert = last_value + 16512 + ((((uint32_t)(compress_code-192))<<16) | ( ((uint32_t)src_buf[1]) <<8 ) | src_buf[2] );
			ret = 3;
		} else if( compress_code < 240 ){
			value_to_insert = last_value + 2113664 + ((((uint32_t)(compress_code-224))<<24) |
																								( ((uint32_t)src_buf[1]) <<16) | ( ((uint32_t)src_buf[2]) <<8 ) | src_buf[3] );
			ret = 4;
		} else {
			assert( compress_code == 255 );
			value_to_insert = bswap_32( ((uint32_t*)(src_buf+1))[0] );
			ret = 5;
		}

		bool hasEnd = ( begin_bytes + 4 ) < entry_size || begin_bits_;
		// recreate the entry
		for( uint32_t j = 0; j < begin_bytes; j++ )
			buf[j] = src_buf[ret++];
		auto splitted = hasEnd ? src_buf[ret++] : 0;
		buf[begin_bytes] = (splitted&mask)|(value_to_insert>>(24+begin_bits_));
		buf[begin_bytes+1] = value_to_insert >> (16+begin_bits_);
		buf[begin_bytes+2] = value_to_insert >> (8+begin_bits_);
		buf[begin_bytes+3] = value_to_insert >> begin_bits_;
		if( hasEnd )
			buf[begin_bytes+4] = (value_to_insert<<(8-begin_bits_)) | (splitted&~mask);

		for( uint32_t j = begin_bytes + 5; j < entry_size; j++ )
			buf[j] = src_buf[ret++];

		last_value = value_to_insert;
		return ret;
	}
};



// ============== Bytes Cutter =============================
struct BlockByteCutter : public IBlockWriter{
	BlockByteCutter( IBlockWriter *disk, uint16_t entry_size, uint16_t begin_bits )
		: disk( disk ), entry_size(entry_size)
		, bytes_begin(begin_bits>>3), mask( 0xff00 >> (begin_bits&7) )
	{	}

	void Write( StreamBuffer & block ) override {
		disk->Write( CompactBuffer( block ) );
	}

	void Close() override { if( disk ) { disk->Close(); disk.reset(); } }

	~BlockByteCutter(){ if(disk) disk->Close(); }
private:
	std::unique_ptr<IBlockWriter> disk;
	const uint16_t entry_size;
	const uint8_t bytes_begin;
	const uint8_t mask;

	inline StreamBuffer & CompactBuffer( StreamBuffer & buffer ){
		uint32_t buf_size = buffer.used();
		uint8_t * buf = buffer.get();

		assert( (buf_size%entry_size) == 0 );

		// extract current bucket
		if( mask ){
			for( uint32_t i = bytes_begin + 1; i < buf_size; i += entry_size ){
				buf[i] = (buf[i-1]&mask) | (buf[i]&(~mask));
			}
		}

		// remove current bucket
		uint32_t i = bytes_begin + 1, j = bytes_begin;
		for( ;i < (buf_size - entry_size); i += entry_size, j += entry_size - 1 )
			memmove( buf + j, buf + i, entry_size - 1 );
		// last tail
		memmove( buf + j, buf + i, entry_size - 1 - bytes_begin );

		return buffer.setUsed( buf_size/entry_size * (entry_size-1) );
	}
};

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
struct BlockByteInserter : public IBlockReader{
	BlockByteInserter( IBlockReader *disk, uint16_t entry_size, uint16_t begin_bits, uint8_t removed_byte )
		: disk( disk ), entry_size(entry_size), begin_bits(begin_bits&7), bytes_begin(begin_bits>>3)
		, removed_byte( removed_byte ), mask( 0xff00 >> (begin_bits&7) )
	{	}

	uint32_t Read( StreamBuffer & block ) override{
		assert( disk );
		disk->Read( read_block );
		assert( read_block.used()%(entry_size-1) == 0 );

		return GrowBuffer( read_block, block );
	}

	bool atEnd() const override { return disk->atEnd(); }

	void Close() override { if(disk){ disk->Close(); disk.reset(); } }


	uint32_t GrowBuffer( StreamBuffer & from, StreamBuffer & block ){
		auto buf_size = from.used();
		if(buf_size == 0 ) {
			block.setUsed( 0 );
			return 0;
		}
		uint32_t full_buf_size = buf_size/(entry_size-1)*entry_size;

		auto buf = block.ensureSize( full_buf_size ).setUsed( full_buf_size ).get();
		auto src = from.get();

		if( !mask ){
			for( uint64_t src_idx = 0, dst_idx = 0, insert_idx = bytes_begin; src_idx < buf_size; ){
				// WARNING this writes after buffer size, and should be used just with safe buffers!!!
				((uint64_t*)(buf+dst_idx))[0] = ((uint64_t*)(src+src_idx))[0];
				dst_idx += 8;
				src_idx += 8;
				if( dst_idx > insert_idx ){
					src_idx -= dst_idx - insert_idx;
					buf[insert_idx] = removed_byte;
					dst_idx = insert_idx + 1;
					insert_idx += entry_size;
				}
			}

			return full_buf_size;
		}

		// The restore works from end to start by historical reason
		// once it was copy inside same buffer.
		// may be it worth it to rewrite to be from begin to end

		// copy last entry
		uint64_t tail_idx = full_buf_size - entry_size + bytes_begin + 1;
		int64_t compacted_idx = buf_size - entry_size + bytes_begin + 1;
		memcpy( buf + tail_idx, src + compacted_idx, entry_size - bytes_begin - 1 );
		for( tail_idx -= entry_size, compacted_idx -= entry_size - 1;
				 compacted_idx >= 0;
				 tail_idx -= entry_size, compacted_idx -= entry_size - 1 )
			memcpy( buf + tail_idx, src + compacted_idx, entry_size - 1 );

		memcpy( buf, src, bytes_begin ); // the begining

		// Now insert bucket
		if( mask ){
			const uint8_t shift = begin_bits&7;
			const uint8_t insert_hi = (removed_byte)>>shift;
			const uint8_t insert_lo = (removed_byte)<<(8-shift);
			for( uint32_t i = bytes_begin; i < full_buf_size; i += entry_size ){
				buf[i] = (buf[i+1]&mask) | insert_hi;
				buf[i+1] = (buf[i+1]&~mask) | insert_lo;
			}
		}else{
			for( uint32_t i = bytes_begin; i < full_buf_size; i += entry_size )
				buf[i] = removed_byte;
		}

		return full_buf_size;
	}

	~BlockByteInserter(){ if(disk) disk->Close();}
private:
	std::unique_ptr<IBlockReader> disk;
	const uint16_t entry_size;
	const uint16_t begin_bits;
	const uint8_t bytes_begin;
	const uint8_t removed_byte;
	const uint8_t mask;
	// this one need to improve problem of heap fragmentation
	// in case of cahce it allow less memory heap operations
	StreamBuffer read_block;
};// end of ByteInserterStream



// ============== Not Freeing =============================
struct BlockNotFreeingWriter: IBlockWriter{
	BlockNotFreeingWriter( IBlockWriter * wstream ) :wstream(wstream){}

	void Write( StreamBuffer & block ) override{
		assert( wstream != nullptr );
		wstream->Write( block );
	}
	void Close() override{ if( wstream != nullptr ) { wstream->Close(); wstream = nullptr; } }

	private:
		IBlockWriter * wstream;
};

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
struct BlockNotFreeingReader : IBlockReader {
	BlockNotFreeingReader( IBlockReader * reader ) : rstream(reader) {}

	uint32_t Read( StreamBuffer & block ) override { return rstream->Read(block);};
	virtual bool atEnd() const override {return rstream->atEnd(); };
	virtual void Close() override {};
private:
	IBlockReader * rstream;
};

// ============== Manager Streams =============================
struct BlockThreadSafeReader : IBlockReader{
	// Thread safe means more than one thread reads from underlying stream
	// It means not delete and no close for underlying stream
	BlockThreadSafeReader( IBlockReader *read_stream, std::mutex & sync_mutex )
		: rstream(read_stream), sync_mutex(sync_mutex) {}

	uint32_t Read( StreamBuffer & block ) override {
		std::lock_guard<std::mutex> lk(sync_mutex);
		return rstream->Read( block );
	}
	bool atEnd() const override { return rstream->atEnd();};
	void Close() override{};

private:
	IBlockReader * rstream;
	std::mutex &sync_mutex;
};

struct BlockOneBuffer : IBlockReader {
	BlockOneBuffer( IBlockReader *disk, StreamBuffer &buf )
		: disk(disk), one_buf(buf.size()) { one_buf.swap(buf); };

	uint32_t Read( StreamBuffer & block ) override {
		if( one_buf.used() > 0 ){
			block.swap( one_buf );
			one_buf.reset();
			return block.used();
		}
		return disk->Read( block );
	}
	bool atEnd() const override { return disk->atEnd();};
	void Close() override{};

private:
	std::unique_ptr<IBlockReader> disk;
	StreamBuffer one_buf;
};


// ==================================================================================
// ==================================================================================
// ==================================================================================
// ==================================================================================
struct IWriteDiskStream{
	virtual void Write( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) = 0;
	virtual void Close() = 0;
	virtual ~IWriteDiskStream() = default;
};

struct IReadDiskStream{
	virtual uint32_t Read( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) = 0;
	virtual bool atEnd() const = 0;
	virtual void Close() = 0;
	virtual ~IReadDiskStream() = default;
};

struct IReadWriteStream : public IReadDiskStream, IWriteDiskStream {
	virtual void Remove() = 0;
};


// This class not thread safe!!!
struct AsyncCopyStreamWriter : public IWriteDiskStream {
	AsyncCopyStreamWriter( IWriteDiskStream * disk ) : disk(disk) {}

	void Write( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override {
		if( io_thread ) io_thread->join();
		if( !buffer || last_buf_size < buf_size )
			buffer.reset( Util::NewSafeBuffer( last_buf_size = buf_size ) );
		memcpy( buffer.get(), buf.get(), buf_size );
		io_thread.reset( new std::thread( [this](uint32_t b_size){ disk->Write(buffer, b_size);}, buf_size ) );
	};

	void Close() override {
		if( io_thread ){
			io_thread->join();
			io_thread.reset();
		}
		disk.reset();// may be close it before reset?
		buffer.reset();
	};

	~AsyncCopyStreamWriter(){ if( io_thread )	io_thread->join(); }
private:
	std::unique_ptr<IWriteDiskStream> disk;
	std::unique_ptr<uint8_t[]> buffer;
	std::unique_ptr<std::thread> io_thread;
	uint32_t last_buf_size = 0;
};
// ====================================================
struct AsyncStreamReader : public IReadDiskStream {
	const uint32_t max_buffer_size;

	AsyncStreamReader( IReadDiskStream * disk, uint32_t buf_size )
		: max_buffer_size(buf_size), disk(disk)
	{	}


	uint32_t Read( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override {
		assert( io_thread == nullptr || buf_size == max_buffer_size );
		if( !started ) ReadNext();
		uint32_t res;
		{
			std::lock_guard<std::mutex> lk(sync_mutex);
			// if no thread than read is done.
			if( io_thread == nullptr || !disk ) return 0;

			io_thread->join();
			delete io_thread;

			buf.swap( buffer );
			res = this->buf_size;
		}
		ReadNext();
		return res;
	};

	bool atEnd() const override { return !disk || (buf_size == 0 && disk->atEnd()); };

	void Close() override {
		std::lock_guard<std::mutex> lk(sync_mutex);
		if( io_thread != nullptr ){
			io_thread->join();
			delete io_thread;
			io_thread = nullptr;
		}
		disk.reset();
		buffer.reset();
	};

	~AsyncStreamReader(){ Close(); }
private:
	std::unique_ptr<IReadDiskStream> disk;
	std::unique_ptr<uint8_t[]> buffer;
	uint32_t buf_size = 0;
	std::thread *io_thread = nullptr;
	std::mutex sync_mutex;
	bool started = false;

private:
	inline void ReadNext(){
		// at this point thread should not exists!!!

		std::lock_guard<std::mutex> lk(sync_mutex);
		started = true;

		// Start next read if need or close the resources
		if( !disk || disk->atEnd() ){
			io_thread = nullptr;
			buf_size = 0;
			buffer.reset();
			disk->Close();
			disk.reset();
		} else
			io_thread = new std::thread( [this](){
				if( !buffer ) buffer.reset( new uint8_t[max_buffer_size] );
				buf_size = disk->Read( buffer, max_buffer_size );
			});

	}
};



struct FileStream : public IReadWriteStream {
	FileStream( std::string fileName ) : disk( new FileDisk(fileName) )	{	}

	void Write( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override {
		disk->Write( write_position, buf.get(), buf_size );
		write_position += buf_size;
	};

	uint32_t Read( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override {
		uint32_t to_read = std::min( (uint64_t)buf_size, write_position - read_position );
		disk->Read( read_position, buf.get(), to_read );
		read_position += to_read;
		return to_read;
	};

	bool atEnd() const override { return write_position == read_position; };

	void Close() override { if( disk ) disk->Close(); };
	void Remove() override { if(disk) disk->Remove( true ); }
	~FileStream() { Remove(); }
private:
	std::unique_ptr<FileDisk> disk;
	uint64_t write_position = 0;
	uint64_t read_position = 0;
};

struct WriteFileStream : public IWriteDiskStream{
	WriteFileStream( std::string fileName ) : disk( new FileDisk(fileName) )	{	}
	WriteFileStream( FileDisk * file ) : disk( file ), external_file(true)	{	}

	uint64_t WritePosition() const { return write_position; }

	void Write( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override {
		disk->Write( write_position, buf.get(), buf_size );
		//disk->Flush();
		write_position += buf_size;
	};

	void Close() override {
		if( disk ) {
			disk->Close();
			if( external_file ) disk.release();
			else	disk.reset();
		} };

	~WriteFileStream() { Close(); }
private:
	std::unique_ptr<FileDisk> disk;
	uint64_t write_position = 0;
	const bool external_file = false;
};

struct ReadFileStream : public IReadDiskStream{
	ReadFileStream( const std::string &fileName, uint64_t file_size )
		: file_size(file_size), disk( new FileDisk(fileName, false) )	{	}

	ReadFileStream( FileDisk *file, uint64_t file_size )
		: file_size(file_size), disk(file ), external_file(true)	{	}

	const uint64_t file_size;

	uint32_t Read( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override {
		uint32_t to_read = std::min( (uint64_t)buf_size, file_size - read_position );
		disk->Read( read_position, buf.get(), to_read );
		read_position += to_read;
		return to_read;
	};

	bool atEnd() const override { return file_size == read_position; };

	void Close() override {
		if(disk) {
			if(external_file) disk.release();
			else {
				disk->Close(); disk.reset(); }
		}
	};

	~ReadFileStream() { Close(); }
private:
	std::unique_ptr<FileDisk> disk;
	uint64_t read_position = 0;
	const bool external_file = false;
};


// ============== Sequenc Compacter ========================
struct SequenceCompacterWriter : public IWriteDiskStream{
	SequenceCompacterWriter( IWriteDiskStream *disk,
													 const uint16_t &entry_size, const uint8_t &begin_bits )
		: disk( disk ) , entry_size(entry_size), begin_bits_( begin_bits&7 )
		, begin_bytes( begin_bits>>3 ), mask( 0xff00 >> begin_bits_ )
	{
	}

	void Write( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override {
		if( !buffer )  buffer.reset( new uint8_t[BUF_SIZE] ); // init buffer on first write
		CompactBuffer( buf.get(), buf_size );
	}

	uint32_t ReleaseBuffer( uint8_t* buf ){
		if( !disk || !buffer_idx ) return 0;
		memcpy( buf, buffer.get(), buffer_idx );
		auto res = buffer_idx;
		buffer_idx = 0;
		return res;
	}

	void Close() override {
		if(disk){
			if( buffer_idx > 0 )
				disk->Write( buffer, buffer_idx );
			disk->Close();
			disk.reset();
		}
		buffer.reset();
	}

	~SequenceCompacterWriter(){Close();}
private:
	std::unique_ptr<IWriteDiskStream> disk;
	uint64_t last_write_value = 0;
	uint32_t buffer_idx = 0;
	std::unique_ptr<uint8_t[]> buffer;
	const uint16_t entry_size;
	const uint8_t begin_bits_;
	const uint8_t begin_bytes;
	const uint8_t mask;

	void CompactBuffer( uint8_t *buf, uint32_t buf_size ){
		assert( (buf_size%entry_size) == 0 );

		auto cur_buf = buffer.get();
#undef setNext
#define setNext( v ) cur_buf[buffer_idx++] = v

		for( uint32_t i = 0; i < buf_size; ){
//			uint64_t next_val = Util::ExtractNum32( buf + i + begin_bytes, begin_bits_ );
			uint64_t next_val = Util::ExtractNum64( buf + i + begin_bytes, begin_bits_, 32 );
			uint64_t diff = next_val - last_write_value;
			if( diff < 128 ){
				// to save one byte of diff
				setNext( diff );
			} else if( diff < 16512 /* 2^14+128 */ ){
				// to save 2 bytes b10xx xxxx xxxx xxxx: x = diff - 128;
				diff -= 128;
				assert( (diff >> 8) < 64 );
				setNext( 128 + (diff>>8) );
				setNext( diff );
			} else if( diff < 2113664 /* 2^21 + 128 + 2^14 */ ){
				// to save 3 bytes ;
				diff -= 16512;
				setNext( 192 + (diff>>16) );
				setNext( diff>>8 );
				setNext( diff );
			} else if( diff < 270549120 /* 2^28 + 2^21 + 2^14 + 128 */){
				// to save 4 bytes
				diff = (diff-2113664) + (((uint32_t)0xe)<<28);
				((uint32_t*)(cur_buf + buffer_idx))[0] = bswap_32( diff );
				buffer_idx += 4;
			}	else {
				// save 255 & full value
				setNext( 255 );
				((uint32_t*)(cur_buf + buffer_idx))[0] = bswap_32( next_val );
				buffer_idx += 4;
			}

			// Process cuted entry
			for( uint32_t j = 0; j < begin_bytes; j++ )
				setNext( buf[i++] );

			// save one more bytes just in case need to do it i.e. there is a next byte
			if( begin_bits_ || (begin_bytes + 4) < entry_size ){
				setNext( (buf[i]&mask) | (buf[i+4]&~mask) );
				i++;
			}
			i += 4; // move saved bytes forward
			for( uint32_t j = begin_bytes + 5; j < entry_size; j++ )
				setNext( buf[i++] );

			// set result
			last_write_value = next_val;

			if( (buffer_idx + 1 + entry_size ) > BUF_SIZE ){ // next value could not fit buffer
				memset( cur_buf + buffer_idx, 0, BUF_SIZE - buffer_idx ); // remove valgring Error :)
				disk->Write( buffer, BUF_SIZE );
				cur_buf = buffer.get(); // the buffer can be changed after writing.
				buffer_idx = 0;
				last_write_value = 0;
			}
		}
	}
};

struct SequenceCompacterReader : public IReadDiskStream{
	SequenceCompacterReader( IReadDiskStream *disk, const uint16_t &entry_size, uint8_t begin_bits, SequenceCompacterWriter * writer = nullptr )
		: disk( disk ), buffer( Util::NewSafeBuffer(BUF_SIZE) )
		, entry_size(entry_size), begin_bits_( begin_bits&7 ), begin_bytes( begin_bits>>3 )
		, mask( 0xff00 >> begin_bits_ )
	{
		if( writer != nullptr )
			bytes_in_buffer = writer->ReleaseBuffer( buffer.get() );
	}

	uint32_t Read( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override{
		return GrowBuffer( buf.get(), buf_size );
	}

	inline bool atEnd() const override { return bytes_in_buffer <= processed_of_buffer && disk->atEnd(); }
	void Close() override { if(disk){ disk->Close(); disk.reset(); }buffer.reset(); }

	~SequenceCompacterReader(){if(disk) disk->Close();}
private:
	std::unique_ptr<IReadDiskStream> disk;
	uint64_t last_value = 0;
	std::unique_ptr<uint8_t[]> buffer;
	uint32_t processed_of_buffer = 0;
	uint32_t bytes_in_buffer = 0;
	const uint16_t entry_size;
	const uint8_t begin_bits_;
	const uint8_t begin_bytes;
	const uint8_t mask;

	// Returns growed buffer size
	uint32_t GrowBuffer( uint8_t * buf, const uint32_t &buf_size ){
		assert( buf_size%entry_size  == 0 ); // buffer size should be aligned to entry size

		uint32_t i = 0;
		while( i < buf_size && !atEnd() ){

			if( ( processed_of_buffer + entry_size + 1 ) > bytes_in_buffer && !disk->atEnd() ){
				// current buffer finished or first hasn't started yet => need to read next
				bytes_in_buffer = disk->Read( buffer, BUF_SIZE ); // read the stream

				assert( bytes_in_buffer == BUF_SIZE || disk->atEnd() );

				processed_of_buffer = 0;
				last_value = 0;
			} else {
				processed_of_buffer += RecreateEntry( last_value, buffer.get() + processed_of_buffer, buf + i );
				i += entry_size;
			}
		}

		return i;
	}
	


	// Recornstruct the entry, update l
	inline uint32_t RecreateEntry( uint64_t &last_value, const uint8_t * src_buf, uint8_t * buf  ) const {
		uint32_t ret = 1, value_to_insert;
		uint8_t compress_code = src_buf[0];

		if( compress_code < 128 )
			value_to_insert = last_value + compress_code;
		else if( compress_code < 192 ){
			value_to_insert = last_value + 128 + ((((uint32_t)(compress_code - 128))<<8) | src_buf[1]);
			ret = 2;
		} else if( compress_code < 224 ){
			value_to_insert = last_value + 16512 + ((((uint32_t)(compress_code-192))<<16) | ( ((uint32_t)src_buf[1]) <<8 ) | src_buf[2] );
			ret = 3;
		} else if( compress_code < 240 ){
			value_to_insert = last_value + 2113664 + ((((uint32_t)(compress_code-224))<<24) | 
																								( ((uint32_t)src_buf[1]) <<16) | ( ((uint32_t)src_buf[2]) <<8 ) | src_buf[3] );
			ret = 4;
		} else {
			assert( compress_code == 255 );
			value_to_insert = bswap_32( ((uint32_t*)(src_buf+1))[0] );
			ret = 5;
		}

		bool hasEnd = ( begin_bytes + 4 ) < entry_size || begin_bits_;
		// recreate the entry
		for( uint32_t j = 0; j < begin_bytes; j++ )
			buf[j] = src_buf[ret++];
		auto splitted = hasEnd ? src_buf[ret++] : 0;
		buf[begin_bytes] = (splitted&mask)|(value_to_insert>>(24+begin_bits_));
		buf[begin_bytes+1] = value_to_insert >> (16+begin_bits_);
		buf[begin_bytes+2] = value_to_insert >> (8+begin_bits_);
		buf[begin_bytes+3] = value_to_insert >> begin_bits_;
		if( hasEnd )
			buf[begin_bytes+4] = (value_to_insert<<(8-begin_bits_)) | (splitted&~mask);

		for( uint32_t j = begin_bytes + 5; j < entry_size; j++ )
			buf[j] = src_buf[ret++];

		last_value = value_to_insert;
		return ret;
	}
};


struct BufferedWriter : IWriteDiskStream{
	// all written out buffer except last would be in size of block_size, but amount of real data inside could be less
	BufferedWriter( IWriteDiskStream* disk, const uint32_t block_size = BUF_SIZE, const uint8_t entry_size = 1 )
		: disk(disk), allocated_size( block_size ), write_size( block_size/entry_size*entry_size )
	{	}

	void Write( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ){
		Write( buf.get(), buf_size );
	};

	inline void Write( const uint8_t * buf, const uint32_t &buf_size ) {
		if( !buffer ) buffer.reset( new uint8_t[allocated_size] ); // late buffer allocation.

		for( uint32_t processed = 0; processed < buf_size; ){
			uint32_t to_process = std::min( write_size - buffer_used, buf_size - processed );
			memcpy( buffer.get() + buffer_used, buf + processed, to_process );
			processed += to_process;
			buffer_used += to_process;
			if( buffer_used == write_size ){
				disk->Write( buffer, allocated_size );
				buffer_used = 0;
			}
		}
	}

	void Close(){
		if( disk ) {
			if( buffer_used ) disk->Write( buffer, buffer_used );
			buffer_used = 0;
			disk.reset(); // this should desctruct i.e. also close it
		}
		if( buffer ) buffer.reset();
	};

	uint32_t ReleaseBuffer( std::unique_ptr<uint8_t[]> &buf ){
		buffer.swap( buf );
		auto res = buffer_used;
		buffer_used = 0;
		return res;
	}

	~BufferedWriter(){ Close(); }

private:
	std::unique_ptr<IWriteDiskStream> disk;
	const uint32_t allocated_size;
	const uint32_t write_size;
	uint32_t buffer_used = 0;
	std::unique_ptr<uint8_t[]> buffer;

	friend struct BufferedReadStream;
};

struct BufferedReadStream : IReadDiskStream {

	BufferedReadStream( IReadDiskStream * disk, BufferedWriter &writer )
		:disk(disk), block_size(writer.allocated_size), read_size(writer.write_size)
		, buffer_size( writer.buffer_used )
		, buffer( buffer_size == 0 ? Util::NewSafeBuffer( block_size ) : writer.buffer.release() )
	{
		writer.buffer_used = 0;
	}

	BufferedReadStream( IReadDiskStream * disk, const uint32_t block_size = BUF_SIZE, const uint8_t entry_size = 1 )
		: disk(disk), block_size(block_size), read_size( block_size/entry_size*entry_size )
		, buffer( Util::NewSafeBuffer(block_size) )
	{	}

	uint32_t Read( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override{
		uint32_t res = 0;

		while( res < buf_size ){
			if( buffer_position >= buffer_size ){
				buffer_size = std::min( read_size, disk->Read( buffer, block_size ) ); // hope that read is OK and amount allinged to entry size
				buffer_position = 0;
			}
			if( buffer_size == 0 ) return res;

			uint32_t to_read = std::min( buffer_size-buffer_position, buf_size - res );
			memcpy( buf.get() + res, buffer.get() + buffer_position, to_read );
			buffer_position += to_read;
			res+=to_read;
		}

		return res;
	}
	bool atEnd() const override { return  disk->atEnd(); }
	void Close() override { disk->Close(); };

private:
	std::unique_ptr<IReadDiskStream> disk;
	const uint32_t block_size;
	const uint32_t read_size;
	uint32_t buffer_size = 0;
	uint32_t buffer_position = 0;
	std::unique_ptr<uint8_t[]> buffer;
};



// creates or cached or simple file stream depends of settings of memory manager
IBlockWriterReader * CreateFileStream( const fs::path &fileName, MemoryManager &memory_manager, uint32_t read_align_bytes = 1 ){
	return memory_manager.CacheEnabled ?
				((IBlockWriterReader*)(new BlockCachedFile( fileName, memory_manager, read_align_bytes )))
			: new BlockedFileStream( fileName, read_align_bytes );
}

// =======================================================
// This stream not garantee same read order as write order
struct BucketStream{
	BucketStream( std::string fileName, MemoryManager &memory_manager, uint16_t bucket_no, uint8_t log_num_buckets,
								uint16_t entry_size, uint16_t begin_bits, bool compact = true, int8_t sequence_start_bit = -1 )
		: fileName(fileName)
		, memory_manager( memory_manager )
		, bucket_no_( bucket_no )
		, log_num_buckets_( log_num_buckets )
		, entry_size_(entry_size)
		, begin_bits_(begin_bits)
		, sequence_start_bit(sequence_start_bit)
		, compact( compact&(log_num_buckets_>7) )
	{	}

	const std::string getFileName() const {return fileName;}


	// The write is NOT thread save and should be synced from outside!!!!
	inline void Write( StreamBuffer & buf ){
		assert( (buf.size() % entry_size_) == 0 );
		if( buf.size() == 0 ) return; // nothing to write

		// std::lock_guard<std::mutex> lk( sync_mutex );
		if( !disk_output ){
			uint32_t read_align_size = entry_size_;
			if( compact )
				read_align_size = sequence_start_bit >= 0 ? 1 : (entry_size_-1);

			bucket_file.reset( CreateFileStream( fileName, memory_manager, read_align_size ) );
			disk_output.reset( new BlockNotFreeingWriter( bucket_file.get() ) );
			if( compact ){
				if( sequence_start_bit >= 0 )
					disk_output.reset( compacter = new BlockSequenceCompacterWriter( disk_output.release(), entry_size_ - 1
																, sequence_start_bit - (sequence_start_bit > begin_bits_ ? 8 : 0 ) ) );
				else
					disk_output.reset( buf_writer = new BlockBufferedWriter( disk_output.release(), read_align_size ) );

				disk_output.reset( new BlockByteCutter( disk_output.release(), entry_size_, begin_bits_ - log_num_buckets_ ) );
			} else
				disk_output.reset( buf_writer = new BlockBufferedWriter( disk_output.release(), read_align_size ) );
		}

		disk_output->Write( buf );
	}


	// This is NOT thread safe
	void FlusToDisk(){
		if( disk_output ) {
			disk_output->Close();
			disk_output.reset();
		}
		compacter = nullptr; // this cleared by resetting disk_output
		buf_writer = nullptr; // this cleared by resetting disk_output
		if( bucket_file )
			((IBlockWriter*)bucket_file.get())->Close();
	}

	// This is thread safe
	IBlockReader * CreateReader( bool isThreadSafe = true ){
		std::lock_guard<std::mutex> lk( sync_mutex );
		if( !bucket_file ) return nullptr; // nothing to create...

		IBlockReader* fin = bucket_file.get();
		if( isThreadSafe )
			fin = new BlockThreadSafeReader( fin, sync_mutex );
		else
			fin = new BlockNotFreeingReader( fin );

		if( compact ){
			if( sequence_start_bit >= 0 ){
				fin = (disk_output&&compacter) ? new BlockSequenceCompacterReader( fin, *compacter )
						: new BlockSequenceCompacterReader( fin, entry_size_ - 1,
											sequence_start_bit - (sequence_start_bit > begin_bits_ ? 8 : 0 ) );
				compacter = nullptr;
			} else if( buf_writer ){
				fin = new BlockOneBuffer( fin, buf_writer->buffer );
				buf_writer = nullptr;
			}

			fin = new BlockByteInserter( fin, entry_size_,
										begin_bits_ - log_num_buckets_, bucket_no_ >> (log_num_buckets_-8) );
		} else if( buf_writer ){
			fin = new BlockOneBuffer( fin, buf_writer->buffer );
			buf_writer = nullptr;
		}

		disk_output.reset();

		return fin;
	}

	// This is NOT thread safe
	uint32_t Read( StreamBuffer &buf ){
		if( !disk_input ) disk_input.reset( CreateReader( false ) );
		return disk_input->Read( buf );
	}

	// It can be called if full reading not finished yet to close and remove resources
	void EndToRead(){
		std::lock_guard<std::mutex> lk( sync_mutex );
		if( disk_input ){
			disk_input->Close();
			disk_input.reset();
		}
		if( bucket_file ) {
			bucket_file.get()->Remove();
			bucket_file.reset();
		}
	}

private:
	const std::string fileName;
	MemoryManager &memory_manager;
	std::unique_ptr<IBlockWriterReader> bucket_file;
	std::unique_ptr<IBlockWriter> disk_output;
	std::unique_ptr<IBlockReader> disk_input;
	const uint16_t bucket_no_;
	const uint8_t log_num_buckets_;
	std::unique_ptr<std::thread> disk_io_thread;
	const uint16_t entry_size_;
	const uint16_t begin_bits_;
	const int8_t sequence_start_bit;
	std::mutex sync_mutex;

	const bool compact;

	BlockSequenceCompacterWriter * compacter = nullptr;
	BlockBufferedWriter * buf_writer = nullptr;
};

struct LastTableRewrited : IReadDiskStream {
	LastTableRewrited( IReadDiskStream *disk, const uint8_t k,const uint16_t entry_size,
											const uint64_t num_entries, const fs::path &filename )
		: k(k), f7_shift( 128 - k ), pos_offset_size( k + kOffsetSize )
		, t7_pos_offset_shift( f7_shift - pos_offset_size ), entry_size(entry_size)
		, num_entries(num_entries), filename_(filename), disk(disk)
	{}

	uint32_t Read( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override {
		if( !bitfield_ ){
			bitfield_.reset( new bitfield( num_entries, filename_ ) );
			index.reset( new bitfield_index( *bitfield_.get() ) );
		}

		auto read = disk->Read( buf, buf_size );

		// clear ram as soon as possible
		if( read == 0 ){
			index.reset();
			bitfield_.reset();
		}

		for( uint64_t buf_ptr = 0; buf_ptr < read; buf_ptr += entry_size ){

			uint8_t * entry = buf.get() + buf_ptr;

			// table 7 is special, we never drop anything, so just build
			// next_bitfield
			uint64_t entry_f7 = Util::SliceInt64FromBytes(entry, 0, k);
			uint64_t entry_pos_offset = Util::SliceInt64FromBytes(entry, k, pos_offset_size);

			uint64_t entry_pos = entry_pos_offset >> kOffsetSize;
			uint64_t entry_offset = entry_pos_offset & ((1U << kOffsetSize) - 1);

			// assemble the new entry and write it to the sort manager

			// map the pos and offset to the new, compacted, positions and
			// offsets
			std::tie(entry_pos, entry_offset) = index->lookup(entry_pos, entry_offset);
			entry_pos_offset = (entry_pos << kOffsetSize) | entry_offset;

			// table 7 is already sorted by pos, so we just rewrite the
			// pos and offset in-place
			uint8_t tmp[16];
			uint128_t new_entry = (uint128_t)entry_f7 << f7_shift;
			new_entry |= (uint128_t)entry_pos_offset << t7_pos_offset_shift;
			Util::IntTo16Bytes( tmp, new_entry );
			memcpy( entry, tmp, entry_size );
		}
		return read;
	};

	bool atEnd() const override { return disk->atEnd(); }
	void Close() override {
		bitfield_.reset();
		index.reset();
		if(disk ) disk->Close();
		disk.reset();
	};

private:
	const uint8_t k;
	uint8_t const f7_shift;
	uint8_t const pos_offset_size;
	uint8_t const t7_pos_offset_shift;
	uint16_t const entry_size;
	const uint64_t num_entries;
	fs::path filename_;

	std::unique_ptr<IReadDiskStream> disk;
	std::unique_ptr<bitfield> bitfield_;
	std::unique_ptr<bitfield_index> index;
};

IReadDiskStream * CreateLastTableReader( FileDisk * file, uint8_t k, uint16_t entry_size,
																				bool withCompaction, uint32_t max_buffer_size = 0 ){
	IReadDiskStream *res = new ReadFileStream( file, file->GetWriteMax() );

	if( withCompaction && k >= 30 )
		res = new SequenceCompacterReader( res, entry_size, k );


	if( max_buffer_size > 0 )
		res = new AsyncStreamReader( res, max_buffer_size );

	return res;
}

IReadDiskStream * CreateLastTableReader( FileDisk * file, uint8_t k, uint16_t entry_size,
																				 uint64_t num_entries, const fs::path &bitfield_filename,
																				 bool withCompaction, uint32_t max_buffer_size = 0 ){
	IReadDiskStream *res = new ReadFileStream( file, file->GetWriteMax() );

	if( withCompaction && k >= 30 )
		res = new SequenceCompacterReader( res, entry_size, k );

	res = new LastTableRewrited( res, k, entry_size, num_entries, bitfield_filename );

	if( max_buffer_size > 0 )
		res = new AsyncStreamReader( res, max_buffer_size );

	return res;
}




struct LastTableScanner : IWriteDiskStream
{
	LastTableScanner( IWriteDiskStream * disk, uint8_t k, uint16_t entry_size, bool full_scan,
										const char * bitfield_filename )
		: disk( disk ), bitmap( full_scan ? new bitfield( ((uint64_t)2)<<k /* POSSIBLE SIZE PROBLEM */ ) : nullptr )
		, filename( bitfield_filename ), k(k), entry_size( entry_size )
		, pos_offset_size( k + kOffsetSize )
	{}

	void Write( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ){
		disk->Write( buf, buf_size );

		table_size += buf_size / entry_size;

		for( int64_t buf_ptr = 0; buf_ptr < buf_size; buf_ptr += entry_size ){
			int64_t entry_pos_offset = Util::SliceInt64FromBytes( buf.get() + buf_ptr, k, pos_offset_size);
			uint64_t entry_pos = entry_pos_offset >> kOffsetSize;
			uint64_t entry_offset = entry_pos_offset & ((1U << kOffsetSize) - 1);
			if( bitmap ){
				bitmap->set( entry_pos );
				bitmap->set( entry_pos + entry_offset );
			}
			else{
				// if( entry_pos > max_entry_index ) max_entry_index = entry_pos;
				if( ( entry_pos + entry_offset ) > max_entry_index ) max_entry_index = entry_pos + entry_offset;
			}
		}
	}

	void Close() {
		disk->Close();
		disk.reset();

		if( bitmap ){
			auto cres = bitmap->MoveToTable7();
			std::cout << "Compacting bitfield of last table: " << ( cres ? " COMPACTED " : " FAILED!!! " ) << std::endl;
			// if( !cres ) bitmap.resize( table_size );
		}
		else
			bitmap.reset( new bitfield( table_size, max_entry_index ) );

		bitmap->FlushToDisk( filename.c_str() );
	}

	~LastTableScanner() { if(disk) Close(); }
private:
	std::unique_ptr<IWriteDiskStream> disk;
	std::unique_ptr<bitfield> bitmap;
	std::string filename;
	const uint8_t k;
	const uint16_t entry_size;
	uint8_t const pos_offset_size;
	uint64_t max_entry_index = 0;
	uint64_t table_size = 0;
};

/* scan_mode: 0 - no scan, 1 - full scan, 2 - quick scan */
IWriteDiskStream * CreateLastTableWriter( FileDisk * file, uint8_t k, uint16_t entry_size,
																					bool withCompaction, bool withScan, bool full_scan ){
	IWriteDiskStream * res = new WriteFileStream( file );

	if( withCompaction && k >= 30 )
		res = new SequenceCompacterWriter( res, entry_size, k );

	if( withScan )
		res = new LastTableScanner( res, k, entry_size, full_scan,
								( file->GetFileName() + ".bitfield.tmp").c_str() );

	res = new AsyncCopyStreamWriter( res );

	return res;
}

struct ReadStreamToDisk : Disk {
	ReadStreamToDisk( IReadDiskStream *strm, uint16_t entry_size)
		: allocated_size( BUF_SIZE/entry_size*entry_size )
		, buffer(new uint8_t[allocated_size] )
		, read_stream( new AsyncStreamReader( strm, allocated_size ) )
	{	}

	uint8_t const* Read(uint64_t begin, uint64_t length) override{
		assert( length <= allocated_size );
		assert( buffer && read_stream );

		while( begin >= (buffer_begin + buf_size) ){
			buffer_begin += buf_size;
			buf_size = read_stream->Read( buffer, allocated_size );
		}

		assert( (begin+length) <= (buffer_begin + buf_size) );
		return buffer.get() + ( begin - buffer_begin );
	};
	void Write(uint64_t begin, const uint8_t *memcache, uint64_t length) override {
		throw InvalidStateException( "Write impossible to read only disk" );
	};

	void Truncate(uint64_t new_size) override {
		throw InvalidStateException( "Truncate impossible to read only disk" );
	};
	std::string GetFileName() override { return "read_disk"; };
	void FreeMemory() override { buffer.reset(); read_stream.reset(); };

private:
	const uint32_t allocated_size;
	uint64_t buffer_begin = 0;
	uint64_t buf_size = 0;
	std::unique_ptr<uint8_t[]> buffer;
	std::unique_ptr<IReadDiskStream> read_stream;
};



#endif  // SRC_CPP_DISK_STREAMS_HPP_
