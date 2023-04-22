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
struct AsyncStreamWriter : public IWriteDiskStream {
	AsyncStreamWriter( IWriteDiskStream * disk, const uint32_t max_buffer_size = 0 )
		: disk(disk), last_buf_size(max_buffer_size) {}

	void Write( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override {
		if( io_thread ) io_thread->join();
		if( !buffer || last_buf_size < buf_size )
			buffer.reset( new uint8_t[last_buf_size = buf_size] );
		buf.swap( buffer );
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

	~AsyncStreamWriter(){ if( io_thread ) io_thread->join(); }
private:
	std::unique_ptr<IWriteDiskStream> disk;
	std::unique_ptr<uint8_t[]> buffer;
	std::unique_ptr<std::thread> io_thread;
	uint32_t last_buf_size = 0;
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
	const uint16_t bufs_count;

	AsyncStreamReader( IReadDiskStream * disk, uint32_t buf_size, uint16_t bufs_count = 10 )
		: max_buffer_size(buf_size), bufs_count(bufs_count), disk(disk)
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
	SequenceCompacterReader( IReadDiskStream *disk, const uint16_t &entry_size, uint8_t begin_bits )
		: disk( disk )
		, entry_size(entry_size), begin_bits_( begin_bits&7 ), begin_bytes( begin_bits>>3 )
		, mask( 0xff00 >> begin_bits_ )
	{
	}

	// this function for support custom buffer growing, defines buffer size to be provided to restore all entries
	static inline uint32_t getMaxRestoredBufferSize( uint32_t entry_size )
	{
		return BUF_SIZE/(entry_size-3)*entry_size + entry_size;
	}

	uint32_t Read( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override{
		return GrowBuffer( buf.get(), buf_size );
	}

	inline bool atEnd() const override { return bytes_in_buffer <= processed_of_buffer && disk->atEnd(); }
	void Close() override { if(disk){ disk->Close(); disk.reset(); }buffer.reset(); }

	// Returns growed buffer size
	uint32_t GrowCustomBuffer( const uint8_t * compacted_buf, uint32_t compacted_buf_size, uint8_t * restore_buf, uint32_t restore_buf_size ) const {

		uint32_t i = 0, processed_of_buffer = 0;
		uint64_t last_value = 0;

		assert( compacted_buf_size <= BUF_SIZE );
		assert( restore_buf_size >= getMaxRestoredBufferSize( entry_size ) );

		if( compacted_buf_size >= BUF_SIZE  ) // in reality it could be equal or less
			compacted_buf_size = compacted_buf_size - entry_size; // it can be new entry started after this point
		restore_buf_size = restore_buf_size/entry_size*entry_size; // align restore buffer to entry size?


		while( i < restore_buf_size && processed_of_buffer < compacted_buf_size ){

			processed_of_buffer += RecreateEntry( last_value, compacted_buf + processed_of_buffer, restore_buf + i );
			i += entry_size;
		}

		assert( processed_of_buffer >= compacted_buf_size ); // check whole buffer is recreated

		return i;
	}

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
				if( !buffer ) buffer.reset( new uint8_t[BUF_SIZE] ); // create buffer if not exists
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


struct ByteCutterStream : public IWriteDiskStream{
	ByteCutterStream( IWriteDiskStream *disk, uint16_t entry_size, uint16_t begin_bits )
		: disk( disk ), entry_size(entry_size)
		, bytes_begin(begin_bits>>3), mask( 0xff00 >> (begin_bits&7) )
	{	}

	void Write( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override{
		disk->Write( buf, CompactBuffer(buf.get(), buf_size) );
	}

	void Close() override { if( disk ) { disk->Close(); disk.reset(); } }

	~ByteCutterStream(){ if(disk) disk->Close(); }
private:
	std::unique_ptr<IWriteDiskStream> disk;
	const uint16_t entry_size;
	const uint8_t bytes_begin;
	const uint8_t mask;

	inline uint32_t CompactBuffer( uint8_t *buf, uint32_t buf_size ){
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

		return buf_size/entry_size * (entry_size-1);
	}
};

struct ByteInserterStream : public IReadDiskStream{
	ByteInserterStream( IReadDiskStream *disk, uint16_t entry_size, uint16_t begin_bits, uint8_t removed_byte )
		: disk( disk ), entry_size(entry_size), begin_bits(begin_bits&7), bytes_begin(begin_bits>>3)
		, removed_byte( removed_byte ), mask( 0xff00 >> (begin_bits&7) )
	{	}

	uint32_t Read( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override{
		uint32_t read = disk->Read( buf, buf_size/entry_size*(entry_size-1) );
		assert( read%(entry_size-1) == 0 );

		return GrowBuffer(buf.get(), read);
	}

	bool atEnd() const override { return disk->atEnd(); }

	void Close() override { if(disk){ disk->Close(); disk.reset(); } }


	uint32_t GrowBuffer( uint8_t * buf, uint32_t buf_size ){
		if(buf_size == 0 ) return 0;
		uint32_t full_buf_size = buf_size/(entry_size-1)*entry_size;

		// copy last entry
		uint64_t tail_idx = full_buf_size - entry_size + bytes_begin + 1;
		int64_t compacted_idx = buf_size - entry_size + bytes_begin + 1;
		memmove( buf + tail_idx, buf + compacted_idx, entry_size - bytes_begin - 1 );
		for( tail_idx -= entry_size, compacted_idx -= entry_size - 1;
				 compacted_idx >= 0;
				 tail_idx -= entry_size, compacted_idx -= entry_size - 1 )
			memmove( buf + tail_idx, buf + compacted_idx, entry_size - 1 );

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

	~ByteInserterStream(){ if(disk) disk->Close();}
private:
	std::unique_ptr<IReadDiskStream> disk;
	const uint16_t entry_size;
	const uint16_t begin_bits;
	const uint8_t bytes_begin;
	const uint8_t removed_byte;
	const uint8_t mask;


};


struct BufferedWriter : IWriteDiskStream{
	BufferedWriter( IWriteDiskStream* disk, const uint32_t &entry_size )
		: disk(disk), allocated_size( BUF_SIZE/entry_size*entry_size )
	{}

	void Write( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ){
		Write( buf.get(), buf_size );
	};

	inline void Write( const uint8_t * buf, const uint32_t &buf_size ){
		if( !buffer ) buffer.reset( new uint8_t[allocated_size] ); // late buffer allocation.

		for( uint32_t processed = 0; processed < buf_size; ){
			uint32_t to_process = std::min( allocated_size - buffer_used, buf_size - processed );
			memcpy( buffer.get() + buffer_used, buf + processed, to_process );
			processed += to_process;
			buffer_used += to_process;
			if( buffer_used == allocated_size ){
				disk->Write( buffer, buffer_used );
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
	uint32_t buffer_used = 0;
	std::unique_ptr<uint8_t[]> buffer;
};

struct NotFreeingWriteStream: IWriteDiskStream{
	NotFreeingWriteStream( IWriteDiskStream * wstream ) :wstream(wstream){}

	void Write( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override{
		wstream->Write( buf, buf_size );
	}
	void Close() override{ }

	private:
		IWriteDiskStream * wstream;
};


/* This class stores in cache as many data as possible before wrtting it to file */
struct CachedFileStream : IReadWriteStream, ICacheConsumer {

	// @max_buf_size - is a size of buffer supposed to be full. if such side provided it is not copied but replaced
	CachedFileStream( const fs::path &fileName, MemoryManager &memory_manager, uint32_t base_buf_size )
		: memory_manager(memory_manager), buf_size(base_buf_size)
	{
		consumer_idx = memory_manager.registerConsumer( this );
		file_name = fileName;
	}

	void Write( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override{
		std::lock_guard<std::mutex> lk(file_sync);

		if( consumer_idx != nullptr // we are caching
				&& memory_manager.request( buf_size, this ) // there is ram to store
				) {
			if( buf_size == this->buf_size ){
				// buffer is allinged write by getting it
				addToCache( buf.release(), buf_size );
				buf.reset( new uint8_t[buf_size] );
			}
			else
			{ // buffer is not alligned write by copy
				uint8_t * new_buf = new uint8_t[buf_size];
				memcpy( new_buf, buf.get(), buf_size );
				addToCache( new_buf, buf_size );
			}
		} else
			DiskWrite( write_position, buf.get(), buf_size );

		write_position += buf_size;
	}

	// warning this is not thread safe function
	uint32_t Read( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override{
		if( consumer_idx != nullptr ){
			// to prevent futher locks on read we unregister started to read stream
			std::lock_guard<std::mutex> lk(file_sync);
			if( consumer_idx != nullptr ){
				memory_manager.unregisterConsumer( consumer_idx );
				consumer_idx = nullptr;
			}
		}

		if( cache != nullptr ){
			// we have cached buffers
			if( read_position == cache->start_pos && buf_size == cache->buf_size ){
				// simplest case
				buf.reset( cache->buf );
				cache->buf = nullptr; // clear without deletion
				moveNextCache();
				read_position += buf_size;
				return buf_size;
			}
		}

		return Read( buf.get(), buf_size );
	};


	bool atEnd() const override{
		return read_position >= write_position;
	};

	uint64_t getUsedCache() const override {
		uint64_t res = 0;

		for( auto c = cache; c != nullptr; c = c->next )
			res += c->buf_size;

		return res;
	}
	void DetachFromCache() override {
		assert( consumer_idx == nullptr || consumer_idx->consumer == this );
		consumer_idx = nullptr;
	}
	void FreeCache() override{
		std::lock_guard<std::mutex> lk(file_sync);
		assert( consumer_idx == nullptr || consumer_idx->consumer == this );

		consumer_idx = nullptr; // this call clears from consumers
		while( cache != nullptr ){
			DiskWrite( cache->start_pos, cache->buf, cache->buf_size );
			moveNextCache();
		}
	}

	void Close() override{ if( disk ) disk->Close(); }

	void Remove() override { if(disk){ disk->Remove( true ); disk.reset(); } }

	~CachedFileStream(){
		if( consumer_idx != nullptr ){
			memory_manager.unregisterConsumer( consumer_idx );
			memory_manager.release( getUsedCache(), this );
		}
		Remove();
		while( cache != nullptr )
			moveNextCache();
	}
private:
	struct CacheEntry{
		uint8_t * buf;
		uint32_t buf_size;
		uint64_t start_pos;

		CacheEntry * next = nullptr;

		CacheEntry( uint8_t * buf, uint32_t size, uint64_t pos ) : buf(buf), buf_size(size), start_pos(pos){}
		void clear(){
			if( buf != nullptr ) {
				delete[]buf;
				buf = nullptr;
			}
		}

		uint32_t leftSize( uint64_t read_pos ) const {
			return (start_pos + buf_size) - read_pos;
		}
		uint8_t* leftBuf( uint64_t read_pos ) const {
			return buf + ( read_pos - start_pos );
		}
		~CacheEntry() { clear(); }
	};

	fs::path file_name;
	MemoryManager &memory_manager;
	ConsumerEntry * consumer_idx;
	const uint32_t buf_size;
	std::unique_ptr<FileDisk> disk;
	uint64_t write_position = 0;
	uint64_t read_position = 0;
	std::mutex file_sync;
	CacheEntry* cache = nullptr;
	CacheEntry* cache_end = nullptr;

	inline void addToCache( uint8_t* buf, uint32_t buf_size ) {
		auto next = new CacheEntry( buf, buf_size, write_position );
		if( cache == nullptr ) cache_end = cache = next;
		else {
			cache_end->next = next;
			cache_end = next;
		}
	}

	inline void moveNextCache(){
		memory_manager.release( cache->buf_size, this );
		auto next = cache->next;
		delete cache;
		cache = next;
	}

	inline void DiskWrite( uint64_t pos, uint8_t * buf, uint32_t buf_size ){
		if( !disk ) {
			disk.reset( new FileDisk(file_name) );
			file_name.clear();
		}

		disk->Write( pos, buf, buf_size );
	}

	uint32_t Read( uint8_t * buf, const uint32_t &buf_size ){

		uint32_t to_read;

		if( cache != nullptr ){
			assert( read_position < (cache->start_pos + cache->buf_size) );
			if( read_position == cache->start_pos ) {
				// this function cannot be called when buf_size equals current cache buffer!

				if( buf_size < cache->buf_size ){
					memcpy( buf, cache->buf, buf_size );
					read_position += buf_size;
					return buf_size;
				}

				// buf_size > cache[cache_idx].buf_size
				to_read = cache->buf_size;
				memcpy( buf, cache->buf, to_read );
				read_position += to_read;
				moveNextCache();

				return to_read + Read( buf + to_read, buf_size - to_read );
			}

			if( read_position > cache->start_pos ) {
				// read from inside cache buffer

				if( buf_size > cache->leftSize(read_position) ){
					to_read = cache->leftSize(read_position);
					memcpy( buf, cache->leftBuf( read_position), to_read );
					read_position += to_read;
					moveNextCache();
					return to_read + Read( buf + to_read, buf_size - to_read );
				}

				memcpy( buf, cache->leftBuf( read_position ), buf_size );
				read_position += buf_size;
				return buf_size;
			}

			// read_position < cache[cache_idx].start_pos
			to_read = std::min( cache->start_pos - read_position, (uint64_t)buf_size );
		}
		else
			to_read = std::min( write_position - read_position, (uint64_t)buf_size );

		if( !to_read ) return 0; // nothing to read

		disk->Read( read_position, buf, to_read );
		read_position += to_read;

		return to_read + (to_read >= buf_size ? 0 : Read( buf+to_read, buf_size - to_read ) );
	}
};

// creates or cached or simple file stream depends of settings of memory manager
IReadWriteStream * CreateFileStream( const fs::path &fileName, MemoryManager &memory_manager, uint32_t base_buf_size ){
	return memory_manager.CacheEnabled ?
				((IReadWriteStream*)(new CachedFileStream( fileName, memory_manager, base_buf_size )))
			: new FileStream( fileName );
}

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
		, size_to_read( this->compact?(sequence_start_bit>=0?
																		 BUF_SIZE // sequenct compacter read by this value
																	 :(BUF_SIZE/(entry_size-1)*(entry_size-1))) // if just bucket_no removed than read amount should be alligned to entry_size-1
																:(BUF_SIZE/entry_size*entry_size) ) // without compaction read amount should be allinged to entry size
		, buffer_size(
				this->compact ? ( (sequence_start_bit>=0 ? SequenceCompacterReader::getMaxRestoredBufferSize(entry_size-1) : size_to_read )/(entry_size-1)*entry_size + entry_size )
											: size_to_read
				)
	{	}

	uint32_t MaxBufferSize() const { return buffer_size; }
	const std::string getFileName() const {return fileName;}


	void Write( const uint8_t * buf, const uint32_t &buf_size ){
		assert( (buf_size % entry_size_) == 0 );
		if( buf_size == 0 ) return; // nothing to write

		// std::lock_guard<std::mutex> lk( sync_mutex );
		if( !disk_output ){
			bucket_file.reset( CreateFileStream( fileName, memory_manager, size_to_read ) );
			disk_output.reset( new NotFreeingWriteStream( bucket_file.get() ) );
			if( compact ){
				if( sequence_start_bit >= 0 )
					disk_output.reset( compacter = new SequenceCompacterWriter( disk_output.release(), entry_size_-1,
																sequence_start_bit - (sequence_start_bit > begin_bits_ ? 8 : 0 ) ) );

				disk_output.reset( new ByteCutterStream( disk_output.release(),
																	entry_size_, begin_bits_ - log_num_buckets_ ) );
			}

			disk_output.reset( new BufferedWriter( disk_output.release(), entry_size_ ) );
		}

		((BufferedWriter*)(disk_output.get()))->Write( buf, buf_size );
	}

	// The write is NOT thread save and should be synced from outside!!!!
	void Write( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ){
		Write( buf.get(), buf_size );
	}


	// this should be NEVER called when object used in threads!!!
	void FlusToDisk(){
		// std::lock_guard<std::mutex> lk( sync_mutex );
		if( disk_output ) {
			disk_output->Close();
			disk_output.reset();
		}
	}

	// This can be called from thread safly without additional syncs.
	uint32_t Read( std::unique_ptr<uint8_t[]> &buf ){
		{
			if( !bucket_file && !disk_input ) return 0;// nothing has been written
			uint32_t read_size = 0;

			{ // reading file part is locked by mutex
				std::lock_guard<std::mutex> lk( sync_mutex );

				if( disk_output ){ // output not flushed to disk or released
					std::unique_ptr<uint8_t[]> obuf;
					auto buf_size = ((BufferedWriter*)disk_output.get())->ReleaseBuffer( obuf );
					if( buf_size ) {
						memcpy( buf.get(), obuf.get(), buf_size );
						return buf_size;
					}
				}

				if( !disk_input ){
					if( sequence_start_bit >= 0 ){
						// create compacter without input stream to grow custom buffer.
						restorer.reset( new SequenceCompacterReader( nullptr, entry_size_-1,
																		sequence_start_bit - (sequence_start_bit > begin_bits_ ? 8 : 0 ) ) );
					}

					inserter.reset( new ByteInserterStream( nullptr,
																		entry_size_, begin_bits_ - log_num_buckets_,
																		bucket_no_ >> (log_num_buckets_-8) ) );


					if( compacter && disk_output ){ // output not flushed to disk
						// we can get last compacted buffer without writting it to disk
						read_size = compacter->ReleaseBuffer( buf.get() );
						compacter = nullptr;
					}

					if( disk_output )
						disk_output.reset(); // close output stream

					disk_input.reset( (IReadDiskStream*)bucket_file.release() );
				}


				if( read_size == 0 )
					read_size = disk_input->Read( buf, size_to_read );
			}

			// this part is not locked and can be done in many threads in parallel
			if( compact ){
				if( sequence_start_bit >= 0 ){
					std::unique_ptr<uint8_t[]> restore_buffer = std::make_unique<uint8_t[]>( buffer_size );
					read_size = restorer->GrowCustomBuffer( buf.get(), read_size, restore_buffer.get(), buffer_size );
					buf.swap( restore_buffer );
				}

				assert( read_size%(entry_size_-1) == 0 ); // is read size aligned
				assert( (read_size/(entry_size_-1)*entry_size_ ) <= buffer_size ); // is buffer has enough space

				read_size = inserter->GrowBuffer( buf.get(), read_size );
			}

			assert( read_size <= buffer_size );

			return read_size;
		}
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
	std::unique_ptr<IReadWriteStream> bucket_file;
	std::unique_ptr<IWriteDiskStream> disk_output;
	std::unique_ptr<IReadDiskStream> disk_input;
	const uint16_t bucket_no_;
	const uint8_t log_num_buckets_;
	std::unique_ptr<std::thread> disk_io_thread;
	const uint16_t entry_size_;
	const uint16_t begin_bits_;
	const int8_t sequence_start_bit;
	std::mutex sync_mutex;

	const bool compact;
	const uint32_t size_to_read;
	const uint32_t buffer_size;

	SequenceCompacterWriter * compacter = nullptr;
	std::unique_ptr<SequenceCompacterReader> restorer;
	std::unique_ptr<ByteInserterStream> inserter;

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
