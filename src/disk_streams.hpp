// Copyright 2018 Chia Network Inc

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
			buffer.reset( new uint8_t[last_buf_size = buf_size] );
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

	~AsyncCopyStreamWriter(){ if( io_thread ) io_thread->join(); }
private:
	std::unique_ptr<IWriteDiskStream> disk;
	std::unique_ptr<uint8_t[]> buffer;
	std::unique_ptr<std::thread> io_thread;
	uint32_t last_buf_size = 0;
};
// ====================================================
struct AsyncStreamReader : public IReadDiskStream {
	AsyncStreamReader( IReadDiskStream * disk, uint32_t buf_size )
		: max_buffer_size(buf_size), disk(disk)
	{
	}

	const uint32_t max_buffer_size;

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
		: disk( disk )
		, entry_size(entry_size), begin_bits_( begin_bits&7 ), begin_bytes( begin_bits>>3 )
		, mask( 0xff00 >> begin_bits_ )
	{
	}

	void Write( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override {
		if( last_buffer_size != buf_size )
			buffer.reset( new uint8_t[(last_buffer_size=buf_size)/entry_size*(entry_size+1)] );
		disk->Write( buffer, CompactBuffer( buf.get(), buf_size ) );
	}

	void Close() override {
		if(disk){ disk->Close(); disk.reset(); }
		buffer.reset();
	}

	~SequenceCompacterWriter(){if(disk) disk->Close();}
private:
	std::unique_ptr<IWriteDiskStream> disk;
	uint64_t last_write_value = 0;
	uint32_t last_buffer_size = 0;
	std::unique_ptr<uint8_t[]> buffer;
	const uint16_t entry_size;
	const uint8_t begin_bits_;
	const uint8_t begin_bytes;
	const uint8_t mask;

	uint32_t CompactBuffer( uint8_t *buf, uint32_t buf_size ){
		assert( (buf_size%entry_size) == 0 );

		uint32_t buf_idx = 0;
		auto cur_buf = buffer.get();
#define setNext( v ) cur_buf[buf_idx++] = v

		for( uint32_t i = 0; i < buf_size; ){
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
				((uint32_t*)(cur_buf + buf_idx))[0] = bswap_32( diff );
				buf_idx += 4;
			}	else {
				// save 255 & full value
				setNext( 255 );
				((uint32_t*)(cur_buf + buf_idx))[0] = bswap_32( next_val );
				buf_idx += 4;
			}

			// Process cuted entry
			for( uint32_t j = 0; j < begin_bytes; j++ )
				setNext( buf[i++] );
			setNext( (buf[i]&mask) | (buf[i+4]&~mask) );
			i += 5;
			for( uint32_t j = begin_bytes + 5; j < entry_size; j++ )
				setNext( buf[i++] );

			// set result
			last_write_value = next_val;
		}

		return buf_idx;
	}
};

struct SequenceCompacterReader : public IReadDiskStream{
	SequenceCompacterReader( IReadDiskStream *disk, const uint16_t &entry_size, uint8_t begin_bits )
		: disk( disk )
		, entry_size(entry_size), begin_bits_( begin_bits&7 ), begin_bytes( begin_bits>>3 )
		, mask( 0xff00 >> begin_bits_ )
	{
	}

	uint32_t Read( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override{
		return GrowBuffer( buf.get(), buf_size );
	}

	inline bool atEnd() const override { return buf_fillness <= processed_of_buffer && disk->atEnd(); }
	void Close() override { if(disk){ disk->Close(); disk.reset(); }buffer.reset(); }

	~SequenceCompacterReader(){if(disk) disk->Close();}
private:
	std::unique_ptr<IReadDiskStream> disk;
	uint64_t last_value = 0;
	std::unique_ptr<uint8_t[]> buffer;
	uint32_t processed_of_buffer = 0;
	uint32_t buf_fillness = 0;
	const uint16_t entry_size;
	const uint8_t begin_bits_;
	const uint8_t begin_bytes;
	const uint8_t mask;

	uint32_t GrowBuffer( uint8_t * buf, const uint32_t &buf_size ){
		assert( buf_size%entry_size  == 0 );

		uint32_t i = 0, value_to_insert;
		while( i < buf_size && !atEnd() ){

			uint8_t compress_code = NextByte();

			if( compress_code < 128 )
				value_to_insert = last_value + compress_code;
			else if( compress_code < 192 ){
				value_to_insert = last_value + 128 + ((((uint32_t)(compress_code - 128))<<8) | NextByte());
			} else if( compress_code < 224 ){
				value_to_insert = last_value + 16512 + ((((uint32_t)(compress_code-192))<<16) | (NextByte()<<8) | NextByte());
			} else if( compress_code < 240 ){
				value_to_insert = last_value + 2113664 + ((((uint32_t)(compress_code-224))<<24) | (NextByte()<<16) | (NextByte()<<8) | NextByte());
			} else {
				assert( compress_code == 255 );
				value_to_insert = (NextByte()<<24) | (NextByte()<<16) | (NextByte()<<8) | NextByte();
			}

			// recreate the entry
			for( uint32_t j = 0; j < begin_bytes; j++ )
				buf[i++] = NextByte();
			auto splitted = NextByte();
			buf[i++] = (splitted&mask)|(value_to_insert>>(24+begin_bits_));
			buf[i++] = value_to_insert >> (16+begin_bits_);
			buf[i++] = value_to_insert >> (8+begin_bits_);
			buf[i++] = value_to_insert >> begin_bits_;
			buf[i++] = (value_to_insert<<(8-begin_bits_)) | (splitted&~mask);
			for( uint32_t j = begin_bytes + 5; j < entry_size; j++ )
				buf[i++] = NextByte();

			last_value = value_to_insert;
		}

		return i;
	}

	inline uint32_t NextByte(){
		if( buf_fillness <= processed_of_buffer ){
			if( !buffer ) buffer.reset( new uint8_t[BUF_SIZE] );
			buf_fillness = disk->Read( buffer, BUF_SIZE );
			processed_of_buffer = 0;
		}

		assert( buf_fillness > 0 );

		return buffer.get()[processed_of_buffer++];
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

	uint32_t CompactBuffer( uint8_t *buf, uint32_t buf_size ){
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

	~ByteInserterStream(){ if(disk) disk->Close();}
private:
	std::unique_ptr<IReadDiskStream> disk;
	const uint16_t entry_size;
	const uint16_t begin_bits;
	const uint8_t bytes_begin;
	const uint8_t removed_byte;
	const uint8_t mask;

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
};


struct BufferedWriter : IWriteDiskStream{
	BufferedWriter( IWriteDiskStream* disk, const uint32_t &entry_size )
		: disk(disk), allocated_size( BUF_SIZE/entry_size*entry_size ), buffer( new uint8_t[allocated_size] )
	{}

	void Write( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ){
		uint32_t processed = 0;
		while( processed < buf_size ){
			uint32_t to_process = std::min( allocated_size - buffer_used, buf_size - processed );
			memcpy( buffer.get() + buffer_used, buf.get() + processed, to_process );
			processed += to_process;
			buffer_used += to_process;
			if( buffer_used == allocated_size ){
				disk->Write( buffer, buffer_used );
				buffer_used = 0;
			}
		}
	};
	void Close(){
		if(buffer_used) disk->Write( buffer, buffer_used );
		buffer_used = 0;
		disk->Close();
		buffer.reset();
		disk.reset();
	};

private:
	std::unique_ptr<IWriteDiskStream> disk;
	const uint32_t allocated_size;
	uint32_t buffer_used = 0;
	std::unique_ptr<uint8_t[]> buffer;
};


struct RamStream : IWriteDiskStream, IReadDiskStream{

	void Write( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override{
		assert( read_position == 0 ); // can't write after read
		for( uint64_t done = 0; done < buf_size; ){
			if( (write_position%BUF_SIZE) == 0 )
				ram.push( new uint8_t[BUF_SIZE] );
			auto to_add = std::min( ram.size()*BUF_SIZE - write_position, buf_size - done );
			memcpy( ram.back(), buf.get(), to_add );
			write_position += to_add;
			done += to_add;
		}
	}

	uint32_t Read( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override{
		for( uint64_t done = 0; done < buf_size; ){
			auto to_add = std::min( std::min( write_position - read_position, BUF_SIZE - (read_position%BUF_SIZE) ), buf_size - done );
			if( to_add == 0 ){
				delete ram.front();
				ram.pop();
				return done;
			}
			memcpy( buf.get(), ram.front(), to_add );
			done += to_add;
			read_position += to_add;
			if( (read_position%BUF_SIZE) == 0 ){
				// Clear ram after read;
				delete ram.front();
				ram.pop();
			}
		}
		return 0;
	};
	bool atEnd() const override { return write_position <= read_position; };
	void Close() override {};

	~RamStream(){
		for( ; !ram.empty(); ram.pop() )
			delete ram.front();
	}
private:
	std::queue<uint8_t*> ram;
	uint64_t write_position = 0;
	uint64_t read_position = 0;
};

// This stream not garantee same read order as write order
struct BucketStream{
	BucketStream( std::string fileName, uint16_t bucket_no, uint8_t log_num_buckets, uint16_t entry_size,
								uint16_t begin_bits, bool compact = true, int8_t sequence_start_bit = -1 )
		: fileName(fileName)
		, bucket_no_( bucket_no )
		, log_num_buckets_( log_num_buckets )
		, entry_size_(entry_size)
		, begin_bits_(begin_bits)
		, sequence_start_bit(sequence_start_bit)
		, buffer_size( BUF_SIZE/entry_size*entry_size )
		, compact( compact&(log_num_buckets_>7) )
	{}

	uint32_t MaxBufferSize() const { return buffer_size; }

	void Write( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ){
		assert( (buf_size % entry_size_) == 0 );
		if( buf_size == 0 ) return;

		std::lock_guard<std::mutex> lk( sync_mutex );
		if( !disk_output ){
			bucket_file.reset( new FileDisk( fileName ) );
			disk_output.reset( new WriteFileStream( bucket_file.get() ) );
			if( compact ){
				if( sequence_start_bit >= 0 )
					disk_output.reset( new SequenceCompacterWriter( disk_output.release(), entry_size_-1,
																sequence_start_bit - (sequence_start_bit > begin_bits_ ? 8 : 0 ) ) );

				disk_output.reset( new ByteCutterStream( disk_output.release(),
																	entry_size_, begin_bits_ - log_num_buckets_ ) );
			}
			disk_output.reset( new AsyncStreamWriter( disk_output.release() ) );
		}

		disk_output->Write( buf, buf_size );
	}


	void EndToWrite(){
		std::lock_guard<std::mutex> lk( sync_mutex );
		if( disk_output ) {
			disk_output->Close();
			disk_output.reset();
		}
	}

	uint32_t Read( std::unique_ptr<uint8_t[]> &buf ){
		{
			std::lock_guard<std::mutex> lk( sync_mutex );
			if( !bucket_file ) return 0;// nothing has been written
			if (!disk_input ){
				disk_input.reset( new ReadFileStream( bucket_file.get(), bucket_file->GetWriteMax() ) );
				if( compact ){
					if( sequence_start_bit >= 0 )
						disk_input.reset( new SequenceCompacterReader( disk_input.release(),
																				entry_size_-1, sequence_start_bit - (sequence_start_bit > begin_bits_ ? 8 : 0 ) ) );

					disk_input.reset( new ByteInserterStream( disk_input.release(),
																		entry_size_, begin_bits_ - log_num_buckets_,
																		bucket_no_ >> (log_num_buckets_-8) ) );
				}
				disk_input.reset( new AsyncStreamReader( disk_input.release(), buffer_size ) );
			}

			return disk_input->Read( buf, buffer_size );
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
			bucket_file->Remove();
			bucket_file.reset();
		}
	}

private:
	const std::string fileName;
	std::unique_ptr<FileDisk> bucket_file;
	std::unique_ptr<IWriteDiskStream> disk_output;
	std::unique_ptr<IReadDiskStream> disk_input;
	const uint16_t bucket_no_;
	const uint8_t log_num_buckets_;
	std::unique_ptr<std::thread> disk_io_thread;
	const uint16_t entry_size_;
	const uint16_t begin_bits_;
	const int8_t sequence_start_bit;
	std::mutex sync_mutex;

	const uint32_t buffer_size;
	const bool compact;

};

struct LastTableRewrited : IReadDiskStream {
	LastTableRewrited( IReadDiskStream *disk, const uint8_t k,const uint16_t entry_size,
											const uint64_t num_entries, const fs::path &filename )
		: k(k), f7_shift( 128 - k ), pos_offset_size( k + kOffsetSize )
		, t7_pos_offset_shift( f7_shift - pos_offset_size ), entry_size(entry_size)
		, num_entries(num_entries), filename_(filename), disk(disk)
	{
	}


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

	if( withCompaction && k >= 20 )
		res = new SequenceCompacterReader( res, entry_size, k );

	if( max_buffer_size > 0 )
		res = new AsyncStreamReader( res, max_buffer_size );

	return res;
}

IReadDiskStream * CreateLastTableReader( FileDisk * file, uint8_t k, uint16_t entry_size,
																				 uint64_t num_entries, const fs::path &bitfield_filename,
																				 bool withCompaction, uint32_t max_buffer_size = 0 ){
	IReadDiskStream *res = new ReadFileStream( file, file->GetWriteMax() );

	if( withCompaction && k >= 20 )
		res = new SequenceCompacterReader( res, entry_size, k );

	res = new LastTableRewrited( res, k, entry_size, num_entries, bitfield_filename );

	if( max_buffer_size > 0 )
		res = new AsyncStreamReader( res, max_buffer_size );

	return res;
}


IWriteDiskStream * CreateLastTableWriter( FileDisk * file, uint8_t k, uint16_t entry_size,
																					bool withCompaction, uint32_t max_buffer_size = 0 ){
	IWriteDiskStream * res = new WriteFileStream( file );
	if( withCompaction && k >= 20 )
		res = new SequenceCompacterWriter( res, entry_size, k );

	res = max_buffer_size > 0 ?
				new AsyncStreamWriter( res, max_buffer_size )
			: (IWriteDiskStream*)new AsyncCopyStreamWriter( res );

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
