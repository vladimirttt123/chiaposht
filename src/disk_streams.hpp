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
#include "disk.hpp"
#include "util.hpp"


// Interfase for streams
struct IDiskStream{

	virtual void Write( uint8_t *buf, const uint32_t &buf_size ) = 0;
	virtual uint32_t Read( uint8_t *buf, const uint32_t &buf_size ) = 0;
	virtual uint64_t LeftToRead() const = 0;

	virtual void Close() = 0;
	virtual void Remove() = 0;
	virtual ~IDiskStream() = default;
};


// NOT threadsafe, NOT asynchrounous, NOT buffered
struct SimpleDiskStream : IDiskStream{

	SimpleDiskStream( const std::string &fileName ) : disk( new FileDisk(fileName) )	{	}

	void Write( uint8_t *buf, const uint32_t &buf_size ){
		disk->Write( disk_write_position, buf, buf_size );
		disk_write_position += buf_size;
	}

	uint32_t Read( uint8_t *buf, const uint32_t &buf_size ){
		uint32_t to_read = std::min( (uint64_t)buf_size, disk_write_position - disk_read_position );
		if( to_read )
			disk->Read( disk_read_position, buf, to_read );
		disk_read_position += to_read;
		return to_read;
	}

	uint64_t LeftToRead() const { return disk_write_position - disk_read_position; }

	void Close() { disk->Close(); }
	void Remove() { disk->Remove(); disk.reset(); }
	~SimpleDiskStream(){ if( disk ) disk->Close(); }

private:
	std::unique_ptr<FileDisk> disk;
	uint64_t disk_write_position = 0;
	uint64_t disk_read_position = 0;
};


struct SequenceCompacterStream : public IDiskStream{
	SequenceCompacterStream( IDiskStream *base_stream, const uint32_t &max_buffer_size,
													 const uint16_t &entry_size, uint8_t begin_bits )
		: disk( base_stream )
		, buf_size( max_buffer_size/entry_size*(entry_size+1) )
		, buffer( new uint8_t[buf_size] )
		, entry_size(entry_size), begin_bits_( begin_bits&7 ), begin_bytes( begin_bits>>3 )
		, mask( 0xff00 >> begin_bits_ )
	{
	}

	void Write( uint8_t *buf, const uint32_t &buf_size ){
		num_entries_written += buf_size/entry_size;
		disk->Write( buffer.get(), CompactBuffer( buf, buf_size ) );
	}

	uint32_t Read( uint8_t *buf, const uint32_t &buf_size ){
		return GrowBuffer( buf, buf_size );
	}

	uint64_t LeftToRead() const { return (num_entries_written - num_entries_read)*entry_size; }
	void Close() { disk->Close(); }
	void Remove() { disk->Remove(); disk.reset(); }

private:
	std::unique_ptr<IDiskStream> disk;
	uint64_t last_value = 0;
	uint64_t last_write_value = 0;
	const uint32_t buf_size;
	std::unique_ptr<uint8_t[]> buffer;
	uint32_t processed_of_buffer = 0;
	uint32_t buf_fillness = 0;
	const uint16_t entry_size;
	const uint8_t begin_bits_;
	const uint8_t begin_bytes;
	uint64_t num_entries_written = 0;
	uint64_t num_entries_read = 0;
	const uint8_t mask;

	uint32_t CompactBuffer( uint8_t *buf, uint32_t buf_size ){
		uint32_t buf_idx = 0;
		auto setNext = [&buf_idx,this]( uint8_t v ){ buffer.get()[buf_idx++] = v;};

		for( uint32_t i = 0; i < buf_size; ){
			uint64_t next_val = Util::ExtractNum64( buf + i + begin_bytes, begin_bits_, 32 );
			uint64_t diff = next_val - last_write_value;
			if( diff <= 252 ){
				// to save one byte of diff
				setNext( diff );
			} else if( diff <= 507 ){
				// to save 2 bytes 253 & diff - 252;
				setNext( 253 );
				setNext( diff - 252 );
			} else if( diff <= 16132 ){
				// to save 3 bytes 254 & diff - 252;
				setNext( 254 );
				setNext( (diff -= 252)>>8 );
				setNext( diff );
			} else {
				// save 255 & full value
				setNext( 255 );
				((uint32_t*)(buffer.get() + buf_idx))[0] = bswap_32( next_val );
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

	uint32_t GrowBuffer( uint8_t * buf, uint32_t buf_size ){
		assert( buf_size%entry_size  == 0 );

		uint32_t i = 0, value_to_insert;
		while( i < buf_size && num_entries_read < num_entries_written ){
			num_entries_read++;

			uint8_t compress_code = NextByte();

			switch( compress_code ){
				case 255: // reconstruct full value
					value_to_insert = (NextByte()<<24) | (NextByte()<<16) | (NextByte()<<8) | NextByte();
				break;

				case 254: // reconstruct from 2 bytes value;
					value_to_insert = last_value + (NextByte()<<8) + NextByte() + 252;
				break;

				case 253: // reconstruct from 1 bytes value;
					value_to_insert = last_value + NextByte() + 252;
				break;

				default: // one byte value
					value_to_insert = last_value + compress_code;
				break;
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
			buf_fillness = disk->Read( buffer.get(), buf_size );
			processed_of_buffer = 0;
		}

		return buffer[processed_of_buffer++];
	}
};


struct ByteRemoverStream : public IDiskStream{
	ByteRemoverStream( IDiskStream *base_stream, uint16_t entry_size, uint16_t begin_bits, uint8_t removed_byte )
		: disk( base_stream ), entry_size(entry_size), begin_bits(begin_bits), bytes_begin(begin_bits>>3)
		, removed_byte( removed_byte ), mask( 0xff00 >> (begin_bits&7) )
	{	}

	void Write( uint8_t *buf, const uint32_t &buf_size ){
		disk->Write( buf, CompactBuffer(buf, buf_size) );
	}

	uint32_t Read( uint8_t *buf, const uint32_t &buf_size ){
		return GrowBuffer( buf, disk->Read( buf, buf_size/entry_size*(entry_size-1) ));
	}

	uint64_t LeftToRead() const { return disk->LeftToRead(); }

	void Close() { disk->Close(); }

	void Remove() { disk->Remove(); }

private:
	std::unique_ptr<IDiskStream> disk;
	const uint16_t entry_size;
	const uint16_t begin_bits;
	const uint8_t bytes_begin;
	const uint8_t removed_byte;
	const uint8_t mask;

	uint32_t CompactBuffer( uint8_t *buf, uint32_t buf_size ){
		assert( (buf_size%entry_size) == 0 );

		// extract current bucket
		if( mask ){
			for( uint32_t i = bytes_begin + 1; i < buf_size; i += entry_size ){
				assert( Util::ExtractNum64( buf + i - 1, begin_bits&7, 8 ) == removed_byte );
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

	uint32_t GrowBuffer( uint8_t * buf, uint32_t buf_size ){
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

// This stream not garantee same read order as write order
struct BucketStream{
	BucketStream( const std::string &fileName, uint16_t bucket_no, uint8_t log_num_buckets, uint16_t entry_size, uint16_t begin_bits, bool compact = true, int8_t sequence_start_bit = -1 )
		: disk( new SimpleDiskStream(fileName) )
		, bucket_no_( bucket_no )
		, log_num_buckets_( log_num_buckets )
		, entry_size_(entry_size)
		, begin_bits_(begin_bits)
		, buffer_size( BUF_SIZE/entry_size*entry_size )
		, buffer( new uint8_t[buffer_size] )
		, compact( compact&(log_num_buckets_>7) )
	{
		if( this->compact ){
			if( sequence_start_bit >= begin_bits )
				disk.reset( new SequenceCompacterStream( disk.release(), buffer_size,
																								 entry_size-1, sequence_start_bit - 8 ) );

			disk.reset( new ByteRemoverStream( disk.release(),
																				 entry_size, begin_bits_ - log_num_buckets,
																				 bucket_no_ >> (log_num_buckets_-8) ) );
		}
	}

	uint32_t MaxBufferSize() const { return buffer_size; }

	void Write( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ){
		assert( (buf_size % entry_size_) == 0 );

		if( buf_size == 0 ) return;
		std::lock_guard<std::mutex> lk( sync_mutex );
		if( disk_io_thread ) disk_io_thread->join();
		buffer.swap( buf );

		assert( compact || Util::ExtractNum64( buffer.get(), begin_bits_-log_num_buckets_, log_num_buckets_ ) == bucket_no_ );

		disk_io_thread.reset( new std::thread( [this](uint32_t buf_size){
			disk->Write( buffer.get(), buf_size );
		}, buf_size) );
	}

	void EndToWrite(){
		if( disk_io_thread ) disk_io_thread->join();
		buffer.reset();
		disk->Close();
		disk_io_thread.reset();
	}

	void StartToRead(){
		if( disk_io_thread ) disk_io_thread->join();
		disk_io_thread.reset();
		if( disk->LeftToRead() > 0 ){
			buffer.reset( new uint8_t[buffer_size] );
			ReadNext();
		}
	}

	uint32_t Read( std::unique_ptr<uint8_t[]> &buf ){
		{
			std::lock_guard<std::mutex> lk( sync_mutex );

			if( !disk_io_thread ) return 0;
			disk_io_thread->join();

			assert( compact || Util::ExtractNum64( buffer.get(), begin_bits_-log_num_buckets_, log_num_buckets_ ) == bucket_no_ );
			buf.swap( buffer );
			auto read = last_read_size;
			assert( compact || Util::ExtractNum64( buf.get(), begin_bits_-log_num_buckets_, log_num_buckets_ ) == bucket_no_ );

			if( disk->LeftToRead() > 0 )
				ReadNext();
			else {
				disk->Remove();
				disk.reset();
				buffer.reset();
				disk_io_thread.reset();
			}

			assert( Util::ExtractNum64( buf.get(), begin_bits_-log_num_buckets_, log_num_buckets_ ) == bucket_no_ );
			return read;
		}
	}

	// It can be called if full reading not finished yet to close and remove resources
	void EndToRead(){
		std::lock_guard<std::mutex> lk( sync_mutex );
		if( disk_io_thread ) disk_io_thread->join();
		if( disk ) disk->Remove();
		disk.reset();
		buffer.reset();
		disk_io_thread.reset();
	}

private:
	std::unique_ptr<IDiskStream> disk;
	const uint16_t bucket_no_;
	const uint8_t log_num_buckets_;
	std::unique_ptr<std::thread> disk_io_thread;
	const uint16_t entry_size_;
	const uint16_t begin_bits_;
	std::mutex sync_mutex;

	const uint32_t buffer_size;
	std::unique_ptr<uint8_t[]> buffer;
	uint32_t last_read_size = 0;
	const bool compact;

	inline void ReadNext(){
		disk_io_thread.reset( new std::thread([this](){
			last_read_size = disk->Read( buffer.get(), buffer_size );}) );
	}
};


#endif  // SRC_CPP_DISK_STREAMS_HPP_
