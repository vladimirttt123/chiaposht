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
	SequenceCompacterStream( IDiskStream &base_stream ) : disk( base_stream ) {
	}

	void Write( uint8_t *buf, const uint32_t &buf_size ){
		disk.Write( buf, buf_size );
	}
	uint32_t Read( uint8_t *buf, const uint32_t &buf_size ){
		return disk.Read( buf, buf_size );
	}

	uint64_t LeftToRead() const { return disk.LeftToRead(); }
	void Close() { disk.Close(); }
	void Remove() { disk.Remove(); }
private:
	IDiskStream &disk;
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
		if( this->compact )
			disk.reset( new ByteRemoverStream( disk.release(),
																				 entry_size, begin_bits_ - log_num_buckets,
																				 bucket_no_ >> (log_num_buckets_-8) ) );
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



// Supposed entry at least 5 bytes long that would be replace by result bytes that with return how many of them
// 5 bytes could be replaced up to 6 bytes
uint8_t SequencedEntryOptimize( const uint8_t * entry, uint16_t start_bit, uint64_t &last_value, uint8_t *result ){
	assert( start_bit < 8 );

	const uint8_t mask = 0xff00>>start_bit;
	uint8_t bytes_to_save = 1;
	uint64_t next_val = Util::ExtractNum64( entry, start_bit, 32 );

	int64_t diff = next_val-last_value;
	if( diff >= 0 && diff <= 252 ){
		// to save one byte of diff
		result[0] = (entry[0]&mask) | (diff>>start_bit);
		if( mask )
			result[1] = (diff<<(8-start_bit)) | (entry[5]&~mask);
	} else if( diff > 252 && diff <= 507 ){
		// to save 2 bytes 253 & diff - 252;
		diff -= 252;
		result[0] = (entry[0]&mask) | (((uint8_t)253)>>start_bit);
		result[1] = (253<<(8-start_bit)) | (((uint8_t)diff)>>start_bit);
		if(mask)
			result[2] = (diff<<(8-start_bit)) | (entry[5]&~mask);
		bytes_to_save = 2;
	} else if( diff > 507 && diff <= 16132 ){
		// to save 3 bytes 254 & diff - 252;
		diff -= 252;
		result[0] = (entry[0]&mask) | (((uint8_t)254)>>start_bit);
		result[1] = (254<<(8-start_bit)) | (diff>>(8+start_bit));
		result[2] = (diff>>start_bit);
		if( mask )
			result[3] = (diff<<(8-start_bit)) | (entry[5]&~mask);

		bytes_to_save = 3;
	} else {
		// to save 5 bytes 255 & next_val
		result[0] = entry[0]|~mask;
		result[1] = entry[0]|mask;
		result[2] = entry[1];
		result[3] = entry[2];
		result[4] = entry[3];
		if( mask ) result[5] = entry[4];
		bytes_to_save = 5;
	}

	// set result
	last_value = next_val;

	return bytes_to_save;
}

uint64_t SequencedEntryReconstruct( const uint8_t * entry, const uint16_t &start_bit, const uint64_t &last_value, uint8_t *result ){
	assert( start_bit < 8 );

	const uint8_t mask = 0xff00>>start_bit;
	auto set_val = [&entry, &mask, &result, &start_bit]( uint64_t val, uint8_t end ){
		result[0] = (entry[0]&mask) | ((val>>(24+start_bit))&~mask);

		return val;
	};
	uint64_t compress_code = Util::ExtractNum64( entry, start_bit, 40 );

	switch( compress_code >> 32 ){
		case 255: // reconstruct full value
			return set_val( compress_code&0xffffffff, entry[5] );

		case 254: // reconstruct from 2 bytes value;
			return set_val( ((compress_code>>16)&0xffff)+252, entry[3] );

		case 253: // reconstruct from 1 bytes value;
			return set_val( ((compress_code>>24)&0xff)+252, entry[2] );

		default: // full value provided
			return set_val( compress_code, entry[1] );
	}
}

#endif  // SRC_CPP_DISK_STREAMS_HPP_
