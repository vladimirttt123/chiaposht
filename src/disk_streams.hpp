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



struct BucketStream{
	BucketStream( const std::string &fileName, uint16_t bucket_no, uint8_t log_num_buckets, uint16_t entry_size, uint16_t bits_begin )
		: disk( new FileDisk(fileName) )
		, bucket_no_( bucket_no )
		, log_num_buckets_( log_num_buckets )
		, entry_size_(entry_size)
		, bits_begin_(bits_begin)
		, buffer_size( BUF_SIZE/entry_size*entry_size )
		, buffer( new uint8_t[buffer_size] )
	{

	}

	uint64_t ReadPosition() const { return disk_read_position; }
	uint64_t WritePosition() const { return disk_write_position; }
	uint32_t MaxBufferSize() const { return buffer_size; }

	void Write( std::unique_ptr<uint8_t[]> &buf, uint32_t buf_size ){
		assert( (buf_size % entry_size_) == 0 );

		if( buf_size == 0 ) return;
		std::lock_guard<std::mutex> lk( sync_mutex );
		if( disk_io_thread ) disk_io_thread->join();
		buffer.swap( buf );

		assert( Util::ExtractNum64( buffer.get(), bits_begin_-log_num_buckets_, log_num_buckets_ ) == bucket_no_ );

		disk_io_thread.reset( new std::thread( [this, buf_size](){
			disk->Write( disk_write_position, buffer.get(), buf_size );
			disk_write_position += buf_size;
		}));
	}

	void EndToWrite(){
		if( disk_io_thread ) disk_io_thread->join();
		buffer.reset();
		disk->Close();
		disk_io_thread.reset();
	}

	void StartToRead(){
		if( disk_io_thread ) disk_io_thread->join();
		if( disk_read_position < disk_write_position ){
			buffer.reset( new uint8_t[buffer_size] );
			disk_io_thread.reset( new std::thread([this](){ReadBuffer();}) );
		}
	}

	uint32_t Read( std::unique_ptr<uint8_t[]> &buf ){
		assert( disk_read_position <= disk_write_position );
		std::lock_guard<std::mutex> lk( sync_mutex );

		if( disk_read_position >= disk_write_position )
			return 0;

		assert( disk_io_thread );
		disk_io_thread->join();

		uint32_t to_read = std::min( buffer_size, (uint32_t)(disk_write_position - disk_read_position) );
		assert( Util::ExtractNum64( buffer.get(), bits_begin_-log_num_buckets_, log_num_buckets_ ) == bucket_no_ );
		buf.swap( buffer );
		assert( Util::ExtractNum64( buf.get(), bits_begin_-log_num_buckets_, log_num_buckets_ ) == bucket_no_ );

		disk_read_position += to_read;
		if( disk_read_position >= disk_write_position ){
			disk->Remove();
			disk.reset();
			buffer.reset();
			disk_io_thread.reset();
		}
		else
			disk_io_thread.reset( new std::thread([this](){ReadBuffer();}) );

		assert( Util::ExtractNum64( buf.get(), bits_begin_-log_num_buckets_, log_num_buckets_ ) == bucket_no_ );

		return to_read;
	}

	// It can be called if full reading not finished yet to close and remove resources
	void EndToRead(){
		if( disk_io_thread ) disk_io_thread->join();
		disk_read_position = disk_write_position;
		if( disk ) disk->Remove();
		buffer.reset();
		disk_io_thread.reset();
	}

private:
	std::unique_ptr<FileDisk> disk;
	const uint16_t bucket_no_;
	const uint8_t log_num_buckets_;
	std::unique_ptr<std::thread> disk_io_thread;
	uint64_t disk_write_position = 0;
	uint64_t disk_read_position = 0;
	const uint16_t entry_size_;
	const uint16_t bits_begin_;
	std::mutex sync_mutex;

	const uint32_t buffer_size;
	std::unique_ptr<uint8_t[]> buffer;

	void ReadBuffer(){
		uint32_t to_read = std::min( buffer_size, (uint32_t)(disk_write_position - disk_read_position) );
		assert( to_read > 0 );
		disk->Read( disk_read_position, buffer.get(), to_read );
		assert( Util::ExtractNum64( buffer.get(), bits_begin_-log_num_buckets_, log_num_buckets_ ) == bucket_no_ );
	}
};

#endif  // SRC_CPP_DISK_STREAMS_HPP_
