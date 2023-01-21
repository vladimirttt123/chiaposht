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

#ifndef SRC_CPP_SORTING_BUCKET_HPP_
#define SRC_CPP_SORTING_BUCKET_HPP_

#include <mutex>
#include "disk.hpp"
#include "util.hpp"
#include "quicksort.hpp"
#include "exceptions.hpp"



struct SortingBucket{
	SortingBucket( const std::string &fileName, uint32_t entry_size, uint32_t bits_begin, uint8_t bucket_bits_count )
		:disk( new FileDisk(fileName) )
		, bucket_bits_count_(bucket_bits_count)
		, entry_size_(entry_size)
		, bits_begin_(bits_begin)
		, statistics( new uint32_t[1<<bucket_bits_count] )
		, disk_buffer_size( (BUF_SIZE/entry_size)*entry_size )
		, disk_buffer( new uint8_t[disk_buffer_size] )
	{
		// clear statistics
		memset( statistics.get(), 0, sizeof(uint32_t)*(1<<bucket_bits_count) );
	}

	inline uint64_t Size() const { return entry_size_*entries_count; }
	inline uint64_t Count() const { return entries_count; }
	inline uint8_t	SubBucketBits() const { return bucket_bits_count_; }
	inline bool isSorted() const { return !!memory_; }

	uint64_t read_time = 0;
	uint64_t sort_time = 0;

	void AddEntry( const uint8_t *entry, uint32_t statistics_bits ){
		assert(disk); // Check not closed
		assert( statistics_bits < ((uint32_t)1<<bucket_bits_count_) );
		assert( Util::ExtractNum( entry, entry_size_, bits_begin_, bucket_bits_count_ ) == statistics_bits );

		statistics[statistics_bits]++;
		memcpy( disk_buffer.get() + BufferSize(), entry, entry_size_ );
		entries_count++;
		if( BufferSize() == 0 ){
			if( !disk ) throw InvalidStateException("Bucket is already closed");
			disk->Write( Size() - disk_buffer_size, disk_buffer.get(), disk_buffer_size );
			disk->Flush();
		}
	}

	void Flush(){
		if( !disk ) return;
		auto size = BufferSize();
		if( size > 0 ){
			disk->Write( Size() - (uint64_t)size, disk_buffer.get(), size );
			disk->Flush();
		}
	}

	void FreeMemory(){
		memory_.reset();
	}

	/* Like destructor totaly removes the bucket including underlying file */
	void Remove(){
		if( !disk ) return; // already removed
		disk->Remove();
		disk.reset();
		statistics.reset();
		disk_buffer.reset();
	}

	const inline uint8_t* get() const{
		assert(memory_);
		return memory_.get();
	}

	inline void FillBuckets( uint8_t* memory, const uint8_t* disk_buffer, const uint64_t &buf_size, uint64_t * bucket_positions, int64_t direction ){
		for( uint32_t buf_ptr = 0; buf_ptr < buf_size; buf_ptr += entry_size_ ){
			uint64_t bucket_bits = Util::ExtractNum( disk_buffer + buf_ptr, entry_size_, bits_begin_, bucket_bits_count_ );
			memcpy( memory + bucket_positions[bucket_bits], disk_buffer + buf_ptr, entry_size_ );
			bucket_positions[bucket_bits] += direction;
		}
	}

	void SortToMemory( uint32_t num_threads = 2 ){
		if( memory_ ) return; // already sorted;

		memory_.reset( new uint8_t[Size()] );
		uint8_t* memory = memory_.get();

		// TODO: It is possible not flush but use last buffer to split it by buckets.
		Flush();

		uint32_t buckets_count = 1<<bucket_bits_count_;
		auto bucket_positions = std::make_unique<uint64_t[]>( buckets_count );
		bucket_positions[0] = 0;

		// Calculate initial buckets positions.
		for( uint32_t i = 1; i < buckets_count; i++ )
			bucket_positions[i] = bucket_positions[i-1] + (statistics[i-1]*entry_size_);

		// Last position should end at count.
		assert( bucket_positions[buckets_count-1]/entry_size_ + statistics[buckets_count-1] == Count() );

		auto start_time = std::chrono::high_resolution_clock::now();
		// Read from file to buckets
		uint64_t read_pos = 0, read_size = Size();
		if( num_threads <= 1 ){
			while( read_pos < read_size ){
				auto buf_size = std::min( (uint64_t)disk_buffer_size, read_size - read_pos );
				disk->Read( read_pos, disk_buffer.get(), buf_size );
				read_pos += buf_size;

				FillBuckets( memory, disk_buffer.get(), buf_size, bucket_positions.get(), entry_size_ );
			}
			assert( bucket_positions[0] == statistics[0]*entry_size_ ); // check first bucket is full
			assert( bucket_positions[buckets_count-1]/entry_size_ == Count() ); // check last bucket is full
		} else {
			// read in 2 threads one fills from start and second from back.
			auto back_bucket_positions = std::make_unique<uint64_t[]>( buckets_count );
			for( uint32_t i = 0; i < buckets_count - 1; i++ )
				back_bucket_positions[i] = bucket_positions[i+1]-entry_size_;
			back_bucket_positions[buckets_count-1] = back_bucket_positions[buckets_count-2] + statistics[buckets_count-1] * entry_size_;

			std::mutex mut;
			auto thread_func = [&read_pos, &read_size, this, &memory, &mut]( uint8_t* disk_buffer, uint64_t* bucket_positions, const uint64_t direction ){
				while( true ){
					uint64_t buf_size;
					{
						std::lock_guard<std::mutex> lk(mut);
						buf_size = std::min( (uint64_t)disk_buffer_size, read_size - read_pos );
						if( buf_size == 0 ) return;
						disk->Read( read_pos, disk_buffer, buf_size );
						read_pos += buf_size;
					}

					FillBuckets( memory, disk_buffer, buf_size, bucket_positions, direction );
				}
			};

			std::thread forward = std::thread( thread_func, this->disk_buffer.get(), bucket_positions.get(), entry_size_ );
			auto back_buf = std::make_unique<uint8_t[]>(disk_buffer_size);
			std::thread backward = std::thread( thread_func, back_buf.get(), back_bucket_positions.get(), -(int64_t)entry_size_ );

			forward.join();
			backward.join();
#ifndef NDEBUG
			for( uint32_t i = 0; i < buckets_count; i++ )
				assert( (int64_t)bucket_positions[i] == (int64_t)back_bucket_positions[i] + entry_size_ );
#endif
			// Fix bucket_positions for sorting step
			bucket_positions[0] = statistics[0]*entry_size_;
			for( uint32_t i = 1; i < buckets_count; i++ )
				bucket_positions[i] = bucket_positions[i-1] + statistics[i]*entry_size_;
		}

		auto end_time = std::chrono::high_resolution_clock::now();
		read_time = (end_time - start_time)/std::chrono::milliseconds(1);

		if( num_threads <= 1 ){
			// Sort first bucket
			QuickSort::Sort( memory, entry_size_, statistics[0], bits_begin_ + bucket_bits_count_ );
			for( uint32_t i = 1; i < buckets_count; i++ ){
				assert( bucket_positions[i] == (bucket_positions[i-1] + statistics[i]*entry_size_) ); // check any bucket is full
				QuickSort::Sort( memory + bucket_positions[i-1], entry_size_, statistics[i], bits_begin_ + bucket_bits_count_ );
			}
		}else {
			// Sort in threads
			auto threads = std::make_unique<std::thread[]>(num_threads);
			for( uint32_t i = 0; i < num_threads; i++ )
				threads[i] = std::thread( [&num_threads, &buckets_count, this, &bucket_positions, &memory]( uint32_t thread_no ){
						for( uint32_t i = thread_no; i < buckets_count; i += num_threads ){
							assert( i == 0 || bucket_positions[i] == (bucket_positions[i-1] + statistics[i]*entry_size_) ); // check any bucket is full
							QuickSort::Sort( memory + (i == 0 ? 0 : bucket_positions[i-1]), entry_size_, statistics[i], bits_begin_ + bucket_bits_count_ );
						}
					}, i );

			for( uint32_t i = 0; i < num_threads; i++ )
				threads[i].join();
		}

		sort_time = (std::chrono::high_resolution_clock::now() - start_time)/std::chrono::milliseconds(1);
#ifndef NDEBUG
		// Check sort
		for( uint64_t i = 1; i < entries_count; i++ )
			assert( Util::MemCmpBits( memory + (i-1) * entry_size_, memory + i*entry_size_, entry_size_, bits_begin_ ) < 0 );
#endif
	}

private:
	std::unique_ptr<FileDisk> disk;
	const uint8_t bucket_bits_count_;
	const uint32_t entry_size_;
	uint32_t bits_begin_;
	std::unique_ptr<uint32_t[]> statistics;
	const uint32_t disk_buffer_size;
	std::unique_ptr<uint8_t[]> disk_buffer;
	uint64_t entries_count = 0;
	std::unique_ptr<uint8_t[]> memory_;

	inline uint32_t BufferSize() const {
		return Size()%disk_buffer_size;
	}
};

#endif // SRC_CPP_SORTING_BUCKET_HPP_
