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

#ifndef SRC_CPP_BUCKETSORT_HPP_
#define SRC_CPP_BUCKETSORT_HPP_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <mutex>

#include "./disk.hpp"
#include "./util.hpp"
#include "./quicksort.hpp"

namespace BucketSort {

		struct BucketInfo{
			uint64_t start_ptr;
			uint64_t end_ptr;
			uint64_t current_ptr;

			bool inline HasSpace() { return current_ptr < end_ptr; }
		};

		inline bool FreeBucektBack( uint8_t * const memory, BucketInfo * const buckets, uint32_t const entry_len, int64_t const min_bucket, int64_t const from_bucket ){
			int64_t empty_bucket = from_bucket - 1;
			while( empty_bucket >= min_bucket && !buckets[empty_bucket].HasSpace() )
				empty_bucket--;

			if( empty_bucket < min_bucket ) return false;

			for( empty_bucket++; empty_bucket <= from_bucket; empty_bucket++ ){
				buckets[empty_bucket].start_ptr = (buckets[empty_bucket-1].end_ptr -= entry_len);
				buckets[empty_bucket].current_ptr -= entry_len; // free one space at current
				if(  buckets[empty_bucket].current_ptr !=  buckets[empty_bucket].start_ptr )
					memcpy( memory + buckets[empty_bucket].start_ptr, memory + buckets[empty_bucket].current_ptr, entry_len ); // move last of bucket to start
			}

			return true;
		}


		inline bool FreeBucektForward( uint8_t * const memory, BucketInfo * const buckets, uint32_t const entry_len, int64_t const max_bucket, int64_t const from_bucket ){
			int64_t empty_bucket = from_bucket + 1;
			while( empty_bucket < max_bucket && !buckets[empty_bucket].HasSpace() )
				empty_bucket++;

			if( empty_bucket >= max_bucket ) return false;

			for( ; empty_bucket > from_bucket; empty_bucket-- ){
				if(  buckets[empty_bucket].current_ptr !=  buckets[empty_bucket].start_ptr )
					memcpy( memory + buckets[empty_bucket].current_ptr, memory + buckets[empty_bucket].start_ptr, entry_len ); // move last of bucket to start
				buckets[empty_bucket-1].end_ptr = (buckets[empty_bucket].start_ptr += entry_len);
				buckets[empty_bucket].current_ptr += entry_len;
			}

			return true;
		}

		/**
			* Eeach thread works on his own destination memeory ( his own buckets) than we do not need thread synchronization.
			* If some value should be inserted outside the working memory it is not inserted and left as rejected in buffer
			* the main thread after all passes the buffer to insert rejects.
		 */
		void BufferThread ( uint8_t *memory, uint64_t memory_len,
												const uint8_t *buffer, uint64_t const  buffer_size, uint32_t const entry_len,
												uint32_t const bits_begin, uint64_t const bucket_length,
												BucketInfo *buckets, uint64_t const start_bucket, uint64_t const end_bucket,
												std::vector<uint64_t> *rejects )
		{
			for( uint64_t i = 0; i < buffer_size; i += entry_len ){

				// First unique bits in the entry give the expected position of it in the sorted array.
				// We take 'bucket_length' bits starting with the first unique one.
				uint64_t bucket_bits = Util::ExtractNum( buffer + i, entry_len, bits_begin, bucket_length);

				// Check is this value should be processed by this thread
				if( bucket_bits >= start_bucket && bucket_bits < end_bucket ){
					auto cur_bucket = buckets + bucket_bits;
					// If Bucket has place insert into it or space can be freed in current thread memory
					if( cur_bucket->HasSpace()
							|| FreeBucektBack( memory, buckets, entry_len, start_bucket, bucket_bits )
							|| FreeBucektForward( memory, buckets, entry_len, end_bucket, bucket_bits ) ){
						assert( cur_bucket->HasSpace() );
						memcpy( memory + cur_bucket->current_ptr, buffer + i, entry_len );
						cur_bucket->current_ptr += entry_len;
					}else {
						//We have a problem - current bucket does not have place
						rejects->push_back( i );
					}
				}
			}
		}

		void SortBucketsThreadMut( uint8_t * const memory, uint32_t const entry_len, uint32_t const bits_begin,
												BucketInfo * const all_buckets, uint64_t const num_buckets, std::mutex * sort_mutex, uint64_t *cur_bucket ){

			while( *cur_bucket < num_buckets ){
				BucketInfo *bucket;
				{
					std::lock_guard<std::mutex> lk( *sort_mutex );
					if( *cur_bucket >= num_buckets ) return; // all done
					bucket = all_buckets + *cur_bucket;
					++ *cur_bucket;
				}

				// Check bucket has content
				if( bucket->start_ptr < bucket->current_ptr )
					QuickSort::Sort( memory + bucket->start_ptr, entry_len, (bucket->current_ptr - bucket->start_ptr)/entry_len, bits_begin );
				// Mark bucket sorted
				bucket->end_ptr = bucket->start_ptr;
			}
		}

		inline uint32_t CheckSort( uint8_t *memory, uint32_t entry_len, uint32_t const bits_begin, uint32_t entries_number ){
			for( uint32_t i = entry_len; i < entries_number*entry_len; i+=entry_len ){
					if( Util::MemCmpBits( memory + i - entry_len, memory + i, entry_len, bits_begin) > 0)
						return 1;

					//if( Util::MemCmpBits( memory + i - entry_len, memory + i, entry_len, bits_begin) == 0 ) std::cout << "equals!!" << std::endl;
				}

			return 0;
		}

		inline void SortToMemory(
				FileDisk &input_disk,
				uint64_t const input_disk_begin,
				uint8_t *const memory,
				uint64_t const	memory_size,
				uint32_t const entry_len,
				uint64_t const num_entries,
				uint32_t const bits_begin,
				uint32_t const num_threads = 2 )
		{
			if( num_entries < 4096 ){
//				std::cout << num_entries << " below minumum entries fall for QS ";
				QuickSort::SortToMemory( input_disk, input_disk_begin, memory, entry_len, num_entries, bits_begin, num_threads );
				return;
			}
				// Align memory size to entry length
				uint64_t const memory_len = entry_len * (memory_size/entry_len);
				uint64_t buf_size = std::max( (uint64_t)1, (uint64_t)num_threads ) * std::min( memory_len, (BUF_SIZE / (uint64_t)entry_len)*entry_len );

				BufferedReader buf_reader( &input_disk, input_disk_begin, buf_size , num_entries * entry_len );

				uint64_t bucket_length = 0;
				while ((1ULL << bucket_length) < memory_len/entry_len ) bucket_length++;
				//bucket_length = std::max( (uint64_t)1, std::min( (uint64_t)10, bucket_length - 10 ) ); // at most 1024 buckets
				bucket_length = std::max( (uint64_t)1, bucket_length - 11 ); // around 2048 entries per bucket
				uint64_t num_buckets = 1ULL << bucket_length;

				uint64_t max_threads = std::max( (uint64_t)1, std::min( num_buckets, (uint64_t)num_threads ) );
				uint64_t buckets_per_thread = num_buckets / max_threads;
				auto const threads = std::make_unique<std::thread[]>( max_threads );
				auto const thread_rejects = std::make_unique<std::vector<uint64_t>[]>( max_threads );

				auto const buckets_u = std::make_unique<BucketInfo[]>( num_buckets );
				auto const buckets = buckets_u.get();
				for( uint64_t i = 0, len = entry_len*(memory_len/num_buckets/entry_len); i < num_buckets; i++ )
				{
					buckets[i].current_ptr = buckets[i].start_ptr = i * len;
					buckets[i].end_ptr = std::min( memory_len, (i+1) * len );
				}

				std::cout << " threads: " << max_threads;
				auto start_time = std::chrono::high_resolution_clock::now();

				while( (buf_size = buf_reader.MoveNextBuffer() ) > 0 )
				{
					// Run threads
					uint64_t start_bucket = 0;
					for( uint64_t i = 0; i < max_threads - 1; i++ ){
						threads.get()[i] = std::thread( BufferThread,
																						memory, memory_len,
																						buf_reader.GetBuffer(), buf_size,
																						entry_len, bits_begin, bucket_length,
																						buckets, start_bucket, start_bucket + buckets_per_thread,
																						thread_rejects.get() + i );
						start_bucket += buckets_per_thread;
					}

					// do work in current thread
					BufferThread( memory, memory_len, buf_reader.GetBuffer(), buf_size,
												entry_len, bits_begin, bucket_length,
												buckets, start_bucket, num_buckets,
												thread_rejects.get() + max_threads );

					// Wait for threads
					for( uint64_t i = 0; i < max_threads - 1; i++ )
						threads.get()[i].join();

					// Process thread rejects
					for( uint64_t i = 0; i < max_threads; i++ ){
						for( uint64_t buf_ptr : thread_rejects.get()[i] ){

							uint64_t bucket_bits = Util::ExtractNum( buf_reader.GetBuffer() + buf_ptr, entry_len, bits_begin, bucket_length );

							if( !FreeBucektBack( memory, buckets, entry_len, 0, bucket_bits )
									&& !FreeBucektForward( memory, buckets, entry_len, num_buckets, bucket_bits ) )
								throw "sort cannot handle"; // this is impossible case because it should be 1.5 times memory we need

							memcpy( memory + buckets[bucket_bits].current_ptr, buf_reader.GetBuffer() + buf_ptr, entry_len );
							buckets[bucket_bits].current_ptr += entry_len;
						}

						// Now rejects can be cleared.
						thread_rejects.get()[i].clear();
					}
				}

				auto end_time = std::chrono::high_resolution_clock::now();
				std::cout << ", read time: " << (end_time - start_time)/std::chrono::milliseconds(1)/1000.0 << "s" << std::flush;

				// Second pass: Sort each bucket
				uint64_t cur_bucket = 0;
				std::mutex sort_mutex;
				for( uint64_t i = 0; i < max_threads; i++ ){
					threads.get()[i] = std::thread( SortBucketsThreadMut, memory, entry_len,
																					// Addion of bucket length can leads to less memory to compare and speed up the process
																					bits_begin + bucket_length,
																					buckets, num_buckets, &sort_mutex, &cur_bucket );
				}

				// Collect sorted buckets to the begining.
				for( uint64_t i = 1, cur_ptr = buckets[0].current_ptr; i < num_buckets; i++ ){
					// Wait for sort is finished: sort thread at the end of the sorting sets end equals start
					while( buckets[i].start_ptr != buckets[i].end_ptr )
						std::this_thread::sleep_for(1ms);

					// The bucket is sorted - move it to the correct position
					for( uint64_t ptr = buckets[i].start_ptr; ptr < buckets[i].current_ptr; ptr += entry_len ){
						memcpy( memory + cur_ptr, memory + ptr, entry_len );

						// check sort order -> next entry bigger than previous one
						assert( cur_ptr == 0 || Util::MemCmpBits(memory + cur_ptr - entry_len, memory + cur_ptr, entry_len, bits_begin ) <= 0 );

						// Move next position
						cur_ptr += entry_len;
					}
				}

				// Wait for threads
				for( uint64_t i = 0; i < max_threads; i++ )
					threads.get()[i].join();

				assert( CheckSort( memory, entry_len, bits_begin, num_entries ) == 0 );
		}

}

#endif  // SRC_CPP_BUCKETSORT_HPP_
