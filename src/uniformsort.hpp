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

#ifndef SRC_CPP_UNIFORMSORT_HPP_
#define SRC_CPP_UNIFORMSORT_HPP_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "./disk.hpp"
#include "./util.hpp"

namespace UniformSort {

		inline uint64_t const BUF_SIZE = 1024*1024;

    inline static bool IsPositionEmpty(const uint8_t *memory, uint32_t entry_len)
    {

			uint32_t i = 0;
			while( (i + 8) <= entry_len ){
					if( ((uint64_t *)(memory+i))[0] != 0 ) return false;
					i += 8;
			}
    	
//    	if( i < entry_len ) return ((((uint64_t *)memory)[i>>3])>>(64-((entry_len-i)<<3))) == 0;
    	
			while( (i + 4) <= entry_len ){
					if( ((uint32_t *)(memory+i))[0] != 0 ) return false;
					i += 4;
			}
    	
			while( i < entry_len )
				if( memory[i++] != 0 ) return false;
	    
//			for (uint32_t i = 0; i < entry_len; i++)
//					if (memory[i] != 0)  return false;

			return true;
    }

		inline uint32_t CheckSort( uint8_t *memory, uint64_t memory_len, uint32_t entry_len, uint32_t const bits_begin, uint32_t entries_number ){
			uint32_t entries_found = 0;
			uint64_t last_entry_found_at = 0;
			for( uint32_t i = 0; i < memory_len; i+=entry_len ){
				if( !IsPositionEmpty(memory + i, entry_len) ){
					entries_found++;
					if( last_entry_found_at > 0 &&
							Util::MemCmpBits( memory + last_entry_found_at, memory + i, entry_len, bits_begin) > 0)
						return 1;
					last_entry_found_at = i;
				}
			}

			return entries_found == entries_number ? 0 : entries_found + 1;
		}

		inline void SwapEntry( uint8_t *mem, uint8_t *buf, uint8_t *space, uint32_t const & entry_len ){
			memcpy( space, mem, entry_len );
			memcpy( mem, buf, entry_len );
			memcpy( buf, space, entry_len );
		}

		/**
			* Eeach thread works on his own destination memeory ( his own buckets) than we do not need thread synchronization.
			* If some value should be inserted outside the working memory it is not inserted and left as rejected in buffer
			* the main thread after all passes the buffer to insert rejects.
		 */
		void BufferThread ( uint8_t *memory, uint64_t memory_len, const uint8_t * buffer, uint64_t buffer_size, uint32_t entry_len, uint32_t const bits_begin, uint64_t bucket_length, uint64_t thread_mask, uint64_t thread_works_on, std::vector<uint64_t> *rejects )
		{
			auto const swap_space_ = std::make_unique<uint8_t[]>(entry_len);
			auto const swap_space = swap_space_.get();
//			uint32_t written_entries = 0;

			for( uint64_t i = 0; i < buffer_size; i += entry_len ){

				// can't sort if some empty value inside
				assert( !IsPositionEmpty(buffer + i, entry_len) );

				// First unique bits in the entry give the expected position of it in the sorted array.
				// We take 'bucket_length' bits starting with the first unique one.
				uint64_t bucket_bits = Util::ExtractNum( buffer + i, entry_len, bits_begin, bucket_length);

				// Check is this value should be processed by this thread
				if( (bucket_bits&thread_mask) == thread_works_on ){
					// Position to insert current value
					uint64_t pos = bucket_bits*entry_len;

					// Find first empty spot
					while( !IsPositionEmpty(memory + pos, entry_len)
								 && ( ( pos/entry_len )&thread_mask ) == thread_works_on // we still in this thread memory
								 && pos < memory_len)
						pos += entry_len;

					if(( ( pos/entry_len )&thread_mask ) != thread_works_on || pos >= memory_len ){
						// The empty spot is overlimit of current thread memory
						rejects->push_back( i );
//						std::cout << "new reject of "<< thread_works_on << " on " << i << " bits " << bucket_bits << ", pos_bits " << ((pos/entry_len)&thread_mask) << std::endl;
					} else {
						// Push the entry in the first free spot.
						memcpy( memory + pos, buffer + i, entry_len );

						// Fix order of inserted
						for( ; pos > bucket_bits*entry_len && Util::MemCmpBits( memory + pos - entry_len, memory + pos, entry_len, bits_begin) > 0; pos -= entry_len )
							SwapEntry( memory + pos - entry_len, memory + pos, swap_space, entry_len );
					}
				}
			}
//			std::cout << "Thread " << thread_works_on << " writes " << written_entries << std::endl;
		}



    inline void SortToMemory(
        FileDisk &input_disk,
        uint64_t const input_disk_begin,
        uint8_t *const memory,
        uint32_t const entry_len,
        uint64_t const num_entries,
				uint32_t const bits_begin,
				uint32_t const num_threads = 2 )
    {
        uint64_t const memory_len = Util::RoundSize(num_entries) * entry_len;
				auto const swap_space = std::make_unique<uint8_t[]>(entry_len);

				uint64_t buf_size = std::max((uint64_t)1,(uint64_t)num_threads)*(BUF_SIZE / (uint64_t)entry_len)*entry_len;
				BufferedReader buf_reader( &input_disk, input_disk_begin, buf_size , num_entries * entry_len );

				uint64_t bucket_length = 0;
        // The number of buckets needed (the smallest power of 2 greater than 2 * num_entries).
        while ((1ULL << bucket_length) < 2 * num_entries) bucket_length++;
        memset(memory, 0, memory_len);

				uint64_t max_threads = num_threads, thread_mask = 0, thread_mask_shift = 0;
				if( max_threads <= 1 ){
					max_threads = 1;
				}else {
					uint64_t mask_bits = 0;
					while( max_threads > 1 ){
						max_threads >>= 1;
						mask_bits++;
					}
					thread_mask_shift = bucket_length - mask_bits;
					max_threads = (1 << mask_bits);
					thread_mask = ((1<<(mask_bits))-1) << thread_mask_shift;
				}

				std::cout << " threads: " << max_threads << std::flush;

				auto const threads = std::make_unique<std::thread[]>( max_threads );
				auto const thread_rejects = std::make_unique<std::vector<uint64_t>[]>( max_threads );

				auto start_time = std::chrono::high_resolution_clock::now();

				while( ( buf_size = buf_reader.MoveNextBuffer() ) > 0 )
				{
					// Run threads
					for( uint64_t i = 1; i < max_threads; i++ ){
						threads.get()[i] = std::thread( BufferThread,
																						memory, memory_len,
																						buf_reader.GetBuffer(), buf_size,
																						entry_len, bits_begin, bucket_length,
																						thread_mask, i << thread_mask_shift, thread_rejects.get() + i );
					}

					// do work in current thread
					BufferThread( memory, memory_len, buf_reader.GetBuffer(), buf_size,
										entry_len, bits_begin, bucket_length, thread_mask, 0, thread_rejects.get() );

					// Wait for threads
					for( uint64_t i = 1; i < max_threads; i++ )
						threads.get()[i].join();

					// Process thread rejects
					for( uint64_t i = 0; i < max_threads; i++ ){
						for( uint64_t buf_ptr : thread_rejects.get()[i] ){
//							std::cout << "Process Reject of thread " << i << " on " << buf_ptr << std::endl; // debug

							// First unique bits in the entry give the expected position of it in the sorted array.
							// We take 'bucket_length' bits starting with the first unique one.
							uint64_t pos = Util::ExtractNum( buf_reader.GetBuffer() + buf_ptr, entry_len, bits_begin, bucket_length) * entry_len;
							uint64_t empty_pos = pos;

							// As long as position is occupied by a previous entry...
							while( !IsPositionEmpty( memory + empty_pos, entry_len) && empty_pos < memory_len)
								empty_pos += entry_len;

							if( empty_pos >= memory_len )
								throw "Cannot sort by uinform sort: may be the input is not uniform distributed";

							// Move entries if they are bigger than current
							while( empty_pos > pos
										 && Util::MemCmpBits( memory + empty_pos - entry_len,  buf_reader.GetBuffer() + buf_ptr, entry_len, bits_begin) > 0 ){
								memcpy( memory + empty_pos, memory + empty_pos - entry_len, entry_len );
								empty_pos -= entry_len;
							}

							// Push the entry in the free spot.
							memcpy(memory + empty_pos,  buf_reader.GetBuffer() + buf_ptr, entry_len);
						}

						// Now rejects can be cleared.
						thread_rejects.get()[i].clear();
					}

//					auto schk = CheckSort( memory, memory_len, entry_len, bits_begin, (buf_reader.GetBufferStartPosition()+buf_reader.BufferSize() - input_disk_begin)/entry_len );
//					switch( schk ){
//						case 0: break;
//						case 1: std::cout << " !sort incorrect order! " << std::endl; break;
//						default: std::cout << " !sort incorrect entries count! expected: "
//															 << (buf_reader.GetBufferStartPosition()+buf_reader.BufferSize())/entry_len
//															 << ", got: " << (schk-1) << std::endl; break;
//					}
					assert( CheckSort( memory, memory_len, entry_len, bits_begin, (buf_reader.GetBufferStartPosition()+buf_reader.BufferSize() - input_disk_begin)/entry_len ) == 0 );
				}

				auto end_time = std::chrono::high_resolution_clock::now();
				std::cout << ", read time: " << (end_time - start_time)/std::chrono::milliseconds(1)/1000.0 << "s" << std::flush;

				//std::cout << " entr_len " << entry_len <<  ", entries " << num_entries << ", swaps " << swaps << ", ";
				assert( CheckSort(memory, memory_len, entry_len, bits_begin, num_entries ) == 0 );

        uint64_t entries_written = 0;
        // Search the memory buffer for occupied entries.
				for (uint64_t pos = 0; entries_written < num_entries && pos < memory_len; pos += entry_len) {
            if (!IsPositionEmpty(memory + pos, entry_len)) {
                // We've found an entry.
								// write the stored entry itself.
								memcpy( memory + entries_written * entry_len, memory + pos, entry_len );
                entries_written++;
            }
        }

				assert(entries_written == num_entries);
    }


}

#endif  // SRC_CPP_UNIFORMSORT_HPP_
