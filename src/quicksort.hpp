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

#ifndef SRC_CPP_QUICKSORT_HPP_
#define SRC_CPP_QUICKSORT_HPP_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "./disk.hpp"
#include "util.hpp"

namespace QuickSort {

		inline static void SortInner( uint8_t *memory, uint64_t memory_len, uint32_t L, uint32_t bits_begin, uint64_t begin, uint64_t end, uint8_t *pivot_space, uint32_t num_threads = 0 );
		inline static void SortInnerUniform( uint8_t *memory, uint64_t memory_len, uint32_t L, uint32_t bits_begin, uint64_t begin, uint64_t end, uint8_t *pivot_space, uint32_t num_threads = 0 );

		inline static void SortPartsThreaded( uint8_t *memory, uint64_t memory_len, uint32_t entry_len, uint32_t bits_begin,
																					uint64_t begin, uint64_t middle, uint64_t middle_end, uint64_t end,
																					uint32_t add_num_threads, uint8_t *pivot_space, bool is_uniform )
		{
			std::thread *tLeft = NULL;
			uint8_t *inner_pivot_space;
			uint32_t th1 = 0, th2 = 0;

			if( add_num_threads > 1 ){
				// Set more threads for bigger part
				th1 = (uint32_t)std::round( (add_num_threads-1)*((middle - begin)/(double)(middle - begin + end - middle_end) ) );
				th2 = add_num_threads - 1 - th1;
			}

			if( begin < middle ){
				if( add_num_threads > 0 ){
					inner_pivot_space = new uint8_t[entry_len];
					tLeft = new std::thread( is_uniform ? SortInnerUniform : SortInner, memory, memory_len, entry_len, bits_begin, begin, middle, inner_pivot_space, th1 );
					add_num_threads--;
				}else
					if( is_uniform )
						SortInnerUniform( memory, memory_len, entry_len, bits_begin, begin, middle, pivot_space, 0 );
					else
						SortInner( memory, memory_len, entry_len, bits_begin, begin, middle, pivot_space, 0 );
			}

			if( middle_end < end ){
				if( is_uniform )
					SortInnerUniform( memory, memory_len, entry_len, bits_begin, middle_end, end, pivot_space, th2 );
				else
					SortInner( memory, memory_len, entry_len, bits_begin, middle_end, end, pivot_space, th2 );
			}

			if( tLeft != NULL ) {
				tLeft->join();
				delete tLeft;
				delete [] inner_pivot_space;
			}
		}

		inline static void SortInner(
        uint8_t *memory,
        uint64_t memory_len,
        uint32_t L,
        uint32_t bits_begin,
        uint64_t begin,
        uint64_t end,
				uint8_t *pivot_space,
				uint32_t num_threads )
    {
        if (end - begin <= 5) {
            for (uint64_t i = begin + 1; i < end; i++) {
                uint64_t j = i;
                memcpy(pivot_space, memory + i * L, L);
                while (j > begin &&
                       Util::MemCmpBits(memory + (j - 1) * L, pivot_space, L, bits_begin) > 0) {
                    memcpy(memory + j * L, memory + (j - 1) * L, L);
                    j--;
                }
                memcpy(memory + j * L, pivot_space, L);
            }
            return;
        }

        uint64_t lo = begin;
        uint64_t hi = end - 1;

        memcpy(pivot_space, memory + (hi * L), L);
        bool left_side = true;

        while (lo < hi) {
            if (left_side) {
                if (Util::MemCmpBits(memory + lo * L, pivot_space, L, bits_begin) < 0) {
                    ++lo;
                } else {
                    memcpy(memory + hi * L, memory + lo * L, L);
                    --hi;
                    left_side = false;
                }
            } else {
                if (Util::MemCmpBits(memory + hi * L, pivot_space, L, bits_begin) > 0) {
                    --hi;
                } else {
                    memcpy(memory + lo * L, memory + hi * L, L);
                    ++lo;
                    left_side = true;
                }
            }
        }
        memcpy(memory + lo * L, pivot_space, L);

				SortPartsThreaded( memory, memory_len, L, bits_begin, begin, lo, lo + 1, end, num_threads, pivot_space, false );
    }


		inline static void SortInnerUniform(
				uint8_t *memory,
				uint64_t memory_len,
				uint32_t L,
				uint32_t bits_begin,
				uint64_t begin,
				uint64_t end,
				uint8_t *pivot_space,
				uint32_t num_threads )
		{
				if( end - begin <= 128 ){
					SortInner( memory, memory_len, L, bits_begin, begin, end, pivot_space, num_threads );
					return;
				}

				uint64_t lo = begin;
				uint64_t hi = end - 1;
				uint64_t compare_byte = bits_begin >> 3;
				// mask is first bit of comparable part.
				// For uniform input it is supposed half of input has 0 and half of inputs has 1
				uint8_t mask = 1 << (7-(bits_begin&7));

				// find start
				while( hi > lo && (memory[ hi * L + compare_byte ] & mask) != 0 )
					hi--;

				if( hi <= lo ){
					SortInner( memory, memory_len, L, bits_begin, begin, end, pivot_space, num_threads );
					return;
				}

				// store record aside
				memcpy(pivot_space, memory + (hi * L), L);
				bool left_side = true;

				while (lo < hi) {
						if (left_side) {
								if( (memory[ lo * L + compare_byte ] & mask) == 0 ) {
										++lo;
								} else {
										memcpy(memory + hi * L, memory + lo * L, L);
										--hi;
										left_side = false;
								}
						} else {
								if( (memory[ hi * L + compare_byte ] & mask) != 0 ) {
										--hi;
								} else {
										memcpy(memory + lo * L, memory + hi * L, L);
										++lo;
										left_side = true;
								}
						}
				}
				memcpy(memory + lo * L, pivot_space, L);
				lo++; // <== last entry was from lo
				assert( (memory[ lo * L + compare_byte ] & mask) != 0 && (memory[ lo * L - L + compare_byte ] & mask) == 0);

				auto ratio = std::abs( (lo - begin)/(double)(end - begin) - 0.5 );
				bool continue_uniform = ratio < 0.1;
				SortPartsThreaded( memory, memory_len, L, bits_begin + 1, begin, lo, lo, end, num_threads, pivot_space, continue_uniform );
		}



    inline void Sort(
        uint8_t *const memory,
        uint32_t const entry_len,
        uint64_t const num_entries,
        uint32_t const bits_begin)
    {
        uint64_t const memory_len = (uint64_t)entry_len * num_entries;
        auto const pivot_space = std::make_unique<uint8_t[]>(entry_len);
				SortInner(memory, memory_len, entry_len, bits_begin, 0, num_entries, pivot_space.get(), 0);
    }



		inline void SortToMemoryUniform(
				FileDisk &input_disk,
				uint64_t const input_disk_begin,
				uint8_t *const memory,
				uint32_t const entry_len,
				uint64_t const num_entries,
				uint32_t const bits_begin,
				uint32_t num_threads = 4
				)
		{
				uint64_t const memory_len = num_entries * entry_len;
				uint64_t buf_size = std::max((uint64_t)1,(uint64_t)num_threads)*(BUF_SIZE / (uint64_t)entry_len)*entry_len;

				BufferedReader buf_reader( &input_disk, input_disk_begin, buf_size, memory_len );

				uint64_t free_space_start = 0, free_space_end = entry_len*(num_entries-1);
				uint64_t compare_byte = bits_begin >> 3;
				// mask is first bit of comparable part.
				// For uniform input it is supposed half of input has 0 and half of inputs has 1
				uint8_t mask = 1 << (7-(bits_begin&7));

				std::cout << " threads: " << num_threads;
				auto start_time = std::chrono::high_resolution_clock::now();

				// Rading the file and do initial pivoting.
				while( (buf_size = buf_reader.MoveNextBuffer()) > 0 ){
					for( uint64_t buf_ptr = 0; buf_ptr < buf_size; buf_ptr += entry_len ){
						if( (buf_reader.GetBuffer()[buf_ptr + compare_byte]&mask) != 0 ) {
							memcpy( memory + free_space_end, buf_reader.GetBuffer() + buf_ptr, entry_len );
							free_space_end -= entry_len;
						} else {
							memcpy( memory + free_space_start, buf_reader.GetBuffer() + buf_ptr, entry_len );
							free_space_start += entry_len;
						}
					}
				}

				auto end_time = std::chrono::high_resolution_clock::now();
				std::cout << ", read time: " << (end_time - start_time)/std::chrono::milliseconds(1)/1000.0 << "s" << std::flush;

				auto const pivot_space = std::make_unique<uint8_t[]>(entry_len);
				// add to middle 1 if last was from lo side
				auto middle = free_space_start/entry_len + ( (memory[free_space_start + compare_byte]&mask) == 0?1:0 );

				assert( (memory[middle*entry_len + compare_byte]&mask) != 0
						&& (memory[middle*entry_len - entry_len + compare_byte]&mask) == 0 );

				// Count current thread
				if( num_threads > 0 ) num_threads--;
				SortPartsThreaded( memory, memory_len, entry_len, bits_begin + 1 /* 1 bit is already sorted */, 0, middle, middle, num_entries, num_threads, pivot_space.get(), true );
		}



		inline void SortToMemory(
				FileDisk &input_disk,
				uint64_t const input_disk_begin,
				uint8_t *const memory,
				uint32_t const entry_len,
				uint64_t const num_entries,
				uint32_t const bits_begin,
				uint32_t num_threads = 4
				)
		{
				uint64_t const memory_len = num_entries * entry_len;
				uint64_t buf_size = std::max((uint64_t)1,(uint64_t)num_threads)*(BUF_SIZE / (uint64_t)entry_len)*entry_len;

				BufferedReader buf_reader( &input_disk, input_disk_begin, buf_size, memory_len );
				uint64_t free_space_start = 0, free_space_end = entry_len*(num_entries-1);

				std::cout << " threads: " << num_threads;
				auto start_time = std::chrono::high_resolution_clock::now();

				while( (buf_size = buf_reader.MoveNextBuffer()) > 0 ){
					for( uint64_t buf_ptr = 0; buf_ptr < buf_size; buf_ptr += entry_len ){
						if( free_space_start != 0 && Util::MemCmpBits( memory, buf_reader.GetBuffer() + buf_ptr, entry_len, bits_begin ) < 0 ) {
							memcpy( memory + free_space_end, buf_reader.GetBuffer() + buf_ptr, entry_len );
							free_space_end -= entry_len;
						} else {
							memcpy( memory + free_space_start, buf_reader.GetBuffer() + buf_ptr, entry_len );
							free_space_start += entry_len;
						}
					}
				}

				auto end_time = std::chrono::high_resolution_clock::now();
				std::cout << ", read time: " << (end_time - start_time)/std::chrono::milliseconds(1)/1000.0 << "s" << std::flush;

				auto const pivot_space = std::make_unique<uint8_t[]>(entry_len);
				auto middle = free_space_start/entry_len;
				if( num_threads > 0 ) num_threads--;
				SortPartsThreaded( memory, memory_len, entry_len, bits_begin, 0, middle, middle, num_entries, num_threads, pivot_space.get(), false );
		}

}


#endif  // SRC_CPP_QUICKSORT_HPP_
