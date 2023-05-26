// Copyright 2020 Chia Network Inc

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//    http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include <algorithm>
#include "bitfield.hpp"

struct bitfield_index
{
    // Cache the number of set bits every kIndexBucket bits.
		// For a bitfield of size 2^32, this means a 2*2^22 + 4*2^16 + 8*2^0 ~ 8.5 MiB index
		// For a bitfield of size 2^38, this means a 2*2^28 + 4*2^22 + 8*2^6 ~ 528 MiB index
		static inline const uint64_t kIndexBucketBits = 10; // i.e. 1024
		static inline const uint64_t kIndexBucket = 1UL<<kIndexBucketBits;

		bitfield_index( uint64_t entries_count ) { setSize( entries_count ); }
		bitfield_index(bitfield &b) { reinit( &b ); }

		void reinit( bitfield *b ){
			bitfield_ = b;
			setSize( b->size() );

			for( uint64_t idx = 0, counter = 0, last_hi = 0, last_mid = 0;
					 idx < uint64_t(bitfield_->size()); idx += kIndexBucket ) {
				assert( (idx >> kIndexBucketBits ) < allocated_size );

				if( (idx&0xffffffff) == 0 )
					last_hi = index_hi[idx>>32] = counter;


				if( (idx&0xffff) == 0 ){
					index_mid[idx>>16] = counter - last_hi;
					last_mid = counter;
				}

				index_lo[idx>>kIndexBucketBits] = counter - last_mid;
				int64_t const left = std::min( uint64_t(bitfield_->size()) - idx, kIndexBucket);
				counter += bitfield_->count(idx, idx + left);
			}
		}

    std::pair<uint64_t, uint64_t> lookup(uint64_t pos, uint64_t offset) const
    {
			uint64_t const bucket_lo = pos >> kIndexBucketBits;
			uint64_t const bucket_mid = pos >> 16;
			uint64_t const bucket_hi = pos >> 32;

			assert(bucket_lo < size);
			assert(pos < uint64_t(bitfield_->size()));
			assert(pos + offset < uint64_t(bitfield_->size()));
			assert(bitfield_->get(pos) && bitfield_->get(pos + offset));

			uint64_t const base = index_hi[bucket_hi] + index_mid[bucket_mid] + index_lo[bucket_lo];

			//assert( base == (uint64_t)bitfield_->count(0, bucket_lo<<kIndexBucketBits ) );

			int64_t const aligned_pos = pos & ~uint64_t(63);

			uint64_t const aligned_pos_count = bitfield_->count(bucket_lo << kIndexBucketBits, aligned_pos);
			uint64_t const offset_count = aligned_pos_count + bitfield_->count(aligned_pos, pos + offset);
			uint64_t const pos_count = aligned_pos_count + bitfield_->count(aligned_pos, pos);

			assert(offset_count >= pos_count);

			return { base + pos_count, offset_count - pos_count };
    }
private:
		static inline uint64_t align_size( uint64_t size ) { return (size-1+(1UL<<kIndexBucketBits))>>kIndexBucketBits; }

		void setSize( uint64_t new_size ){
			size = align_size( new_size );
			if( size > allocated_size ){
				index_lo.reset( new uint16_t[allocated_size = size] );
				index_mid.reset( new uint32_t[1 + (new_size>>16)] );
				index_hi.reset( new uint64_t[1 + (new_size>>32)] );
			}
		}

		bitfield * bitfield_ = nullptr;
		uint64_t allocated_size = 0, size = 0;
		std::unique_ptr<uint64_t[]> index_hi;
		std::unique_ptr<uint32_t[]> index_mid;
		std::unique_ptr<uint16_t[]> index_lo;
};

