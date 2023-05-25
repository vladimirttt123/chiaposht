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
    // For a bitfield of size 2^32, this means a 32 MiB index
		static inline const uint64_t kIndexBucket = 1024;
		static inline uint64_t align_size( uint64_t size ) { return (size-1+kIndexBucket)/kIndexBucket; }

		bitfield_index( uint64_t entries_count ) : bitfield_(nullptr)
			, allocated_size( align_size( entries_count ) ), size(0)
			, index_( new uint64_t[allocated_size] )
		{

		}

		bitfield_index(bitfield &b) : bitfield_(&b)
			, allocated_size( align_size( bitfield_->size() ) ), size(allocated_size)
			, index_( new uint64_t[allocated_size] )
    {
			reinit( &b );
		}

		void reinit( bitfield *b ){
			bitfield_ = b;
			size = align_size( bitfield_->size() );
			if( size > allocated_size )
				index_.reset( new uint64_t[allocated_size = size] );

			for( uint64_t idx = 0, counter = 0; idx < uint64_t(bitfield_->size()); idx += kIndexBucket) {
					assert( idx/kIndexBucket < allocated_size );

					index_[idx/kIndexBucket] = counter;
					int64_t const left = std::min( uint64_t(bitfield_->size()) - idx, kIndexBucket);
					counter += bitfield_->count(idx, idx + left);
			}
		}

    std::pair<uint64_t, uint64_t> lookup(uint64_t pos, uint64_t offset) const
    {
        uint64_t const bucket = pos / kIndexBucket;

				assert(bucket < size);
				assert(pos < uint64_t(bitfield_->size()));
				assert(pos + offset < uint64_t(bitfield_->size()));
				assert(bitfield_->get(pos) && bitfield_->get(pos + offset));

        uint64_t const base = index_[bucket];

        int64_t const aligned_pos = pos & ~uint64_t(63);

				uint64_t const aligned_pos_count = bitfield_->count(bucket * kIndexBucket, aligned_pos);
				uint64_t const offset_count = aligned_pos_count + bitfield_->count(aligned_pos, pos + offset);
				uint64_t const pos_count = aligned_pos_count + bitfield_->count(aligned_pos, pos);

        assert(offset_count >= pos_count);

        return { base + pos_count, offset_count - pos_count };
    }
private:
		bitfield * bitfield_;
		uint64_t allocated_size, size;
		std::unique_ptr<uint64_t[]> index_;
};

