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

#include <memory>
#include "disk.hpp"


struct bitfield
{
		explicit bitfield(int64_t size)
        : buffer_(new uint64_t[(size + 63) / 64])
				, size_((size + 63) / 64)
				, file_( nullptr)
				, b_file_( nullptr )
    {
        clear();
    }

		bitfield( const bitfield &other, int64_t start_bit, int64_t size )
				: buffer_(new uint64_t[(size + 63) / 64])
				, size_((size + 63) / 64)
				, file_( nullptr)
				, b_file_( nullptr )
		{
			assert((start_bit % 64) == 0);
			assert( size >= 0 );
			assert( (start_bit + size) <= other.size() );

			if( other.file_ == nullptr )
				memcpy( (uint8_t*)buffer_.get(), ((uint8_t*)other.buffer_.get()) + (start_bit>>3), size_ << 3 );
			else
				other.file_->Read( start_bit>>3, (uint8_t*)buffer_.get(), size_<<3 );
		}

		inline void set(int64_t const bit)
    {
				if( b_file_ != nullptr )
					throw InvalidStateException( "Cannot set in RO bitfield" );
				assert(bit / 64 < size_);
				buffer_[bit / 64] |= uint64_t(1) << (bit & 63);
    }

		inline void set( const uint64_t *bits, const uint32_t &count ){
			if( b_file_ != nullptr )
				throw InvalidStateException( "Cannot set in RO bitfield" );
			for( uint32_t i = 0; i < count; i++ ){
				assert( bits[i] / 64 < size_);
				buffer_[bits[i] / 64] |= uint64_t(1) << (bits[i] & 63);
			}

		}

		inline bool get(int64_t const bit) const
    {
				const auto pos = bit >> 6;
				assert( pos < size_);
				return ( (b_file_ == nullptr ? buffer_[pos]:((uint64_t*)b_file_->Read(pos<<3, 8))[0])
								 & (uint64_t(1) << (bit % 64))) != 0;
    }


		inline void clear()
    {
			if( b_file_ != nullptr )
				throw InvalidStateException( "Cannot clear RO bitfield" );
			std::memset(buffer_.get(), 0, size_ * 8);
		}

    int64_t size() const { return size_ * 64; }
		uint64_t memSize() const { return file_ == nullptr? ((uint64_t)size_ << 3) : 0; }
		static inline uint64_t memSize( uint64_t bits ) { return bits >> 3; }

    int64_t count(int64_t const start_bit, int64_t const end_bit) const
    {
        assert((start_bit % 64) == 0);
        assert(start_bit <= end_bit);

				int64_t ret = 0;
				int const tail = end_bit % 64;

				if( file_ == nullptr ){
					uint64_t const* start = buffer_.get() + start_bit / 64;
					uint64_t const* end = buffer_.get() + end_bit / 64;
					while (start != end) {
							ret += Util::PopCount(*start);
							++start;
					}
					if (tail > 0) {
							uint64_t const mask = (uint64_t(1) << tail) - 1;
							ret += Util::PopCount(*end & mask);
					}
				}else{
					for( int64_t start = start_bit >> 6, end = end_bit >> 6; start < end; start++ ){
						ret += Util::PopCount( ((uint64_t*)b_file_->Read(start<<3, 8))[0] );
					}
					if (tail > 0) {
							uint64_t const mask = (uint64_t(1) << tail) - 1;
							ret += Util::PopCount( ((uint64_t*)b_file_->Read((end_bit>>6)<<3, 8))[0] & mask);
					}
				}
        return ret;
    }

		void FreeMemory()
    {
				buffer_.reset();
				if( b_file_ != nullptr ){
					b_file_->FreeMemory();
					delete b_file_;
					b_file_ = nullptr;
					file_->Remove();
					delete file_;
					file_ = nullptr;
				}
				size_ = 0;
    }

		void FlushToDisk( const fs::path &filename ){
			auto const length = memSize();
			file_ = new FileDisk( filename );
			file_->Write( 0, (uint8_t*)buffer_.get(), length );
			file_->Flush();
			b_file_ = new BufferedDisk( file_, length );
			// free memory
			buffer_.reset();
		}

		inline bool is_readonly() const { return file_ != nullptr; }
private:
		std::unique_ptr<uint64_t[]> buffer_;
    // number of 64-bit words
    int64_t size_;

		FileDisk * file_;
		BufferedDisk * b_file_;
};


struct bitfieldReader
{
	bitfieldReader( const bitfield &src ) : src_(src)	{	}

	void setLimits( uint64_t start_bit, uint64_t size ){
		if( src_.is_readonly() ){
			this->start = start_bit%64;
			reader.reset( new bitfield( src_, (start_bit>>6) <<6, size + this->start ) );
		}
		else {
			this->start = start_bit;
		}
	}

	inline bool get( int64_t const & bit ) const {
		return src_.is_readonly()?reader->get(bit+start):src_.get( bit + start );
	}

	private:
		const bitfield& src_;
		std::unique_ptr<bitfield> reader;
		uint64_t start;
};


