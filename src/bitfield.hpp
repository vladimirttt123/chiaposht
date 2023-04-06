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
		explicit bitfield( int64_t size )
        : buffer_(new uint64_t[(size + 63) / 64])
				, size_((size + 63) / 64)
    {
        clear();
    }

		explicit bitfield( int64_t size, int64_t table_7_max_entry )
				: size_((size + 63) / 64), table_7_max_entry( table_7_max_entry )
		{}

		// Restore from file
		bitfield( int64_t size, const fs::path &filename, bool with_file_remove = true )
			: size_((size + 63) / 64)
		{
			FileDisk file = FileDisk( filename, false );
			file.Read( 0, (uint8_t*)(&table_7_max_entry), 8 );
			if( table_7_max_entry < 0 ){
				buffer_.reset( new uint64_t[(size + 63) / 64] );
				file.Read( 8, (uint8_t*)buffer_.get(), size_*8 );
			}
			if( with_file_remove ) file.Remove();
		}

		// copy partially from other
		bitfield( const bitfield &other, int64_t start_bit, int64_t size )
				: size_((size + 63) / 64)
		{
			assert((start_bit % 64) == 0);
			assert( size >= 0 );
			assert( (start_bit + size) <= other.size() );

			if( other.table_7_max_entry >= 0 )
				table_7_max_entry = other.table_7_max_entry - start_bit;
			else{
				buffer_.reset(new uint64_t[(size + 63) / 64]);
				if( other.file_ == nullptr )
					memcpy( (uint8_t*)buffer_.get(), ((uint8_t*)other.buffer_.get()) + (start_bit>>3), size_ << 3 );
				else
					other.file_->Read( 8 + (start_bit>>3), (uint8_t*)buffer_.get(), size_<<3 );
			}
		}



		inline void set(int64_t const bit)
    {
				if( b_file_ != nullptr || table_7_max_entry >= 0 )
					throw InvalidStateException( "Cannot set in RO bitfield" );
				assert(bit / 64 < size_);
				buffer_[bit / 64] |= uint64_t(1) << (bit & 63);
    }

		inline void set( const uint64_t *bits, const uint32_t &count ){
			if( b_file_ != nullptr || table_7_max_entry >= 0 )
				throw InvalidStateException( "Cannot set in RO bitfield" );
			for( uint32_t i = 0; i < count; i++ ){
				assert( bits[i] / 64 < (uint64_t)size_);
				buffer_[bits[i] / 64] |= uint64_t(1) << (bits[i] & 63);
			}

		}

		inline bool get(int64_t const bit) const
    {
			if( table_7_max_entry >= 0 )
				return bit <= table_7_max_entry;

			const auto pos = bit >> 6;
			assert( pos < size_);
			return ( (b_file_ == nullptr ? buffer_[pos]:((uint64_t*)b_file_->Read( 8 + (pos<<3), 8))[0])
							 & (uint64_t(1) << (bit % 64))) != 0;
    }


		inline void clear()
    {
			if( b_file_ != nullptr || table_7_max_entry >= 0 )
				throw InvalidStateException( "Cannot clear RO bitfield" );
			std::memset(buffer_.get(), 0, size_ * 8);
		}

		// size is the max number of elements could be stored
		inline int64_t size() const { return size_ * 64; }
		inline uint64_t memSize() const { return ( file_ == nullptr && table_7_max_entry < 0 )? ((uint64_t)size_ << 3) : 0; }
		static inline uint64_t memSize( uint64_t bits ) { return bits >> 3; }

		inline int64_t count(int64_t const start_bit, int64_t const end_bit) const
    {
        assert((start_bit % 64) == 0);
        assert(start_bit <= end_bit);

				if( table_7_max_entry >= 0 ){
					auto to = end_bit > (table_7_max_entry +1) ? (table_7_max_entry+1) : end_bit;
					return to >= start_bit ? (to-start_bit) : 0;
				}

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
						ret += Util::PopCount( ((uint64_t*)b_file_->Read( 8 +(start<<3), 8))[0] );
					}
					if (tail > 0) {
							uint64_t const mask = (uint64_t(1) << tail) - 1;
							ret += Util::PopCount( ((uint64_t*)b_file_->Read( 8 + ((end_bit>>6)<<3), 8))[0] & mask);
					}
				}
        return ret;
    }

		// Table 7 bitmap is very specific it is full up to some entry and empty after it
		// than such bitmap can be described by one number only.
		bool MoveToTable7(){
			if( table_7_max_entry >= 0 ) return true;
			if( is_readonly() ) return false;

			int64_t i = size_ - 1;
			while( i >= 0 && buffer_.get()[i] == 0)
				i--;
			uint64_t idx = i * 64;
			i--;
			while( i >= 0 && buffer_.get()[i] == 0xffffffffffffffff )
				i--;
			if( i >= 0 ) return false;

			bool isIn = true;
			int64_t max_is = 0;
			for( uint64_t j = 0; j < 65; j++ ){
				if( isIn != get(idx + j) ){
					if( isIn ){
						max_is = idx + j - 1;
						isIn = false;
					}
					else {
						return false;
					}
				}
			}

			table_7_max_entry = max_is;
			buffer_.reset();
			return true;
		}

		void FreeMemory( bool with_file_remove = true )
    {
				buffer_.reset();
				if( b_file_ != nullptr ){
					b_file_->FreeMemory();
					b_file_.reset();
					if( with_file_remove )
						file_->Remove();
					file_.reset();
				}
				size_ = 0;
    }

		void FlushToDisk( const fs::path &filename ){
			if( file_ ) return;
			auto const length = memSize();
			file_.reset( new FileDisk( filename ) );
			file_->Write( 0, (uint8_t*)(&table_7_max_entry), 8 );
			if( table_7_max_entry < 0 )
				file_->Write( 8, (uint8_t*)buffer_.get(), length );
			file_->Flush();
			b_file_.reset( new BufferedDisk( file_.get(), length ) );
			// free memory
			buffer_.reset();
		}

		inline bool is_readonly() const { return file_ != nullptr; }
		inline bool is_table_7() const { return table_7_max_entry >= 0; }
private:
		std::unique_ptr<uint64_t[]> buffer_;
    // number of 64-bit words
    int64_t size_;

		std::unique_ptr<FileDisk> file_;
		std::unique_ptr<BufferedDisk> b_file_;

		int64_t table_7_max_entry = -1;
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

	// Evaluates correct in limits
	int64_t count(int64_t const start_bit, int64_t const end_bit) const 	{
		const bitfield * bf = src_.is_readonly()?reader.get():&src_;
		int64_t start64 = (start+start_bit)&(~(uint64_t)63);
		return bf->count(start64, end_bit + start) - bf->count( start64, start+start_bit );
	}

	inline bool get( int64_t const & bit ) const {
		return src_.is_readonly()?reader->get(bit+start):src_.get( bit + start );
	}

	private:
		const bitfield& src_;
		std::unique_ptr<bitfield> reader;
		uint64_t start;
};


