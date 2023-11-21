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

#ifdef __GNUC__
#include <atomic>
#endif

struct bitfield
{
		struct ThreadWriter{
			static const uint32_t BUF_PER_MASK = 1024;

			inline void set(int64_t const &bit){
				uint64_t idx = bit >> 6;
				uint32_t i = idx%syncs_num;
				uint64_t ptr = i*BUF_PER_MASK + counts[i];
				idxs[ptr] = idx;
				masks[ptr] = uint64_t(1) << (bit & 63);
				counts[i]++;
				if( counts[i] >= BUF_PER_MASK ){
					process( i );
					counts[i] = 0;
				}
			}


			ThreadWriter( bitfield & src )
				: src(src), syncs_num(src.syncs_num)
				, counts( new uint16_t[syncs_num] )
				, idxs( new uint64_t[syncs_num*BUF_PER_MASK] )
				, masks( new uint64_t[syncs_num*BUF_PER_MASK] )
			{
				memset( counts.get(), 0, 2*syncs_num );
			}

			~ThreadWriter(){
				// reset rests
				for( uint32_t i = 0; i < syncs_num; i++ )
					if( counts[i] > 0 )
						process( i );
			}
		private:
			bitfield & src;
			const uint16_t syncs_num;

			std::unique_ptr<uint16_t[]> counts;
			std::unique_ptr<uint64_t[]> idxs;
			std::unique_ptr<uint64_t[]> masks;


			void inline process( uint32_t idx ){
				std::lock_guard<std::mutex> lk(src.thread_syncs[idx]);
				for( uint64_t i = 0, ptr = idx*BUF_PER_MASK; i < counts[idx]; i++, ptr++ ){
					src.buffer_.get()[idxs[ptr]] |= masks[ptr];
				}
			}
		};



		explicit bitfield( int64_t size )
//				: buffer_(new uint64_t[(size + 63) / 64])
				: size_((size + 63) / 64), allocated_size(size_)
				, buffer_( Util::allocate<uint64_t>( size_ ) )
		{
			clear();
    }

		explicit bitfield( int64_t size, int64_t table_7_max_entry )
				: size_((size + 63) / 64), buffer_( NULL, [](uint64_t*d){delete[]d;} )
				, table_7_max_entry( table_7_max_entry )
		{}

		// Restore from file
		bitfield( int64_t size, const fs::path &filename, bool with_file_remove = true )
				: size_((size + 63) / 64), buffer_( NULL, [](uint64_t*d){delete[]d;} )
		{
			file_.reset( new FileDisk( filename, false ) );
			file_->Read( 0, (uint8_t*)(&table_7_max_entry), 8 );
			if( table_7_max_entry < 0 ){
				b_file_.reset( new BufferedDisk( file_.get(), size + 8 ) );
			}
			else {
				if( with_file_remove )
					file_->Remove();
				file_.reset();
			}
		}

		// copy partially from other
		bitfield( const bitfield &other, int64_t start_bit, int64_t size )
				: size_((size + 63) / 64), buffer_(NULL,[](uint64_t*d){delete[]d;})
		{
			assert((start_bit % 64) == 0);
			assert( size >= 0 );
			assert( (start_bit + size) <= other.size() );

			if( other.table_7_max_entry >= 0 )
				table_7_max_entry = other.table_7_max_entry - start_bit;
			else{
				buffer_.reset(new uint64_t[ allocated_size = size_ ]);
				if( other.file_ == nullptr )
					memcpy( (uint8_t*)buffer_.get(), ((uint8_t*)other.buffer_.get()) + (start_bit>>3), size_ << 3 );
				else
					other.file_->Read( 8 + (start_bit>>3), (uint8_t*)buffer_.get(), size_<<3 );
			}
		}


		//---------------------------------------------------
		void PrepareToThreads( uint32_t num_threads ){
			assert( !thread_syncs );// assume single call

			syncs_num = 1 << (uint16_t)(1 + std::log2( num_threads ));

			thread_syncs.reset( new std::mutex[syncs_num] );
		}
		// -----------------------------------

		inline void set( int64_t const bit )
    {
//				if( b_file_ != nullptr || table_7_max_entry >= 0 )
//					throw InvalidStateException( "Cannot set in RO bitfield" );
				assert( !b_file_ ); // not flushed to disk
				assert( table_7_max_entry < 0 ); // not table_7 bitfield
				assert(bit / 64 < size_);
				buffer_.get()[bit / 64] |= uint64_t(1) << (bit & 63);
    }

		inline void setTS( int64_t const bit )
		{
//				if( b_file_ != nullptr || table_7_max_entry >= 0 )
//					throw InvalidStateException( "Cannot set in RO bitfield" );
				assert( !b_file_ ); // not flushed to disk
				assert( table_7_max_entry < 0 ); // not table_7 bitfield
				assert(bit / 64 < size_);
#ifdef __GNUC__
				__atomic_fetch_or( buffer_.get() + (bit / 64), uint64_t(1) << (bit & 63), __ATOMIC_RELAXED );
#else // __GNUC__
				std::lock_guard<std::mutex> lk(single_set_sync);
				buffer_.get()[bit / 64] |= uint64_t(1) << (bit & 63);
#endif // __GNUC__
		}


		inline bool get(int64_t const bit) const
    {
			if( table_7_max_entry >= 0 )
				return bit <= table_7_max_entry;

			const auto pos = bit >> 6;
			assert( pos < size_);
			return ( (b_file_ == nullptr ? buffer_.get()[pos]:((uint64_t*)b_file_->Read( 8 + (pos<<3), 8))[0])
							 & (uint64_t(1) << (bit % 64))) != 0;
    }


		inline void clear()
    {
			if( b_file_ != nullptr || table_7_max_entry >= 0 )
				throw InvalidStateException( "Cannot clear RO bitfield" );
			std::memset(buffer_.get(), 0, size_ * 8);
		}

		inline void reset( int64_t new_size ){
			b_file_.reset(); // clear file
			table_7_max_entry = -1; // clear table 7 bitfield
			size_ = (new_size + 63) / 64; // set new size
			if( size_ > allocated_size ) // reallocate mem if need
				//buffer_.reset( new uint64_t[ allocated_size = size_ ] );
				Util::allocate<uint64_t>( allocated_size = size_).swap( buffer_ );

			clear();
		}

		// size is the max number of elements could be stored
		inline int64_t size() const { return size_ * 64; }
		inline uint64_t memSize() const { return allocated_size*8; }
		static inline uint64_t memSize( uint64_t bits ) { return ((bits + 63) / 64)*8; }

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
			FreeMemory();
			return true;
		}

		void RemoveFile(){
			if( b_file_ ){
				b_file_->FreeMemory();
				b_file_.reset();
			}
			if( file_ ){
				file_->Remove();
				file_.reset();
			}
		}

		int64_t FreeMemory()
    {
			int64_t res = allocated_size;
			buffer_.reset();
			size_ = allocated_size = 0;
			return res * 8;
    }

		// return amount of freed ram
		int64_t FlushToDisk( const fs::path &filename, bitfield *transfer_ram_to = nullptr ){
			if( file_ ) return 0;
			auto const length = memSize();
			file_.reset( new FileDisk( filename ) );
			file_->Write( 0, (uint8_t*)(&table_7_max_entry), 8 );
			if( table_7_max_entry < 0 )
				file_->Write( 8, (uint8_t*)buffer_.get(), length );
			file_->Flush();
			b_file_.reset( new BufferedDisk( file_.get(), length + 8 ) );

			int64_t res = 0;
			if( transfer_ram_to == nullptr || transfer_ram_to->allocated_size > 0 ){
				buffer_.reset(); // free memory
				res = allocated_size * 8;
			}
			else {
				buffer_.swap( transfer_ram_to->buffer_ );
				transfer_ram_to->allocated_size = allocated_size;
			}
			allocated_size = 0;

			return res;
		}

		inline bool is_readonly() const { return file_ != nullptr; }
		inline bool is_table_7() const { return table_7_max_entry >= 0; }
private:
		int64_t size_, allocated_size = 0;
		std::unique_ptr<uint64_t,void(*)(uint64_t*)> buffer_;
    // number of 64-bit words

		std::unique_ptr<FileDisk> file_;
		std::unique_ptr<BufferedDisk> b_file_;

		int64_t table_7_max_entry = -1;

		std::mutex sync_mutex;

		uint16_t syncs_num;
		std::unique_ptr<std::mutex[]> thread_syncs; // used for thread writer
#ifndef __GNUC__
		std::mutex single_set_sync;
#endif
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


