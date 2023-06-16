// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//    http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef SRC_CPP_LAST_TABLE_FILE_HPP_
#define SRC_CPP_LAST_TABLE_FILE_HPP_

#include "pos_constants.hpp"
#include "bitfield.hpp"
#include "bitfield_index.hpp"
#include "disk.hpp"
#include "stream_buffer.hpp"
#include <mutex>
#include <atomic>

using namespace std::chrono_literals; // for time units;

struct LastTableWriter {
	LastTableWriter( FileDisk * file, uint8_t k, uint16_t entry_size,
									 bool withCompaction, bool withScan, bool full_scan )
		: k(k), entry_size(entry_size), pos_offset_size( k + kOffsetSize ),
			with_scan(withScan), with_compaction(withCompaction), disk(file)
	{
		if( withScan && full_scan )
			bitmap.reset(  new bitfield( ((uint64_t)2)<<k /* POSSIBLE SIZE PROBLEM */ ) );
	}

	const inline uint64_t getWritePosition() { return write_position; }

	void Write( const uint64_t &write_pos, uint8_t * buf, const uint32_t &buf_size ){
		if( with_scan ) scan( buf, buf_size );

		if( with_compaction ){
			uint32_t to_write = compact(buf, buf_size);

			while( to_write > 0 ){
				if( file_sync.try_lock() ){
					if( write_pos == write_position ){
						disk->Write( disk_position, buf, to_write );
						disk_position += to_write;
						write_position += buf_size;
						to_write = 0; // writted
					}
					file_sync.unlock();
				}
				if( to_write > 0 )
					std::this_thread::sleep_for( 10ns );
			}
		} else {
			// no compaction
			std::lock_guard<std::mutex> lk(file_sync);
			disk->Write( write_pos, buf, buf_size );
		}
	}

	void Close(){

		disk->Close();

		if( with_scan ){
			if( bitmap ){
				auto cres = bitmap->MoveToTable7();
				std::cout << "Compacting bitfield of last table: " << ( cres ? " COMPACTED " : " FAILED!!! " ) << std::endl;
				// if( !cres ) bitmap.resize( table_size );
			}
			else
				bitmap.reset( new bitfield( write_position/entry_size, max_entry_index ) );

			bitmap->FlushToDisk( (disk->GetFileName() + ".bitfield.tmp").c_str() );
		}
	}
private:
	const uint8_t k;
	const uint16_t entry_size;
	const uint8_t  pos_offset_size;
	const bool with_scan, with_compaction;

	uint64_t write_position = 0, disk_position = 0;
	std::mutex file_sync;
	FileDisk * disk;
	std::unique_ptr<bitfield> bitmap;
	std::atomic_uint_fast64_t max_entry_index = 0;

	void scan( uint8_t * buf, const uint32_t &buf_size ){

		for( int64_t buf_ptr = 0; buf_ptr < buf_size; buf_ptr += entry_size ){
			int64_t entry_pos_offset = Util::SliceInt64FromBytes( buf + buf_ptr, k, pos_offset_size);
			uint64_t entry_pos = entry_pos_offset >> kOffsetSize;
			uint64_t entry_offset = entry_pos_offset & ((1U << kOffsetSize) - 1);
			if( bitmap ){
				// TODO lock bitmap instead for non GNU systems.
				bitmap->setTS( entry_pos );
				bitmap->setTS( entry_pos + entry_offset );
			}
			else{
				uint64_t top = entry_pos + entry_offset, cur = max_entry_index.load(std::memory_order_relaxed );
				while( top > cur && !max_entry_index.compare_exchange_weak( cur, top, std::memory_order_relaxed ) )
					cur = max_entry_index.load(std::memory_order_relaxed );
			}
		}
	}

	uint32_t compact( uint8_t * buf, const uint32_t &buf_size ){
		return buf_size;
	}
};


struct LastTableReader : Disk {

	LastTableReader( FileDisk * file, const uint8_t &k, const uint16_t &entry_size,
									 const uint64_t &num_entries,	bool withCompaction )
		: k(k), entry_size(entry_size), f7_shift( 128 - k ), pos_offset_size( k + kOffsetSize )
		, t7_pos_offset_shift( f7_shift - pos_offset_size ), num_entries(num_entries)
		, disk(file), buffer( BUF_SIZE/entry_size*entry_size )
	{

	}

	uint8_t const* Read( uint64_t begin, uint64_t length) override{
		if( !bitfield_ ){ // first read
			bitfield_.reset( new bitfield( num_entries, disk->GetFileName() + ".bitfield.tmp" ) );
			index.reset( new bitfield_index( *bitfield_.get() ) );
			buffer.ensureSize( BUF_SIZE/entry_size*entry_size );
		}

		assert( length == entry_size ); // alway reads by one entry
		assert( begin >= buffer_start_pos );
		while( begin >= buffer_start_pos + buffer.used() ){
			buffer_start_pos += buffer.used();
			buffer.setUsed( std::min( BUF_SIZE/entry_size*entry_size, num_entries*entry_size - buffer_start_pos ) );

			assert( buffer.used() > 0 ); // we should be here when it is something to read;

			disk->Read( buffer_start_pos, buffer.get(), buffer.used() );


			// Rewrite data
			for( uint64_t buf_ptr = 0; buf_ptr < buffer.used(); buf_ptr += entry_size ){

				uint8_t * entry = buffer.get() + buf_ptr;

				// table 7 is special, we never drop anything, so just build
				// next_bitfield
				uint64_t entry_f7 = Util::SliceInt64FromBytes(entry, 0, k);
				uint64_t entry_pos_offset = Util::SliceInt64FromBytes(entry, k, pos_offset_size);

				uint64_t entry_pos = entry_pos_offset >> kOffsetSize;
				uint64_t entry_offset = entry_pos_offset & ((1U << kOffsetSize) - 1);

				// assemble the new entry and write it to the sort manager

				// map the pos and offset to the new, compacted, positions and
				// offsets
				std::tie(entry_pos, entry_offset) = index->lookup(entry_pos, entry_offset);
				entry_pos_offset = (entry_pos << kOffsetSize) | entry_offset;

				// table 7 is already sorted by pos, so we just rewrite the
				// pos and offset in-place
				uint8_t tmp[16];
				uint128_t new_entry = (uint128_t)entry_f7 << f7_shift;
				new_entry |= (uint128_t)entry_pos_offset << t7_pos_offset_shift;
				Util::IntTo16Bytes( tmp, new_entry );
				memcpy( entry, tmp, entry_size );
			}
		}

		return buffer.get() + begin - buffer_start_pos;
	};
	void Write(uint64_t begin, const uint8_t *memcache, uint64_t length) override {
		throw InvalidStateException( "Write impossible to read only disk" );
	};

	void Truncate(uint64_t new_size) override {
		throw InvalidStateException( "Truncate impossible to read only disk" );
	};

	std::string GetFileName() override { return disk->GetFileName(); };
	void FreeMemory() override { buffer.reset(); index.reset(); bitfield_.reset(); };

private:
	const uint8_t k;
	const uint16_t entry_size;
	const uint8_t f7_shift;
	const uint8_t pos_offset_size;
	const uint8_t t7_pos_offset_shift;
	const uint64_t num_entries;

	FileDisk * disk;
	std::unique_ptr<bitfield> bitfield_;
	std::unique_ptr<bitfield_index> index;
	StreamBuffer buffer;
	uint64_t buffer_start_pos = 0;
};

#endif // SRC_CPP_LAST_TABLE_FILE_HPP_
