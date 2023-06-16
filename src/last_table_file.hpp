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

#include "entry_sizes.hpp"
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
		: k(k), entry_size(entry_size), entry_size_bits( k + k + kOffsetSize /* with bitfield only version */)
		,	pos_offset_size( k + kOffsetSize )
		, with_scan(withScan), with_compaction(withCompaction), disk(file)
	{
		if( withScan && full_scan )
			bitmap.reset(  new bitfield( ((uint64_t)2)<<k /* POSSIBLE SIZE PROBLEM */ ) );
	}

	const inline uint64_t getWritePosition() { return write_position; }

	void Write( const uint64_t &write_pos, uint8_t * buf, const uint32_t &buf_size ){
		if( with_scan ) scan( buf, buf_size );

		if( with_compaction ){
			StreamBuffer cbuf;

			compact(buf, buf_size, cbuf);

			assert( cbuf.used() > 0 );

			while( cbuf.used() > 0 ){
				if( file_sync.try_lock() ){
					if( write_pos == write_position ){
						disk->Write( disk_position, cbuf.get(), cbuf.used() );
						disk_position += cbuf.used();
						write_position += buf_size;
						cbuf.setUsed( 0 ); // writted
					}
					file_sync.unlock();
				}
				if( cbuf.used() > 0 )
					std::this_thread::sleep_for( 10ns );
			}
		} else {
			// no compaction
			std::lock_guard<std::mutex> lk(file_sync);
			disk->Write( write_pos, buf, buf_size ); // just write on mentioned position
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
	const uint64_t entry_size_bits;
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

	static void inline addBits( StreamBuffer & cbuf, ParkBits & bits ){
		bits.ToBytes( cbuf.getEnd() );
		cbuf.addUsed( (bits.GetSize() +7)/8 );
	}

	void compact( uint8_t * buf, const uint32_t &buf_size, StreamBuffer &cbuf ){

		if( (entry_size_bits%8) == 0 ){
			// no compaction
			cbuf.ensureSize(buf_size).setUsed( buf_size );
			memcpy( cbuf.get(), buf, buf_size );
		}else{
			assert( k + kOffsetSize < 64 ); // this works up to this limit only since kOffsetSize == 10 than up to k54
			// compact here
			ParkBits parkStart, parkMid, parkEnd;
			const uint32_t num_entries = buf_size / entry_size;

			cbuf.ensureSize( 4 + (num_entries * k + 7)/8 + (num_entries*(k+kOffsetSize)+7)/8 ).setUsed( 4 );
			memcpy( cbuf.get(), &num_entries, 4 );

			for( uint32_t buf_ptr = 0; buf_ptr < buf_size; buf_ptr += entry_size ){
				if( buf_ptr > 0 && ((buf_ptr/entry_size)%8) == 0 ){
					addBits( cbuf, parkStart );
					addBits( cbuf, parkEnd );
					parkStart.Clear();
					parkEnd.Clear();
				}

				uint64_t start = Util::SliceInt64FromBytes( buf + buf_ptr, k );
				uint64_t end = Util::SliceInt64FromBytes( buf + buf_ptr, k, k + kOffsetSize );

				parkStart.AppendValue( start, k );
				parkEnd.AppendValue( end, k + kOffsetSize );
			}

			if( parkStart.GetSize() > 0 ){
				addBits( cbuf, parkStart );
				addBits( cbuf, parkEnd );
			}
		}
	}
};




struct LastTableReader : Disk {

	LastTableReader( FileDisk * file, const uint8_t &k, const uint16_t &entry_size,
									 const uint64_t &num_entries,	bool withCompaction )
		: k(k), entry_size(entry_size), entry_size_bits( k + k + kOffsetSize /* with bitfield only version */)
		, f7_shift( 128 - k ), pos_offset_size( k + kOffsetSize )
		, t7_pos_offset_shift( f7_shift - pos_offset_size ), num_entries(num_entries), with_compaction(withCompaction)
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

			if( with_compaction && (entry_size_bits%8) )
				UncompactNext();
			else{
				buffer.setUsed( std::min( BUF_SIZE/entry_size*entry_size, num_entries*entry_size - buffer_start_pos ) );

				assert( buffer.used() > 0 ); // we should be here when it is something to read;

				disk->Read( buffer_start_pos, buffer.get(), buffer.used() );
			}

			RewriteBuffer();
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
	const uint64_t entry_size_bits;
	const uint8_t f7_shift;
	const uint8_t pos_offset_size;
	const uint8_t t7_pos_offset_shift;
	const uint64_t num_entries;
	const bool with_compaction;

	FileDisk * disk;
	std::unique_ptr<bitfield> bitfield_;
	std::unique_ptr<bitfield_index> index;
	StreamBuffer buffer, cbuf;
	uint64_t buffer_start_pos = 0;
	uint64_t disk_read_pos = 0;
	uint32_t next_block_num_entries = 0;


	int32_t inflateOne( uint8_t * from, uint8_t * to, uint8_t count ){
		ParkBits parkStart( from, (k*count + 7)/8, (k*count + 7)/8*8 );
		ParkBits parkEnd( from + (k*count+7)/8, ((k+kOffsetSize) * count + 7)/8, ((k+kOffsetSize)*count + 7)/8*8 );

		for( uint32_t i = 0; i < count; i++ ){
			Bits entry;
			entry.AppendValue( parkStart.SliceBitsToInt( i*k, i*k + k ), k );
			entry.AppendValue( parkEnd.SliceBitsToInt( i*(k+kOffsetSize), (i+1)*(k+kOffsetSize) ), k+kOffsetSize );
			entry.ToBytes( to + i * entry_size );
		}
		return (k*count + 7)/8 + ((k+kOffsetSize) * count + 7)/8;
	}

	void UncompactNext(){
		if( disk_read_pos == 0 ){
			disk->Read( 0, (uint8_t*)(&next_block_num_entries), 4 );
			disk_read_pos += 4;
		}

		// define size to read
		bool hasNext = buffer_start_pos / entry_size + next_block_num_entries < num_entries;
		uint32_t size_to_read = Util::ByteAlign( next_block_num_entries*k ) / 8
				+ Util::ByteAlign( next_block_num_entries*(k+kOffsetSize) )/8
				+ (hasNext?4:0);

		cbuf.ensureSize( size_to_read ).setUsed( size_to_read - (hasNext?4:0) );
		disk->Read( disk_read_pos, cbuf.get(), size_to_read );
		disk_read_pos += size_to_read;

		buffer.ensureSize( next_block_num_entries*entry_size ).setUsed( next_block_num_entries*entry_size );
		uint8_t* buf_at = cbuf.get();
		for( uint32_t i = 0; i < next_block_num_entries; i += 8 )
			buf_at += inflateOne( buf_at, buffer.get() + i * entry_size, std::min( 8U, next_block_num_entries-i ) );


		if( hasNext )
			memcpy( &next_block_num_entries, cbuf.get() + size_to_read - 4, 4 );
		else next_block_num_entries = 0;
	}

	void RewriteBuffer(){
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
};

#endif // SRC_CPP_LAST_TABLE_FILE_HPP_
