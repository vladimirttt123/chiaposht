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
#include "threading.hpp"

#include <mutex>
#include <atomic>

using namespace std::chrono_literals; // for time units;


struct LastTableWriter {
	LastTableWriter( FileDisk * file, uint8_t k, uint16_t entry_size,
									 bool withCompaction, bool withScan, bool full_scan )
		: k(k), entry_size(entry_size),	pos_offset_size( k + kOffsetSize )
		, with_scan(withScan), with_compaction(withCompaction), disk(file)
	{
		if( withScan && full_scan )
			bitmap.reset(  new bitfield( ((uint64_t)2)<<k /* POSSIBLE SIZE PROBLEM */ ) );

		//debug_disk = new FileDisk( file->GetFileName() + ".src.tmp" );
	}

	const inline uint64_t getWritePosition() { return write_position; }

	void Write( const uint64_t &write_pos, uint8_t * buf, const uint32_t &buf_size ){
		if( with_scan ) scan( buf, buf_size );

		if( with_compaction ){
//			{
//				std::lock_guard<std::mutex> lk(file_sync);
//				debug_disk->Write( write_pos, buf, buf_size ); // just write on mentioned position
//			}

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
		//debug_disk->Close();

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
//	FileDisk * debug_disk;
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

		assert( k + kOffsetSize < 64 ); // this works up to this limit only since kOffsetSize == 10 than up to k54
		// compact here
		ParkBits parkStart, parkMid, parkEnd;
		const uint32_t num_entries = buf_size / entry_size;

		// pre evaluation for header
		uint64_t last_val = Util::SliceInt64FromBytes( buf, k, k + kOffsetSize );
		int64_t max_diff = 0, min_diff = 0;
		for( uint32_t buf_ptr = entry_size; buf_ptr < buf_size; buf_ptr += entry_size ){
			uint64_t cur = Util::SliceInt64FromBytes( buf + buf_ptr, k, k + kOffsetSize );
			int64_t diff = cur - last_val;
			if( diff > max_diff ) max_diff = diff;
			if( diff < min_diff ) min_diff = diff;
			last_val = cur;
		}

		const uint8_t bits_need = std::log2( max_diff - min_diff ) + 1;
		last_val = Util::SliceInt64FromBytes( buf, k, k + kOffsetSize );
		min_diff = -min_diff;
		assert( std::log2( min_diff ) + 1 < bits_need );

		const uint8_t first_val_bytes = (k + kOffsetSize + 7 )/8;
		const uint32_t headers_size = 4 /*num_entries*/ + 1 /*bits_need*/
				+ first_val_bytes + (bits_need+7)/8 ;
		const uint32_t compacted_size = headers_size - 4 // variable header
				+ (k + bits_need)*(num_entries/8) // size of all full parks
				+ (k*(num_entries%8) + 7)/8 // not full start park
				+ (bits_need*(num_entries%8) + 7)/8; // not full endpark

		// write header
		cbuf.ensureSize( compacted_size + 4 ).setUsed( headers_size );
		memcpy( cbuf.get(), &num_entries, 4 ); // number of entry
		cbuf.get()[4] = bits_need;
		for( uint64_t i = 0; i < first_val_bytes; i++ )
			cbuf.get()[i + 5] = (last_val>>(i*8))&0xff;
		for( uint64_t i = 0; i < (bits_need+7U)/8U; i++ )
			cbuf.get()[i+5+first_val_bytes] = (min_diff>>(i*8))&0xff;


		for( uint32_t buf_ptr = 0; buf_ptr < buf_size; buf_ptr += entry_size ){
			if( buf_ptr > 0 && ((buf_ptr/entry_size)%8) == 0 ){
				addBits( cbuf, parkStart );
				addBits( cbuf, parkEnd );
				parkStart.Clear();
				parkEnd.Clear();
			}

			uint64_t start = Util::SliceInt64FromBytes( buf + buf_ptr, k );
			uint64_t end = Util::SliceInt64FromBytes( buf + buf_ptr, k, k + kOffsetSize );
			uint64_t end_diff = end - last_val + min_diff;
			last_val = end;
			assert( end_diff < (1UL<<bits_need) );

			parkStart.AppendValue( start, k );
			parkEnd.AppendValue( end_diff, bits_need );
		}

		if( parkStart.GetSize() > 0 ){
			addBits( cbuf, parkStart );
			addBits( cbuf, parkEnd );
		}

		assert( compacted_size + 4 == cbuf.used() );
	}
};


//=======================================================
//======================= READER ========================
//=======================================================
struct LastTableReader : Disk {

	LastTableReader( FileDisk * file, const uint8_t &k, const uint16_t &entry_size,
									 const uint64_t &num_entries,	bool withCompaction, uint16_t num_threads )
		: k(k), entry_size(entry_size)
		, f7_shift( 128 - k ), pos_offset_size( k + kOffsetSize )
		, t7_pos_offset_shift( f7_shift - pos_offset_size )
		, num_entries(num_entries), with_compaction(withCompaction)
		, max_threads(std::max( 1U, (uint32_t)num_threads )), disk(file)
		, threads( new ThreadData[max_threads] )
	{	}

	inline uint16_t EntrySize() const { return entry_size; }
	inline uint64_t Count() const { return num_entries; }

	uint8_t const* Read( uint64_t begin, uint64_t length) override{
		if( !bitfield_ ){ // first read
			bitfield_.reset( new bitfield( num_entries, disk->GetFileName() + ".bitfield.tmp" ) );
			index.reset( new bitfield_index( *bitfield_.get() ) );
			StartThread(); // start first thread
		}

		assert( length == entry_size ); // alway reads by one entry
		assert( (begin%entry_size) == 0 ); // all reads alligned to entry size
		assert( begin >= cur_buffer_start_pos );// forward only read
		assert( begin < num_entries*entry_size ); // read inside table
		assert( cur_buffer.used() > 0 ? (begin > 0) : (cur_buffer_start_pos == 0 && begin == 0 ) ); // zeros on begining of read only

		while( cur_buffer_start_pos + cur_buffer.used() <= begin ){
			assert( cur_buffer_start_pos + cur_buffer.used() == begin ); // it should read every entry

			// need to get next from threads...
			for( uint16_t i = 0; i < max_threads && threads[i].thread != nullptr; i++ ){
				if( threads[i].buffer_start_position == begin ) {
					if( !threads[i].ready && begin > 0 ) StartThread(); // start more threads if we need to wait for not first buffer
					Sem::Wait( &threads[i].sem_done ); // now get semaphore
					threads[i].buffer.swap( cur_buffer.setUsed( 0 ) ); // get buffer from thread
					cur_buffer_start_pos = threads[i].buffer_start_position;
					threads[i].ready = false;
					Sem::Post( &threads[i].sem_run ); // we got value from thread free it to run
				}
			}
			if( cur_buffer_start_pos + cur_buffer.used() <= begin ) { // we not found thread working on next
				if( begin > 0 ) StartThread();
				// it can happen when first thread havn't started yet and need to wait
				std::this_thread::sleep_for( 10ns );
			}
		}

		return cur_buffer.get() + begin - cur_buffer_start_pos;
	};


	void Write(uint64_t begin, const uint8_t *memcache, uint64_t length) override {
		throw InvalidStateException( "Write impossible to read only disk" );
	};

	void Truncate(uint64_t new_size) override {
		throw InvalidStateException( "Truncate impossible to read only disk" );
	};

	std::string GetFileName() override { return disk->GetFileName(); };
	void FreeMemory() override { cur_buffer.reset(); index.reset(); bitfield_.reset(); };

	~LastTableReader(){

	}
private: // *********************

	struct ThreadData {
		Sem::type sem_done, sem_run;
		std::thread * thread = nullptr;

		StreamBuffer buffer, cbuf;
		uint64_t buffer_start_position = 0;
		bool ready = false;

		ThreadData(){
			sem_run = Sem::Create();
			sem_done = Sem::Create();
			Sem::Post(&sem_run); // ready to start
		}
		~ThreadData(){
			if( thread != nullptr ) {
				thread->join();
				delete thread;
			}
			Sem::Destroy( sem_run );
			Sem::Destroy( sem_done );
		}
	};

	const uint8_t k;
	const uint16_t entry_size;
	const uint8_t f7_shift;
	const uint8_t pos_offset_size;
	const uint8_t t7_pos_offset_shift;
	const uint64_t num_entries;
	const bool with_compaction;
	const uint16_t max_threads;

	FileDisk * disk;
	std::mutex file_sync;
//	FileDisk * debug_disk = nullptr;
	std::unique_ptr<bitfield> bitfield_;
	std::unique_ptr<bitfield_index> index;
	StreamBuffer cur_buffer;
	uint64_t cur_buffer_start_pos = 0;
	uint64_t disk_read_pos = 0;
	uint64_t disk_unpacked_pos = 0;
	uint8_t next_block_prefix[5];
	std::unique_ptr<ThreadData[]> threads;


	void StartThread(){
		for( uint16_t i = 0; i < max_threads; i++ ){
			if( threads[i].thread == nullptr ){
				if( with_compaction ){
					threads[i].thread = new std::thread( [this]( ThreadData *dat ){
							while( disk_unpacked_pos < num_entries*entry_size ){
								Sem::Wait( &dat->sem_run );
								if( !UncompactNext( dat->buffer, dat->buffer_start_position, dat->cbuf ) )
									return; // nothing more to read;
								RewriteBuffer( dat->buffer );
								dat->ready = true;
								Sem::Post( &dat->sem_done );
							}
						}, &threads[i] );
				} else {
					threads[i].thread = new std::thread( [this]( ThreadData * dat ){
							const uint64_t buf_size = BUF_SIZE/entry_size*entry_size;
							while( true ){
								Sem::Wait( &dat->sem_run );
								{	// Read next buffer
									std::lock_guard<std::mutex> lk(file_sync);
									if( disk_read_pos  >= num_entries*entry_size )
										return ; // nothing more to read

									dat->buffer_start_position = disk_read_pos;
									dat->buffer.ensureSize( buf_size ).setUsed( std::min( buf_size, num_entries*entry_size - disk_read_pos ) );

									assert( dat->buffer.used() > 0 ); // we should be here when it is something to read;

									disk->Read( disk_read_pos, dat->buffer.get(), dat->buffer.used() );
									disk_read_pos += dat->buffer.used();
								}

								RewriteBuffer( dat->buffer );
								dat->ready = true;
								Sem::Post( &dat->sem_done );
							}
						}, &threads[i] );
				}

				return; // thread started
			}
		}
	}

	inline int32_t inflateOne( uint8_t * from, uint8_t * to, uint8_t count,
														 const uint8_t &bits_need, const uint64_t &min_diff, uint64_t &last_val ){
		const uint32_t parkStartSize = (k*count + 7)/8;
		const uint32_t parkEndSize = (bits_need*count+7)/8;
		ParkBits parkStart( from, parkStartSize, parkStartSize*8 );
		ParkBits parkEnd( from + parkStartSize, parkEndSize, parkEndSize*8 );

		for( uint32_t i = 0; i < count; i++ ){
			Bits entry;
			entry.AppendValue( parkStart.SliceBitsToInt( i*k, i*k + k ), k );
			uint64_t diff = parkEnd.SliceBitsToInt( i*bits_need, (i+1)*bits_need );
			uint64_t end = last_val + diff - min_diff;
			entry.AppendValue( end, k+kOffsetSize );
			last_val = end;
			entry.ToBytes( to + i * entry_size );
		}
		return parkStartSize + parkEndSize;
	}

	bool UncompactNext( StreamBuffer &buffer, uint64_t &buffer_start_position, StreamBuffer &cbuf ){
		uint32_t block_num_entries;
		uint8_t bits_need;

		{ // Read from file
			std::lock_guard<std::mutex> lk(file_sync);
			if( disk_unpacked_pos >= num_entries*entry_size ) return false;

			buffer_start_position = disk_unpacked_pos; // reset buffer start position as soon as possible

			if( disk_read_pos == 0 ){
				disk->Read( 0, (uint8_t*)(&next_block_prefix), 5 );
				disk_read_pos += 5;
			}

			block_num_entries = *((uint32_t*)next_block_prefix);
			bits_need = next_block_prefix[4];

			// define size to read
			bool hasNext = disk_unpacked_pos / entry_size + block_num_entries < num_entries;
			uint32_t size_to_read = (bits_need+7)/8 // min diff size
														+ (k+kOffsetSize+7)/8 // first value
														+ (k + bits_need)*(block_num_entries/8) // size of all full parks
														+ (k*(block_num_entries%8) + 7)/8 // not full start park
														+ (bits_need*(block_num_entries%8) + 7)/8 // not full endpark
														+ (hasNext?5:0);

			cbuf.ensureSize( size_to_read ).setUsed( size_to_read - (hasNext?5:0) );
			disk->Read( disk_read_pos, cbuf.get(), size_to_read );
			disk_read_pos += size_to_read;
			disk_unpacked_pos += block_num_entries*entry_size;

			if( hasNext )
				memcpy( next_block_prefix, cbuf.get() + size_to_read - 5, 5 );
			else memset( next_block_prefix, 0, 5);
		}



		const uint8_t first_val_bytes = (k + kOffsetSize + 7 )/8; // length of first end value
		uint64_t min_diff = 0, last_val = 0;
		for( uint64_t i = 0; i < first_val_bytes; i++ )
			last_val |= ((uint64_t)cbuf.get()[i])<<(i*8);
		for( uint64_t i = 0; i < (bits_need+7U)/8U; i++ )
			min_diff |= ((uint64_t)cbuf.get()[i+first_val_bytes])<<(i*8);

		buffer.ensureSize( block_num_entries*entry_size ).setUsed( block_num_entries*entry_size );
		uint8_t* buf_at = cbuf.get() + first_val_bytes + (bits_need+7)/8;
		for( uint32_t i = 0; i < block_num_entries; i += 8 )
			buf_at += inflateOne( buf_at, buffer.get() + i * entry_size, std::min( 8U, block_num_entries-i ),
														bits_need, min_diff, last_val );

		assert( buf_at == cbuf.get() + cbuf.used() );
		return true;
	}

	void RewriteBuffer( StreamBuffer & buffer ){
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

// Temporal implementation to allow threads in Phase3 computational pass 1
// need to be done better in particular reading source file
struct LastTableBucketReader{
	LastTableReader &reader;
	const uint64_t bucket_size, file_size;

	LastTableBucketReader( LastTableReader&reader, uint64_t bucekt_size )
			:reader(reader), bucket_size((bucekt_size-MEM_SAFE_BUF_SIZE)/reader.EntrySize()*reader.EntrySize() )
			, file_size( reader.Count()*reader.EntrySize() ), bucket(Util::allocate<uint8_t>(this->bucket_size)){}

	inline uint64_t Count() const { return reader.Count(); }
	inline uint16_t EntrySize() const { return reader.EntrySize(); }
	inline uint64_t CurrentBucketSize() const { return bucket_used; }
	inline uint64_t CurrentBucketCount() const { return bucket_used/reader.EntrySize();}
	inline const uint8_t* CurrentBucketBuffer(uint64_t up_to_global_pos = 0){
		return bucket.get();
	}
	inline void EnsureSortingStarted(){ SwitchNextBucket(); }

	inline void SwitchNextBucket(){
		for( bucket_used = 0; bucket_used < bucket_size && file_read_position < file_size ;bucket_used+=reader.EntrySize(), file_read_position += reader.EntrySize() )
			memcpy( bucket.get() + bucket_used, reader.Read( file_read_position, reader.EntrySize() ), reader.EntrySize() );
	}

private:
	std::unique_ptr<uint8_t, Util::Deleter<uint8_t>>  bucket;
	uint64_t bucket_used = 0, file_read_position = 0;
};



struct TableFileWriter {
	const uint16_t entry_leftover_bits, entry_size_bytes, entry_size_full_bytes;

	TableFileWriter( FileDisk &disk, uint16_t entry_size_bits )
			: entry_leftover_bits(entry_size_bits%8), entry_size_bytes( (entry_size_bits+7)/8 ),
			entry_size_full_bytes( entry_size_bits/8 ), disk( disk )
	{	}

	inline void Write( uint64_t begin, const uint8_t *mem, uint64_t length ){

		if( entry_size_bytes == entry_size_full_bytes )
			disk.Write( begin, mem, length );
		else{
			if( begin != write_position )
				throw InvalidStateException("incorrect begin of write");

			write_position += length;

			uint8_t buf[length+8];
			uint32_t buf_position = 0;

			for( uint64_t i = 0; i < length; i += entry_size_bytes ){
				memcpy( buf + buf_position, mem+i, entry_size_full_bytes ); // copy full bytes
				buf_position += entry_size_full_bytes; // update buf position
				leftovers = (leftovers<<entry_leftover_bits) | (mem[i+entry_size_full_bytes]>>(8-entry_leftover_bits));
				if( leftovers_count == 7 ){
					buf_position += leftovers_to_buf( buf + buf_position );
					assert( leftovers == 0 );
				} else {
					leftovers_count++;
				}
			}

			disk.Write( file_write_position, buf, buf_position );
			file_write_position += buf_position;
		}
	}

	~TableFileWriter(){
		if( leftovers_count > 0 ) {
			file_write_position += entry_size_full_bytes*(8-leftovers_count); // allign writing to file
			leftovers <<= (8-leftovers_count)*entry_leftover_bits;

			uint8_t buf[8];
			disk.Write(file_write_position, buf, leftovers_to_buf( buf ) );
		}
	}

private:
	FileDisk &disk;
	uint8_t leftovers_count = 0;
	uint64_t write_position = 0, file_write_position = 0, leftovers = 0;

	inline uint64_t leftovers_to_buf( uint8_t*buf ){
		for( uint i = 0; i < entry_leftover_bits; i++, leftovers >>= 8 )
			buf[i] = leftovers&0xff;

		leftovers_count = 0;

		return entry_leftover_bits;
	}
};

struct TableFileReader{
	const uint16_t entry_leftover_bits, entry_size_bytes, entry_size_full_bytes;
	const uint64_t block_size, max_position;

	TableFileReader( FileDisk &disk, uint16_t entry_size_bits, uint64_t entries_no = 0 )
			: entry_leftover_bits(entry_size_bits%8), entry_size_bytes( (entry_size_bits+7)/8 ),
			entry_size_full_bytes( entry_size_bits/8 ), block_size( entry_size_bits ),
			max_position( entries_no * entry_size_bytes ), disk( disk )
	{	}

	inline void Read( uint64_t begin, uint8_t *buf, uint64_t length ){
		if( entry_size_bytes == entry_size_full_bytes ) {
			disk.Read( begin, buf, length );
			return;
		}

		assert( (begin%(entry_size_bytes*8)) == 0 );

		uint64_t first_entry = begin / entry_size_bytes;
		uint64_t first_entry_block_no = first_entry/8;
		uint64_t read_pos = first_entry_block_no*block_size;
		uint64_t read_blocks = (length / entry_size_bytes + 7 )/8;
		uint64_t read_blocks_size = read_blocks*block_size;
		uint8_t inner_buf[read_blocks_size];
		disk.Read( read_pos, inner_buf, read_blocks_size );

		for( uint64_t i = 0; i < read_blocks; i++ ){
			uint8_t * block_buf = inner_buf + i*block_size;

			uint64_t leftovers = 0;
			for( uint64_t l = 1; l <= entry_leftover_bits; l++ )
				leftovers = (leftovers<<8) | block_buf[block_size-l];
			uint8_t leftovers_bytes[8];
			for( int l = 7; l >= 0; l--, leftovers>>=entry_leftover_bits )
				leftovers_bytes[l] = (leftovers&0xff)<<(8-entry_leftover_bits);

			for( uint64_t buf_pos = i*8*entry_size_bytes, e = 0; e < 8 && buf_pos < length; e++, buf_pos+=entry_size_bytes ){
				memcpy( buf + buf_pos, block_buf + e*entry_size_full_bytes, entry_size_full_bytes );
				buf[buf_pos+entry_size_full_bytes] = leftovers_bytes[e];
			}
		}
	}
	inline uint32_t Read( uint8_t * buf, uint64_t length ){
		uint32_t real_length = std::min( length, max_position - read_position );
		if( real_length )
			Read( read_position, buf, real_length );
		read_position += real_length;
		return real_length;
	}

private:
	FileDisk &disk;
	uint64_t read_position = 0;
};

#endif // SRC_CPP_LAST_TABLE_FILE_HPP_
