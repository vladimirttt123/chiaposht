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

#ifndef SRC_CPP_PHASE3_HPP_
#define SRC_CPP_PHASE3_HPP_

#include "encoding.hpp"
#include "entry_sizes.hpp"
#include "phase2.hpp"
#include "pos_constants.hpp"
#include "sort_manager.hpp"
#include "progress.hpp"
#include "park_writer.hpp"

#include <chrono>
#include <condition_variable>
using namespace std::chrono_literals; // for operator""min;

// Results of phase 3. These are passed into Phase 4, so the checkpoint tables
// can be properly built.
struct Phase3Results {
		// Pointers to each table start byet in the final file
		std::vector<uint64_t> final_table_begin_pointers;
		// Number of entries written for f7
		uint64_t final_entries_written;
		uint32_t right_entry_size_bits;

		uint32_t header_size;
		std::unique_ptr<SortManager> table7_sm;
};

struct NewPosTable1Reader{
	NewPosTable1Reader( const uint8_t k, FilteredDisk & disk ) : k(k), disk(disk) {	}

	// The left entries are in the new format: (sort_key, new_pos), except for table 1: (y, x).
	// Only k bits, since this is x
	inline uint64_t ReadNext() {	return Util::SliceInt64FromBytes(disk.ReadNext(), k);	}
private:
	const uint8_t k;
	FilteredDisk & disk;
};

struct NewPosSortManagerReader{
	NewPosSortManagerReader( const uint8_t k, const uint32_t right_sort_key_size, SortManager * sm)
			: k(k), right_sort_key_size(right_sort_key_size), sm(sm){}

	inline uint64_t ReadNext() {
		auto left_entry_disk_buf = sm->ReadEntry(left_reader);
		left_reader += sm->EntrySize();
		// k+1 bits in case it overflows ( seems this comment is incorrect )
		return Util::SliceInt64FromBytes(left_entry_disk_buf, right_sort_key_size, k);
	}
private:
	const uint8_t k;
	const uint32_t right_sort_key_size;
	SortManager * sm;
	uint64_t left_reader = 0;
};

struct P3Entry{
	uint64_t sort_key = 0, pos = 0, new_pos /*pos+offset*/= 0;

	inline void Set( uint8_t const* entry_buf, const uint8_t k ){
		sort_key = Util::SliceInt64FromBytes( entry_buf, /*sort_key_size*/k);
		pos = Util::SliceInt64FromBytes( entry_buf, /*sort_key_size*/k, /*pos_size*/k );
		new_pos = pos + /*offset*/ Util::SliceInt64FromBytes( entry_buf, /*sort_key_size*/k + /*pos_size*/k, kOffsetSize);
	}

	static inline uint64_t Pos( uint8_t const* entry_buf, const uint8_t k ){
		return Util::SliceInt64FromBytes( entry_buf, /*sort_key_size*/k, /*pos_size*/k );
	}
};

// This used to process entries from 2 tables in computational pass one of phase 3
struct PassOneBlockProcessor {
	const uint8_t k, left_entry_size_bytes, right_entry_size_bytes, right_sort_key_size, line_point_size, p2_entry_size_bytes;
	const uint64_t POSITION_LIMIT;

	PassOneBlockProcessor( const uint8_t k, bool is_table_1, const uint8_t left_entry_size_bytes,
												const uint8_t right_entry_size_bytes, const uint8_t p2_entry_size_bytes )
			: k(k), left_entry_size_bytes(left_entry_size_bytes), right_entry_size_bytes(right_entry_size_bytes)
			, right_sort_key_size(is_table_1?0:k), line_point_size( 2 * k - 1 ), p2_entry_size_bytes( p2_entry_size_bytes )
			, POSITION_LIMIT((uint64_t)1 << k)
	{}

	// returns processed size in bytes
	uint64_t Process( const uint64_t left_disk_start_idx,	const uint64_t left_table_end_idx, const uint8_t *left_disk,
					const uint64_t right_disk_size_bytes, const uint8_t * right_disk, SortManager::ThreadWriter &R_sort_manager) {

		assert( (right_disk_size_bytes%right_entry_size_bytes) == 0 );
		assert( P3Entry::Pos( right_disk, k ) >= left_disk_start_idx );


		uint8_t entry_buf[right_entry_size_bytes + 8/*safe distance*/];

		P3Entry entry;
		for( uint64_t current_pos = 0; current_pos < right_disk_size_bytes; current_pos += p2_entry_size_bytes ){
			// The right entries are in the format from backprop, (sort_key, pos, offset)
			entry.Set( right_disk + current_pos, k );


			if( entry.new_pos >= left_table_end_idx )
				return current_pos; // not enough left disk

			assert( entry.pos >= left_disk_start_idx );

			uint64_t left_new_pos_1 = Util::SliceInt64FromBytes( left_disk + (entry.pos - left_disk_start_idx) * left_entry_size_bytes, right_sort_key_size, k );
			uint64_t left_new_pos_2 = Util::SliceInt64FromBytes( left_disk + (entry.new_pos - left_disk_start_idx) * left_entry_size_bytes, right_sort_key_size, k );

			// A line point is an encoding of two k bit values into one 2k bit value.
			uint128_t line_point = Encoding::SquareToLinePoint( left_new_pos_1, left_new_pos_2 );

			if( left_new_pos_1 > POSITION_LIMIT || left_new_pos_2 > POSITION_LIMIT ) {
				std::cout << "left or right positions too large" << std::endl;
				std::cout << (line_point > ((uint128_t)1 << (2 * k)));
				if ((line_point > ((uint128_t)1 << (2 * k)))) {
					std::cout << "L, R: " << left_new_pos_1 << " " << left_new_pos_2
										<< std::endl;
					std::cout << "Line point: " << line_point << std::endl;
					abort();
				}
			}

			Bits to_write = Bits(line_point, line_point_size);
			to_write.AppendValue( entry.sort_key, /*right_sort_key_size*/k );

			to_write.ToBytes( entry_buf );
			R_sort_manager.Add( entry_buf );
		}

		return right_disk_size_bytes; // all done
	}
};

template<class T>
void PassOneRegular( const uint8_t k, bool is_table_1, T &left_disk,
										 const uint32_t p2_entry_size_bytes, const uint32_t right_sort_key_size,
										 SortManager &right_disk, SortManager &R_sort_manager, uint32_t num_threads ){

	assert( num_threads > 1 );
	//assert( num_threads = 1 ); // for debug

	left_disk.EnsureSortingStarted();
	right_disk.EnsureSortingStarted();

	assert( left_disk.CurrentBucketEnd() > kCachedPositionsSize*left_disk.EntrySize() ); // bucket should be bigger than size of cached position

	std::unique_ptr<SortManager::ThreadWriter> writers[num_threads];
	for( uint32_t i = 0; i < num_threads; i++ )
		writers[i].reset( new SortManager::ThreadWriter(R_sort_manager) );
	std::unique_ptr<std::thread> threads[num_threads];

	PassOneBlockProcessor block_processor( k, is_table_1, left_disk.EntrySize(), right_disk.EntrySize(), p2_entry_size_bytes );

	const uint64_t left_cache_size = kCachedPositionsSize*left_disk.EntrySize();
	const uint64_t chunk_size = std::max( BUF_SIZE/right_disk.EntrySize()*right_disk.EntrySize(),
																			 2UL*left_cache_size /* need to threads not overlap */ );

	std::atomic_uint64_t right_bucket_from_pos = 0, right_bucket_last_position = right_disk.CurrentBucketEnd();

	for( int64_t right_entries_todo = right_disk.Count(); right_entries_todo >= 0; ){ // for all right entries
		// start threads...
		for( uint32_t t = 0; t < num_threads; t++ )
			threads[t].reset( new std::thread([ &chunk_size, &right_bucket_from_pos, &left_disk, &right_disk,
																				&right_bucket_last_position, &block_processor, &k]( SortManager::ThreadWriter *wri ){
				const uint64_t left_bucket_start_idx = left_disk.CurrentBucketStart()/left_disk.EntrySize();
				const uint64_t left_bucket_end_idx = left_disk.CurrentBucketEnd()/left_disk.EntrySize();
				while(true){
					uint64_t right_bucket_start = right_bucket_from_pos.fetch_add( chunk_size, std::memory_order::relaxed );
					if( right_bucket_start >= right_disk.CurrentBucketSize() ) return; // right disk need switch

					uint64_t right_chunk_size_bytes = std::min( chunk_size, right_disk.CurrentBucketSize() - right_bucket_start );

					// Cannot process chunk with first entry position too close to end of the left bucket but not the last one.
					uint64_t processed_to = !left_disk.CurrentBucketIsLast()
																					&& (P3Entry::Pos( right_disk.CurrentBucketBuffer() + right_bucket_start, k ) + kCachedPositionsSize) >= left_bucket_end_idx ? 0 :
																			block_processor.Process( left_bucket_start_idx, left_bucket_end_idx, left_disk.CurrentBucketBuffer(),
																											right_chunk_size_bytes, right_disk.CurrentBucketBuffer() + right_bucket_start, *wri );

					assert( processed_to + right_bucket_start <= right_disk.CurrentBucketSize() );

					if( processed_to < right_chunk_size_bytes ){ // chunk is not finished need left bucket switch
						// set last processed postion outside thread
						uint64_t cur = right_bucket_last_position.load(std::memory_order::relaxed), top_value = processed_to + right_bucket_start;
						while( cur > top_value && !right_bucket_last_position.compare_exchange_weak( cur, top_value, std::memory_order::relaxed ) )/*empty loop*/;

						return;
					}
				}
			}, writers[t].get() ) );

		// wait for threads
		for( uint32_t t = 0; t < num_threads; t++ )
			threads[t]->join();

		auto switch_right_to_next_bucket = [&right_entries_todo, &right_disk, &right_bucket_from_pos, &right_bucket_last_position]( uint64_t pos ){
			if( pos != right_disk.CurrentBucketSize() ) return false;

			right_entries_todo -= right_disk.CurrentBucketCount();
			if( right_entries_todo > 0){
				right_disk.SwitchNextBucket(); // moving to next bucket
				right_bucket_from_pos = 0;
				right_bucket_last_position = right_disk.CurrentBucketSize();
			}
			return true;
		};

		// now need to finished not finished by threads.
		if( !switch_right_to_next_bucket( right_bucket_last_position.load(std::memory_order::relaxed) ) ){	// need to switch left disk to next bucket.
			uint8_t left_disk_cache[ 2UL*left_cache_size ];
			memcpy( left_disk_cache, left_disk.CurrentBucketBuffer() + left_disk.CurrentBucketSize() - left_cache_size, left_cache_size );
			left_disk.SwitchNextBucket(); // moving to next bucket
			// copy the rest
			uint32_t to_copy = std::min( left_cache_size, left_disk.CurrentBucketSize() );
			memcpy( left_disk_cache + left_cache_size, left_disk.CurrentBucketBuffer(), to_copy );

			right_bucket_from_pos = right_bucket_last_position.load(std::memory_order::relaxed);
			assert( right_bucket_from_pos < right_disk.CurrentBucketSize() );
			do{
				right_bucket_from_pos += block_processor.Process( left_disk.CurrentBucketStart()/left_disk.EntrySize() - kCachedPositionsSize,
																													(left_disk.CurrentBucketStart()+to_copy)/left_disk.EntrySize(), left_disk_cache,
																													right_disk.CurrentBucketSize() - right_bucket_from_pos,
																													right_disk.CurrentBucketBuffer() + right_bucket_from_pos, *writers[0].get() );

				switch_right_to_next_bucket( right_bucket_from_pos.load(std::memory_order::relaxed) );
			} while( P3Entry::Pos( right_disk.CurrentBucketBuffer() + right_bucket_from_pos, k ) < left_disk.CurrentBucketStart()/left_disk.EntrySize() );
			right_bucket_last_position = right_disk.CurrentBucketSize();
		}
	}
}


template<class T>
struct LeftTableCache{
	const uint8_t entry_size;
	const uint32_t cache_size;
	const uint64_t table_count;
	T &reader;
	uint64_t start_idx = 0, disk_idx = 0;
	std::unique_ptr<uint8_t[]> disk_cache;

	LeftTableCache( T & reader, uint8_t entry_size, uint64_t table_count )
			: entry_size(entry_size), cache_size(kCachedPositionsSize*entry_size), table_count(table_count),
									 reader(reader), disk_cache( new uint8_t(cache_size) ) {	}

	void ReadToIndexPlus( uint8_t *buf ){
		// start from copy cache
		uint64_t size = std::min( (uint64_t)kCachedPositionsSize, disk_idx)*entry_size;
		memcpy( buf, disk_cache.get(), size );
		buf += size;

		for( ;disk_idx < start_idx-1 && disk_idx < table_count; disk_idx++, buf += entry_size )
			memcpy( buf, reader.ReadNext(), entry_size );

		uint8_t *dc = disk_cache.get();
		for( uint64_t upto = disk_idx + std::min( (uint64_t)kCachedPositionsSize, table_count - disk_idx);
				 disk_idx < upto; disk_idx++, dc += entry_size, buf += entry_size ){
			const uint8_t * next = reader.ReadNext();
			memcpy( buf, next, entry_size );
			memcpy( dc, next, entry_size );
		}
	}
};

template<class T>
void PassOneSemiThreaded( const uint8_t k, const uint64_t left_table_count, T& left_disk,
													const uint32_t p2_entry_size_bytes, const uint32_t right_sort_key_size,
													const uint32_t right_entry_size_bytes, const uint64_t right_table_size,
													Disk& right_disk, SortManager * R_sort_manager, uint32_t num_threads ){

	PassOneBlockProcessor block_processor( k, left_disk.EntrySize(), right_entry_size_bytes, p2_entry_size_bytes );

	LeftTableCache left_cache( left_disk, left_disk.EntrySize(), left_table_count );
	uint64_t right_disk_pos = 0;

	std::unique_ptr<std::thread> threads[num_threads];
	std::mutex left_read_mutex, right_read_mutex;

	for( uint32_t i = 0; i < num_threads; i++ ){ // satrt threads
		threads[i].reset( new std::thread( [ &R_sort_manager, &left_disk, &left_cache, &right_read_mutex,
																				 &right_disk, &right_disk_pos, &k, &right_table_size, &left_read_mutex,
																				 &block_processor, &p2_entry_size_bytes ](){
			SortManager::ThreadWriter writer( *R_sort_manager );
			uint64_t left_max_entries = ( HUGE_MEM_PAGE_SIZE - MEM_SAFE_BUF_SIZE )/left_disk.EntrySize();
			auto left_buf = Util::allocate<uint8_t>( left_max_entries*left_disk.EntrySize() );
			uint64_t right_max_size = (HUGE_MEM_PAGE_SIZE-MEM_SAFE_BUF_SIZE)/ p2_entry_size_bytes*p2_entry_size_bytes;
			auto right_buf = Util::allocate<uint8_t>( right_max_size );
			uint64_t left_buf_start_idx, left_disk_end_idx;
			while( right_disk_pos < right_table_size ){
				uint64_t right_buf_size = 0;
				{
					std::lock_guard<std::mutex> m_guard( right_read_mutex );
					if( right_disk_pos >= right_table_size ) return; // we done
					left_buf_start_idx = left_cache.start_idx;

					// the max idx we can read up to
					uint64_t left_max_idx = left_cache.start_idx + left_max_entries - kCachedPositionsSize - 1;

					do{
						memcpy( right_buf.get() + right_buf_size, right_disk.Read( right_disk_pos, p2_entry_size_bytes ), p2_entry_size_bytes );
						left_disk_end_idx = P3Entry::Pos( right_buf.get() + right_buf_size, k );
						right_buf_size += p2_entry_size_bytes;
						right_disk_pos++;
					} while( right_buf_size < right_max_size && right_disk_pos < right_table_size && left_disk_end_idx < left_max_idx );

					// now we can update left disk cache data and free write lock after catching left lock
					left_cache.start_idx = left_disk_end_idx;
					left_read_mutex.lock();
				}

				// reading left disk
				left_cache.ReadToIndexPlus( left_buf.get() );
				left_read_mutex.unlock();

				auto tmp = block_processor.Process( left_buf_start_idx, left_disk_end_idx + kCachedPositionsSize, left_buf.get(), right_buf_size, right_buf.get(), writer );
				assert( tmp == right_buf_size );
			}
		} ) );
	}

	// wait for threads finished
	for( uint32_t i = 0; i < num_threads; i++ )
		threads[i]->join();
}

template<class T>
void PassOneSingleThreaded( const uint8_t k, const uint64_t left_table_count,
						 T& left_disk_new_pos_reader,
						 const uint32_t p2_entry_size_bytes, const uint32_t right_sort_key_size,
						 const uint32_t right_entry_size_bytes, const uint64_t right_table_size,
						 Disk& right_disk, SortManager * R_sort_manager ){

	const uint64_t POSITION_LIMIT = (uint64_t)1 << k;
	const uint8_t line_point_size = 2 * k - 1;
	StreamBuffer entry_buf(right_entry_size_bytes);

	uint64_t left_new_pos[kCachedPositionsSize];

	P3Entry entry;
	for( uint64_t current_pos = 0, left_entreis_read = 0, right_table_count = right_table_size/p2_entry_size_bytes;
			 current_pos < right_table_count; current_pos++ ){
		// The right entries are in the format from backprop, (sort_key, pos, offset)
		entry.Set( right_disk.Read( current_pos * p2_entry_size_bytes, p2_entry_size_bytes ), k );


		for( uint64_t left_max = std::min( entry.pos + kCachedPositionsSize, left_table_count ); left_max > left_entreis_read; left_entreis_read++ )
			left_new_pos[left_entreis_read%kCachedPositionsSize] = left_disk_new_pos_reader.ReadNext();

		uint64_t &left_new_pos_1 = left_new_pos[ entry.pos % kCachedPositionsSize ];
		uint64_t &left_new_pos_2 = left_new_pos[ entry.new_pos % kCachedPositionsSize ];

		// A line point is an encoding of two k bit values into one 2k bit value.
		uint128_t line_point = Encoding::SquareToLinePoint( left_new_pos_1, left_new_pos_2 );

		if( left_new_pos_1 > POSITION_LIMIT || left_new_pos_2 > POSITION_LIMIT ) {
			std::cout << "left or right positions too large" << std::endl;
			std::cout << (line_point > ((uint128_t)1 << (2 * k)));
			if ((line_point > ((uint128_t)1 << (2 * k)))) {
				std::cout << "L, R: " << left_new_pos_1 << " " << left_new_pos_2
									<< std::endl;
				std::cout << "Line point: " << line_point << std::endl;
				abort();
			}
		}

		Bits to_write = Bits(line_point, line_point_size);
		to_write.AppendValue( entry.sort_key, right_sort_key_size );

		to_write.ToBytes( entry_buf.get() );
		R_sort_manager->AddToCache( entry_buf.setUsed(right_entry_size_bytes) );
	}
}

// Compresses the plot file tables into the final file. In order to do this, entries must be
// reorganized from the (pos, offset) bucket sorting order, to a more free line_point sorting
// order. In (pos, offset ordering), we store two pointers two the previous table, (x, y) which
// are very close together, by storing  (x, y-x), or (pos, offset), which can be done in about k
// + 8 bits, since y is in the next bucket as x. In order to decrease this, We store the actual
// entries from the previous table (e1, e2), instead of pos, offset pointers, and sort the
// entire table by (e1,e2). Then, the deltas between each (e1, e2) can be stored, which require
// around k bits.

// Converting into this format requires a few passes and sorts on disk. It also assumes that the
// backpropagation step happened, so there will be no more dropped entries. See the design
// document for more details on the algorithm.
Phase3Results RunPhase3(
		uint8_t k,
		FileDisk &tmp2_disk /*filename*/,
		Phase2Results &res2,
		const uint8_t *id,
		const std::string &tmp_dirname,
		const std::string &filename,
		uint32_t header_size,
		MemoryManager &memory_manager,
		std::unique_ptr<SortStatisticsStorage> full_stats[],
		uint32_t num_buckets,
		uint32_t log_num_buckets,
		const uint8_t flags,
		uint32_t num_threads)
{
		uint8_t const line_point_size = 2 * k - 1;

		std::vector<uint64_t> final_table_begin_pointers(12, 0);
		final_table_begin_pointers[1] = header_size;

		uint8_t table_pointer_bytes[8];
		Util::IntToEightBytes(table_pointer_bytes, final_table_begin_pointers[1]);
		tmp2_disk.Write(header_size - 10 * 8, table_pointer_bytes, 8);

		uint32_t new_pos_entry_size_bytes = 0;

		uint32_t next_stats_idx = 0;
		std::unique_ptr<SortManager> L_sort_manager;
		std::unique_ptr<SortManager> R_sort_manager;

		// Sort key is k bits for all tables. For table 7 it is just y, which
		// is k bits, and for all other tables the number of entries does not
		// exceed 0.865 * 2^k on average.
		const uint32_t right_sort_key_size = k;

		const uint32_t p2_entry_size_bytes = EntrySizes::GetKeyPosOffsetSize(k);

		// Iterates through all tables, starting at 1, with L and R pointers.
		// For each table, R entries are rewritten with line points. Then, the right table is
		// sorted by line_point. After this, the right table entries are rewritten as (sort_key,
		// new_pos), where new_pos is the position in the table, where it's sorted by line_point,
		// and the line_points are written to disk to a final table. Finally, table_i is sorted by
		// sort_key. This allows us to compare to the next table.
		double progress_percent[] = {0.64, 0.7, 0.76, 0.82, 0.88, 0.94};
		for (int table_index = 1; table_index < 7; table_index++) {
				Timer table_timer;
				Timer computation_pass_1_timer;
				std::cout << "Compressing tables " << table_index << " and " << (table_index + 1) << std::endl;
				std::cout << "Progress update: " << std::setprecision(2) << progress_percent[table_index - 1] << std::endl;

				Disk& right_disk = res2.disk_for_table(table_index + 1);

				const uint32_t right_entry_size_bytes = EntrySizes::GetMaxEntrySize(k, table_index + 1, false);


				if (table_index > 1)
					L_sort_manager->FreeMemory();

				// We read only from this SortManager during the second pass, so all
				// memory is available
				R_sort_manager = std::make_unique<SortManager>(
						memory_manager,			*full_stats[next_stats_idx].get(),
						num_buckets,
						right_entry_size_bytes,
						tmp_dirname,				filename + ".p3.t" + std::to_string(table_index + 1),
						0/* begin_bits */,	0/* stripe_size*/,
						k/* plot size */,		3/* Phase */,
						table_index,				num_threads,
						(flags&NO_COMPACTION)==0, (flags&PARALLEL_READ)!=0 );
				next_stats_idx = !full_stats[1] ? 0 : (next_stats_idx + 1)%2;

				if( table_index == 1){
					if( num_threads <= 1 ){
						NewPosTable1Reader t1(k, res2.table1 );
						PassOneSingleThreaded( k, /*left table count:*/ res2.table_sizes[table_index], t1,
										p2_entry_size_bytes, right_sort_key_size,
										right_entry_size_bytes, /*right_table_size:*/ p2_entry_size_bytes * res2.table_sizes[table_index + 1],
										right_disk, R_sort_manager.get() );
					} else {
						PassOneRegular( k, true, res2.table1, p2_entry_size_bytes, right_sort_key_size,
													 (SortManager&)right_disk, *R_sort_manager.get(), num_threads );
					// 	PassOneSemiThreaded( k, /*left table count:*/ res2.table_sizes[table_index], res2.table1,
					// 											 p2_entry_size_bytes, right_sort_key_size,
					// 											 right_entry_size_bytes, /*right_table_size:*/ p2_entry_size_bytes * res2.table_sizes[table_index + 1],
					// 											 right_disk, R_sort_manager.get(), /*num_threads*/ 1 );
					}

					// Remove no longer needed file
					res2.table1.FreeMemory();

				}else{
					if( num_threads <= 1 || table_index == 6  ){
						NewPosSortManagerReader sr( k, right_sort_key_size, L_sort_manager.get() );
						PassOneSingleThreaded( k, /*left table count:*/ res2.table_sizes[table_index], sr,
										p2_entry_size_bytes, right_sort_key_size,
										right_entry_size_bytes, /*right_table_size:*/ p2_entry_size_bytes * res2.table_sizes[table_index + 1],
										right_disk, R_sort_manager.get() );
					} else {
						PassOneRegular( k, false, *L_sort_manager.get(), p2_entry_size_bytes, right_sort_key_size,
													 (SortManager&)right_disk, *R_sort_manager.get(), num_threads );
					}
				}

				// Flush cache so all entries are written to buckets
				R_sort_manager->FlushCache( !full_stats[1] );
				R_sort_manager->FreeMemory();

				computation_pass_1_timer.PrintElapsed("\tFirst computation pass time:");

				Timer computation_pass_2_timer;

				if (table_index > 1)
					L_sort_manager.reset(); // Make sure all files are removed


				// In the second pass we read from R sort manager and write to L sort
				// manager, and they both handle table (table_index + 1)'s data. The
				// newly written table consists of (sort_key, new_pos). Add one extra
				// bit for 'new_pos' to the 7-th table as it may have more than 2^k
				// entries.
				new_pos_entry_size_bytes = cdiv(2 * k + (table_index == 6 ? 1 : 0), 8);

				L_sort_manager = std::make_unique<SortManager>(
						memory_manager,		*full_stats[next_stats_idx].get(),
						num_buckets,
						new_pos_entry_size_bytes,
						tmp_dirname,					filename + ".p3s.t" + std::to_string(table_index + 1),
						0/* bits_begin */,		0/* strip_size */,
						k/* plot size */,			4/* Phase */,
						table_index + 1,			num_threads,
						(flags&NO_COMPACTION)==0, (flags&PARALLEL_READ)!=0 );
				next_stats_idx = !full_stats[1] ? 0 : (next_stats_idx + 1)%2;

				// Now we will write on of the final tables, since we have a table sorted by line point.
				// The final table will simply store the deltas between each line_point, in fixed space
				// groups(parks), with a checkpoint in each group.
				uint8_t const sort_key_shift = 128 - right_sort_key_size;
				uint8_t const index_shift = sort_key_shift - (k + (table_index == 6 ? 1 : 0));
				const uint32_t park_size_bytes = EntrySizes::CalculateParkSize( k, table_index );

				// At this point we know how many parks will be written,
				// than we can evalueate next table pointer
				const uint64_t entries_to_be_written = R_sort_manager->Count();
				const uint64_t R_sort_manager_size_bytes = entries_to_be_written*R_sort_manager->EntrySize();
				final_table_begin_pointers[table_index + 1] = final_table_begin_pointers[table_index] +
						( entries_to_be_written/kEntriesPerPark + ( (entries_to_be_written%kEntriesPerPark) ? 1 : 0) ) * park_size_bytes;
#define NO_MUTEX_ALG
#ifdef NO_MUTEX_ALG
				std::atomic_uint64_t reading_position = 0, processed_total = 0, bucket_end = 0;

				std::mutex write_mutex;
				auto parking_thread = [ &R_sort_manager, &line_point_size, &right_sort_key_size,
																&sort_key_shift, &index_shift, &L_sort_manager, &write_mutex,
																&tmp2_disk, &k, &final_table_begin_pointers, &table_index,
															 &park_size_bytes, &reading_position, &R_sort_manager_size_bytes, &processed_total, &bucket_end](){

					uint32_t sort_buf_size = kEntriesPerPark*R_sort_manager->EntrySize();
					uint8_t switch_bucket_buffer[sort_buf_size + MEM_SAFE_BUF_SIZE];
					ParkWriterTS parker( &tmp2_disk, &write_mutex,k, table_index,
															 final_table_begin_pointers[table_index], park_size_bytes );
					auto park_deltas = std::make_unique<std::vector<uint8_t>>();
					auto park_stubs = std::make_unique<std::vector<uint64_t>>();
					uint128_t last_line_point = 0;

					SortManager::ThreadWriter sort_writer = SortManager::ThreadWriter( *L_sort_manager );
					uint8_t const* right_reader_entry_buf;

					auto switch_to_next_bucket = [&processed_total, &R_sort_manager, &last_line_point, &line_point_size, &right_reader_entry_buf, &bucket_end]( uint64_t &pos ){
							last_line_point = Util::SliceInt128FromBytes( right_reader_entry_buf - R_sort_manager->EntrySize(), 0, line_point_size );
							uint64_t old;
							while( (old = processed_total.load(std::memory_order::relaxed) ) != pos )
								processed_total.wait( old, std::memory_order::relaxed );

							// no threads working with ram -> bucket can be switched
							R_sort_manager->SwitchNextBucket(); // throws exception if no additional buckets
							bucket_end.store( R_sort_manager->CurrentBucketEnd(), std::memory_order::relaxed );
							bucket_end.notify_all();
						};

					for( uint64_t position = reading_position.fetch_add( sort_buf_size, std::memory_order::relaxed);
							 position < R_sort_manager_size_bytes; position = reading_position.fetch_add( sort_buf_size, std::memory_order::relaxed) ) {

						sort_buf_size = std::min( (uint64_t)sort_buf_size, R_sort_manager_size_bytes-position ); // define real buffer size
						right_reader_entry_buf = R_sort_manager->CurrentBucketBuffer( std::min( R_sort_manager->CurrentBucketEnd(), position+sort_buf_size ) ) + position - R_sort_manager->CurrentBucketStart();

						if( ( position + sort_buf_size ) > /*>=?*/ R_sort_manager->CurrentBucketEnd() ){ // check the buffer inside current bucket
							if( position == R_sort_manager->CurrentBucketEnd() ){
								switch_to_next_bucket( position );
								right_reader_entry_buf = R_sort_manager->CurrentBucketBuffer( position + sort_buf_size );
							} else if( position < R_sort_manager->CurrentBucketEnd() ){
								uint64_t prev_size = R_sort_manager->CurrentBucketEnd() - position;
								assert( prev_size > 0 );
								assert( prev_size < sort_buf_size );
								memcpy( switch_bucket_buffer, right_reader_entry_buf, prev_size ); // store unprocessed prev bucket part
								switch_to_next_bucket( position );
								// add to processing bufer part from new bucket
								memcpy( switch_bucket_buffer + prev_size, R_sort_manager->CurrentBucketBuffer( position + sort_buf_size ), sort_buf_size - prev_size );
								right_reader_entry_buf = switch_bucket_buffer; // set processing buffer to cached
								assert( (position + sort_buf_size) <= R_sort_manager->CurrentBucketEnd() );
								assert( position < R_sort_manager->CurrentBucketStart() );
								assert( R_sort_manager->checkSort( right_reader_entry_buf, sort_buf_size ) );
							} else {
								uint64_t old;
								while(  position > ( old = bucket_end.load( std::memory_order::relaxed ) ) ) // wait for bucket switch
									bucket_end.wait( old, std::memory_order::relaxed ); // small sleep to wait

								assert( (position + sort_buf_size) <= R_sort_manager->CurrentBucketEnd() );
								assert( position > R_sort_manager->CurrentBucketStart() );
								right_reader_entry_buf = R_sort_manager->CurrentBucketBuffer( position + sort_buf_size ) + position - R_sort_manager->CurrentBucketStart();
							}
						}

						if( position > R_sort_manager->CurrentBucketStart() )
							last_line_point = Util::SliceInt128FromBytes( right_reader_entry_buf - R_sort_manager->EntrySize(), 0, line_point_size );

						// Every EPP entries, writes a park
						uint128_t checkpoint_line_point = Util::SliceInt128FromBytes(right_reader_entry_buf, 0, line_point_size);

						for( uint64_t index = position/R_sort_manager->EntrySize(), last_index = index + sort_buf_size/R_sort_manager->EntrySize()
								 ; index < last_index; index++, right_reader_entry_buf += R_sort_manager->EntrySize() ){

							// Right entry is read as (line_point, sort_key)
							uint128_t line_point = Util::SliceInt128FromBytes(right_reader_entry_buf, 0, line_point_size);
							uint64_t sort_key =
									Util::SliceInt64FromBytes(right_reader_entry_buf, line_point_size, right_sort_key_size);

							// Write the new position (index) and the sort key
							uint128_t to_write = (uint128_t)sort_key << sort_key_shift;
							to_write |= (uint128_t)index << index_shift;

							sort_writer.Add( to_write );

							uint128_t big_delta = line_point - last_line_point;

							// Since we have approx 2^k line_points between 0 and 2^2k, the average
							// space between them when sorted, is k bits. Much more efficient than storing each
							// line point. This is diveded into the stub and delta. The stub is the least
							// significant (k-kMinusStubs) bits, and largely random/incompressible. The small
							// delta is the rest, which can be efficiently encoded since it's usually very
							// small.

							uint64_t stub = big_delta & ((1ULL << (k - kStubMinusBits)) - 1);
							uint64_t small_delta = big_delta >> (k - kStubMinusBits);

							assert(small_delta < 256);

							if( (index % kEntriesPerPark) != 0 ) {
								park_deltas->push_back(small_delta);
								park_stubs->push_back(stub);
							}
							last_line_point = line_point;
						}

						processed_total.fetch_add( sort_buf_size, std::memory_order::relaxed );
						processed_total.notify_all();
						assert( processed_total <= R_sort_manager->CurrentBucketEnd() );
						parker.Write( position / R_sort_manager->EntrySize() / kEntriesPerPark, checkpoint_line_point, park_deltas, park_stubs );
					}
				};

				// Start threads
				R_sort_manager->EnsureSortingStarted();
				bucket_end.store( R_sort_manager->CurrentBucketEnd(), std::memory_order::relaxed );
#else // NO_MUTEX_ALG
				uint128_t global_last_line_point = 0;
				std::mutex read_mutex, write_mutex;
				auto parking_thread = [&right_entry_size_bytes, &R_sort_manager, &read_mutex,
															 &line_point_size, &right_sort_key_size, &global_last_line_point,
															 &sort_key_shift, &index_shift, &L_sort_manager, &write_mutex,
															 &tmp2_disk, &k, &final_table_begin_pointers, &table_index,
															 &park_size_bytes](){

					uint32_t sort_buf_size = kEntriesPerPark*right_entry_size_bytes;
					std::unique_ptr<uint8_t[]> sort_buf( Util::NewSafeBuffer( sort_buf_size ) );
					ParkWriterTS parker( &tmp2_disk, &write_mutex,k, table_index,
															final_table_begin_pointers[table_index], park_size_bytes );
					uint64_t index;
					auto park_deltas = std::make_unique<std::vector<uint8_t>>();
					auto park_stubs = std::make_unique<std::vector<uint64_t>>();
					uint128_t last_line_point = 0;

					SortManager::ThreadWriter sort_writer = SortManager::ThreadWriter( *L_sort_manager );

					while( true ){
						{
							std::lock_guard<std::mutex> lk(read_mutex);
							index = R_sort_manager->GetReadPosition()/right_entry_size_bytes;
							sort_buf_size = R_sort_manager->Read( sort_buf.get(), sort_buf_size );
							if( sort_buf_size == 0 ) return;
							last_line_point = global_last_line_point;
							global_last_line_point = Util::SliceInt128FromBytes(
									sort_buf.get() + sort_buf_size - right_entry_size_bytes, 0, line_point_size );
						}
						uint128_t checkpoint_line_point = 0;
						uint64_t park_index = index / kEntriesPerPark;

						uint8_t *right_reader_entry_buf = sort_buf.get();
						for( uint64_t last_index = index + sort_buf_size/right_entry_size_bytes
								 ; index < last_index; index++, right_reader_entry_buf += right_entry_size_bytes ){

							// Right entry is read as (line_point, sort_key)
							uint128_t line_point = Util::SliceInt128FromBytes(right_reader_entry_buf, 0, line_point_size);
							uint64_t sort_key =
									Util::SliceInt64FromBytes(right_reader_entry_buf, line_point_size, right_sort_key_size);

							// Write the new position (index) and the sort key
							uint128_t to_write = (uint128_t)sort_key << sort_key_shift;
							to_write |= (uint128_t)index << index_shift;

							sort_writer.Add( to_write );

							// Every EPP entries, writes a park
							if (index % kEntriesPerPark == 0)
								checkpoint_line_point = line_point;

							uint128_t big_delta = line_point - last_line_point;

							// Since we have approx 2^k line_points between 0 and 2^2k, the average
							// space between them when sorted, is k bits. Much more efficient than storing each
							// line point. This is diveded into the stub and delta. The stub is the least
							// significant (k-kMinusStubs) bits, and largely random/incompressible. The small
							// delta is the rest, which can be efficiently encoded since it's usually very
							// small.

							uint64_t stub = big_delta & ((1ULL << (k - kStubMinusBits)) - 1);
							uint64_t small_delta = big_delta >> (k - kStubMinusBits);

							assert(small_delta < 256);

							if( (index % kEntriesPerPark) != 0 ) {
								park_deltas->push_back(small_delta);
								park_stubs->push_back(stub);
							}
							last_line_point = line_point;
						}

						parker.Write( park_index, checkpoint_line_point, park_deltas, park_stubs );
					}
				};
#endif // NO_MUTEX_ALG
				if( num_threads <= 1 ){
					parking_thread();
				}
				else {
					uint32_t t_num = 1;// num_threads*2; // double the number of threads for this part

					std::unique_ptr<std::thread> threads[t_num];
					for( uint32_t t = 0; t < t_num; t++ )
						threads[t].reset( new std::thread(parking_thread) );

					for( uint32_t t = 0; t < t_num; t++ )
						threads[t]->join();
				}

				R_sort_manager.reset();
				L_sort_manager->FlushCache( table_index != 6 && !full_stats[1] );

				computation_pass_2_timer.PrintElapsed("\tSecond computation pass time:");
				std::cout << "\tWrote " << L_sort_manager->Count() << " entries" << std::endl;

				Encoding::ANSFree(kRValues[table_index - 1]);


				uint64_t final_table_writer = header_size - 8 * (10 - table_index);
				Util::IntToEightBytes(table_pointer_bytes, final_table_begin_pointers[table_index + 1]);
				tmp2_disk.Write(final_table_writer, (table_pointer_bytes), 8);

				table_timer.PrintElapsed("Total compress table time:");

				right_disk.FreeMemory();
				if (flags & SHOW_PROGRESS) { progress(3, table_index, 6); }
		}

		// no need free memory because it is going to be used stight away
		// L_sort_manager->FreeMemory();

		// These results will be used to write table P7 and the checkpoint tables in phase 4.
		return Phase3Results{
				final_table_begin_pointers,
				L_sort_manager->Count(),
				new_pos_entry_size_bytes * 8,
				header_size,
				std::move(L_sort_manager)};
}

#endif  // SRC_CPP_PHASE3_HPP
