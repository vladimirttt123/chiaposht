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

#ifndef SRC_CPP_PHASE2_HPP_
#define SRC_CPP_PHASE2_HPP_

#include <mutex>

#include "disk.hpp"
#include "entry_sizes.hpp"
#include "sort_manager.hpp"
#include "bitfield.hpp"
#include "bitfield_index.hpp"
#include "progress.hpp"
#include "filtereddisk.hpp"

struct Phase2Results
{
		Disk& disk_for_table(int const table_index)
    {
        if (table_index == 1) return table1;
        else if (table_index == 7) return table7;
        else return *output_files[table_index - 2];
    }
    FilteredDisk table1;
    BufferedDisk table7;
    std::vector<std::unique_ptr<SortManager>> output_files;
    std::vector<uint64_t> table_sizes;
};


inline void ScanTable( FileDisk* const disk, int table_index, const int64_t &table_size, int16_t const &entry_size,
											 bitfield &current_bitfield, bitfield &next_bitfield, const uint32_t &num_threads,
											 uint8_t const pos_offset_size, uint8_t const k ){
	Timer scan_timer;

	next_bitfield.clear();

	// read_index is the number of entries we've processed so far (in the
	// current table) i.e. the index to the current entry. This is not used
	// for table 7

	int64_t read_cursor = 0;
	auto max_threads = std::max((uint32_t)1, num_threads);
	auto threads = std::make_unique<std::thread[]>( max_threads );
	std::mutex read_mutex, union_mutex;
	const uint32_t processed_buf_size = 16384;
	// Run the threads
	for( uint32_t i = 0; i < max_threads; i++ ){
		threads[i] = std::thread( [table_index, table_size, entry_size, pos_offset_size, k]
															(FileDisk* disk, int64_t *read_cursor, std::mutex *read_mutex,
																std::mutex *union_mutex,  const bitfield * current_bitfield, bitfield *next_bitfield){
			// ensure buffer size is even.
			uint64_t processed[processed_buf_size + (processed_buf_size&1)];
			uint32_t processed_count = 0;
			const int64_t read_bufsize = (BUF_SIZE/entry_size)*entry_size; // allign size to entry length
			auto buffer = std::make_unique<uint8_t[]>(read_bufsize);
			bitfieldReader cur_bitfield( *current_bitfield );
			while( true ){
				int64_t buf_size = 0, buf_start = 0;
				{	// Read next buffer
					const std::lock_guard<std::mutex> lk(*read_mutex);
					buf_start = *read_cursor;
					buf_size = std::min( read_bufsize, (table_size*(int64_t)entry_size) - *read_cursor );
					if( buf_size > 0 )
						disk->Read( *read_cursor, buffer.get(), buf_size );
					*read_cursor += buf_size;
					cur_bitfield.setLimits( buf_start/entry_size, buf_size/entry_size );
				}
				// Nothing to read => end of work
				if( buf_size == 0 ){
					if( processed_count > 0 ){
						const std::lock_guard<std::mutex> lk(*union_mutex);
						for( uint32_t i = 0; i < processed_count; i++ )
							next_bitfield->set( processed[i] );
					}
					return;
				}

				// Do normal work
				for( int64_t buf_ptr = 0, entry_pos_offset = 0; buf_ptr < buf_size; buf_ptr += entry_size ){
					if (table_index == 7) {
							// table 7 is special, we never drop anything, so just build
							// next_bitfield
							entry_pos_offset = Util::SliceInt64FromBytes( buffer.get() + buf_ptr, k, pos_offset_size);
					} else {
							if( !cur_bitfield.get( buf_ptr/entry_size ) )
							{
									// This entry should be dropped.
									continue;
							}
							entry_pos_offset = Util::SliceInt64FromBytes( buffer.get() + buf_ptr, 0, pos_offset_size);
					}

					// Check if buffer is full or half full and we can get mutex to flush results into bitfield
					if( processed_count >= processed_buf_size
							|| (processed_count > (processed_buf_size>>1) && union_mutex->try_lock() ) ) {
						// Check we need the lock
						if( processed_count >= processed_buf_size ) union_mutex->lock();
						for( uint32_t i = 0; i < processed_count; i++ )
							next_bitfield->set( processed[i] );
						union_mutex->unlock();
						processed_count = 0;
					}

					uint64_t entry_pos = entry_pos_offset >> kOffsetSize;
					uint64_t entry_offset = entry_pos_offset & ((1U << kOffsetSize) - 1);
					// mark the two matching entries as used (pos and pos+offset)
					processed[processed_count++] = entry_pos;
					processed[processed_count++] = entry_pos + entry_offset;

				}
			}
		}, disk, &read_cursor, &read_mutex, &union_mutex, &current_bitfield, &next_bitfield );
	}

	// Wait for job done
	for( uint32_t i = 0; i < num_threads; i++ )
		threads[i].join();

	std::cout << "scanned table " << table_index << std::endl;
	scan_timer.PrintElapsed("scanned time = ");
}

inline void SortTable7Thread( const uint8_t *buffer, const uint64_t buf_size,
															int16_t const &entry_size, bitfield_index const *index,
															uint64_t from_ptr, uint64_t to_ptr, uint8_t *result,
															uint8_t const pos_offset_size, uint8_t const k ){

	uint8_t const f7_shift = 128 - k;
	uint8_t const t7_pos_offset_shift = f7_shift - pos_offset_size;

	for( uint64_t buf_ptr = from_ptr; buf_ptr < to_ptr; buf_ptr += entry_size ){

		uint8_t const* entry = buffer + buf_ptr;

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
		memcpy( result + buf_ptr, tmp, entry_size );
	}
}


inline void SortTable7( FileDisk* const disk, bitfield_index const &index,
												const int64_t &table_size, int16_t const &entry_size,
												uint32_t num_threads,
												uint8_t const pos_offset_size, uint8_t const k){

	num_threads = std::max( (uint32_t)1, num_threads );
	const uint64_t MIN_SIZE_PER_THREAD = 256; // NEED TO THINK ABOUT REAL VALUE

	auto threads = std::make_unique<std::thread[]>(num_threads);
	uint64_t buf_size = ((uint64_t)num_threads)*(BUF_SIZE/entry_size)*entry_size;
	auto reader = BufferedReader( disk, 0, buf_size , table_size*entry_size );
	auto result = std::make_unique<uint8_t[]>(buf_size + 16);

	while( (buf_size = reader.MoveNextBuffer() ) > 0 ){

		uint64_t from_ptr = 0, size_per_thread = (buf_size/num_threads/entry_size)*entry_size;

		// Run threads
		if( size_per_thread > MIN_SIZE_PER_THREAD )
			for( uint64_t t = 0; t < num_threads-1; t++, from_ptr += size_per_thread )
				threads[t] = std::thread( SortTable7Thread, reader.GetBuffer(), buf_size, entry_size, &index,
																	from_ptr, from_ptr+size_per_thread, result.get(), pos_offset_size, k );

		// Do some work in current thread
		SortTable7Thread( reader.GetBuffer(), buf_size, entry_size, &index, from_ptr, buf_size, result.get(),
											pos_offset_size, k );

		// Wait for threads
		if( size_per_thread > MIN_SIZE_PER_THREAD )
			for( uint32_t t = 0; t < num_threads-1; t++ )
				threads[t].join();

		disk->Write( reader.GetBufferStart(), result.get(), buf_size );
	}
}


inline int64_t SortRegularTableThread( const uint8_t *buffer, const uint64_t buf_size,
																			 int16_t const &entry_size, int16_t const &new_entry_size,
																			 int64_t start_write_counter, const bitfieldReader *current_bitfield,
																			 bitfield_index const index, uint32_t const log_num_buckets,
																			 uint64_t thread_no, uint64_t num_threads, uint8_t * const result,
																			 uint8_t const pos_offset_size, uint8_t const k
																			 ,SortManager * sort_manager
																			 ) {

	uint8_t const write_counter_shift = 128 - k;
	uint8_t const pos_offset_shift = write_counter_shift - pos_offset_size;

	int64_t write_counter = start_write_counter;
	for( uint64_t buf_ptr = 0; buf_ptr < buf_size; buf_ptr += entry_size ){

		// skipping
		if( !current_bitfield->get( buf_ptr/entry_size ) ) continue;
		// check is it for current thread.
		if( (write_counter%num_threads) == thread_no ){
			uint8_t const* entry = buffer + buf_ptr;

			uint64_t entry_pos_offset = Util::SliceInt64FromBytes( entry, 0, pos_offset_size );
			uint64_t entry_pos = entry_pos_offset >> kOffsetSize;
			uint64_t entry_offset = entry_pos_offset & ((1U << kOffsetSize) - 1);

			// assemble the new entry and write it to the sort manager

			// map the pos and offset to the new, compacted, positions and offsets
			std::tie(entry_pos, entry_offset) = index.lookup(entry_pos, entry_offset);
			entry_pos_offset = (entry_pos << kOffsetSize) | entry_offset;

			// The new entry is slightly different. Metadata is dropped, to
			// save space, and the counter of the entry is written (sort_key). We
			// use this instead of (y + pos + offset) since its smaller.
			uint8_t bytes[16];
			uint128_t new_entry = (uint128_t)write_counter << write_counter_shift;
			new_entry |= (uint128_t)entry_pos_offset << pos_offset_shift;
			Util::IntTo16Bytes(bytes, new_entry);

			memcpy( result + (write_counter-start_write_counter)*new_entry_size, bytes, new_entry_size );
	//		sort_manager->AddToCache( bytes );
		}
		++write_counter;
	}

	return write_counter;
}
// Backpropagate takes in as input, a file on which forward propagation has been done.
// The purpose of backpropagate is to eliminate any dead entries that don't contribute
// to final values in f7, to minimize disk usage. A sort on disk is applied to each table,
// so that they are sorted by position.
Phase2Results RunPhase2(
    std::vector<FileDisk> &tmp_1_disks,
    std::vector<uint64_t> table_sizes,
    uint8_t const k,
    const uint8_t *id,
    const std::string &tmp_dirname,
    const std::string &filename,
    uint64_t memory_size,
    uint32_t const num_buckets,
    uint32_t const log_num_buckets,
		uint8_t const flags,
		uint32_t num_threads )
{
		num_threads = std::max((uint32_t)1,num_threads);
    // After pruning each table will have 0.865 * 2^k or fewer entries on
    // average
    uint8_t const pos_size = k;
		uint8_t const pos_offset_size = pos_size + kOffsetSize;
    uint8_t const new_entry_size = EntrySizes::GetKeyPosOffsetSize(k);

    std::vector<uint64_t> new_table_sizes(8, 0);
    new_table_sizes[7] = table_sizes[7];

    // Iterates through each table, starting at 6 & 7. Each iteration, we scan
    // the current table twice. In the first scan, we:

    // 1. drop entries marked as false in the current bitfield (except table 7,
    //    where we don't drop anything, this is a special case)
    // 2. mark entries in the next_bitfield that non-dropped entries have
    //    references to

    // The second scan of the table, we update the positions and offsets to
    // reflect the entries that will be dropped in the next table.

    // At the end of the iteration, we transfer the next_bitfield to the current bitfield
    // to use it to prune the next table to scan.

		auto current_bitfield = std::make_unique<bitfield>(0);

    std::vector<std::unique_ptr<SortManager>> output_files;

    // table 1 and 7 are special. They are passed on as plain files on disk.
    // Only table 2-6 are passed on as SortManagers, to phase3
    output_files.resize(7 - 2);

    // note that we don't iterate over table_index=1. That table is special
    // since it contains different data. We'll do an extra scan of table 1 at
    // the end, just to compact it.
    double progress_percent[] = {0.43, 0.48, 0.51, 0.55, 0.58, 0.61};
    for (int table_index = 7; table_index > 1; --table_index) {

        std::cout << "Backpropagating on table " << table_index << std::endl;
        std::cout << "Progress update: " << progress_percent[7 - table_index] << std::endl;

				int64_t const table_size = table_sizes[table_index];
				int16_t const entry_size = cdiv(k + kOffsetSize + (table_index == 7 ? k : 0), 8);

				// 2 instances of bitfield can be larger than memory_size
				// for example recomended minimum for k32 is 900MiB and
				// it is acound 2^32 entries, that implice  2^29 bytes of bitfield, or around 512MiB
				// that twice is 1GiB or larger than 900KiB
				// in such case we can operate with one bitfield only and the current one write
				// to disk and use FilteredDisk class
				if( table_index < 7 && std::max(bitfield::memSize(table_size), current_bitfield->memSize())*2 > memory_size ){
					std::cout << "Warning: Cannot fit 2 bitfields in buffer, flushing one to disk" << std::endl;
					current_bitfield->flush_to_disk( tmp_1_disks[table_index].GetFileName() + ".bitfield.tmp" );
				}
				auto next_bitfield = std::make_unique<bitfield>( std::max( current_bitfield->size(), table_size ) );

				ScanTable( &tmp_1_disks[table_index], table_index, table_size, entry_size,
									 *current_bitfield.get(), *next_bitfield.get(), num_threads, pos_offset_size, k );

        std::cout << "sorting table " << table_index << std::endl;
        Timer sort_timer;

        // read the same table again. This time we'll output it to new files:
        // * add sort_key (just the index of the current entry)
        // * update (pos, offset) to remain valid after table_index-1 has been
        //   compacted.
        // * sort by pos
        //
        // As we have to sort two adjacent tables at the same time in phase 3,
        // we can use only a half of memory_size for SortManager. However,
        // table 1 is already sorted, so we can use all memory for sorting
        // table 2.

				// as we scan the table for the second time, we'll also need to remap
				// the positions and offsets based on the next_bitfield.
				bitfield_index const index(*next_bitfield.get());

				if( table_index == 7 )
					SortTable7( &tmp_1_disks[table_index], index, table_size, entry_size,
											num_threads, pos_offset_size, k );
				else{
					auto sort_manager = std::make_unique<SortManager>(
							table_index == 2 ? memory_size : memory_size / 2,
							num_buckets,
							log_num_buckets,
							new_entry_size,
							tmp_dirname,
							filename + ".p2.t" + std::to_string(table_index),
							uint32_t(k),
							0,
							strategy_t::quicksort_last,
							num_threads );

					uint64_t buf_size = num_threads*(BUF_SIZE/entry_size)*entry_size;
					int64_t write_counter = 0;
					BufferedReader reader( &tmp_1_disks[table_index], 0, buf_size, table_size*entry_size );
					auto threads = std::make_unique<std::thread[]>( num_threads - 1 );
					std::thread write_thread;
					bitfieldReader cur_bitfield = bitfieldReader( *current_bitfield );
					while( (buf_size = reader.MoveNextBuffer()) > 0 ){
						uint8_t * result = new uint8_t[buf_size/entry_size*new_entry_size];
						cur_bitfield.setLimits( reader.GetBufferStart()/entry_size, buf_size/entry_size );

						// Run threads
						for( uint64_t t = 0; t < num_threads - 1; t++ )
							threads[t] = std::thread(SortRegularTableThread, reader.GetBuffer(), buf_size, entry_size,
																			 new_entry_size, write_counter, &cur_bitfield, index,
																			 log_num_buckets, t, num_threads, result, pos_offset_size, k, sort_manager.get() );

						auto new_write_counter = SortRegularTableThread( reader.GetBuffer(), buf_size, entry_size,
																														 new_entry_size, write_counter, &cur_bitfield, index,
																														 log_num_buckets, num_threads-1, num_threads, result, pos_offset_size,
																														 k, sort_manager.get() );

						for( uint64_t t = 0; t < num_threads - 1; t++ )
							threads[t].join();

						// Write results to sort manager without wait in independent thread
						if( write_counter > 0 ) write_thread.join();
						write_thread = std::thread( [new_entry_size](SortManager* srt_mngr, uint8_t * res, int64_t num_entries ){
								for( int64_t i = 0; i < num_entries; i++ )
									srt_mngr->AddToCache( res + i*new_entry_size );
								delete [] res;
						}, sort_manager.get(), result, new_write_counter - write_counter );


						write_counter = new_write_counter;
					}

					if( write_counter > 0 ) write_thread.join();

					// clear disk caches and memory
					sort_manager->FlushCache();
					sort_manager->FreeMemory();

					output_files[table_index - 2] = std::move(sort_manager);
					new_table_sizes[table_index] = write_counter;
			}
			sort_timer.PrintElapsed("sort time = ");
			current_bitfield.swap(next_bitfield);
			next_bitfield->free_memory();

			// The files for Table 1 and 7 are re-used, overwritten and passed on to
			// the next phase. However, table 2 through 6 are all written to sort
			// managers that are passed on to the next phase. At this point, we have
			// to delete the input files for table 2-6 to save disk space.
			// This loop doesn't cover table 1, it's handled below with the
			// FilteredDisk wrapper.
			if (table_index != 7) {
					tmp_1_disks[table_index].Truncate(0);
			}
			if (flags & SHOW_PROGRESS) {
					progress(2, 8 - table_index, 6);
			}
    }

    // lazy-compact table 1 based on current_bitfield

    int const table_index = 1;
    int64_t const table_size = table_sizes[table_index];
    int16_t const entry_size = EntrySizes::GetMaxEntrySize(k, table_index, false);

    // at this point, table 1 still needs to be compacted, based on
    // current_bitfield. Instead of compacting it right now, defer it and read
    // from it as-if it was compacted. This saves one read and one write pass
		new_table_sizes[table_index] = current_bitfield->count(0, table_size);
    BufferedDisk disk(&tmp_1_disks[table_index], table_size * entry_size);

    std::cout << "table " << table_index << " new size: " << new_table_sizes[table_index] << std::endl;

		// to free memory we flush the bitfield to disk
		current_bitfield->flush_to_disk( tmp_1_disks[table_index].GetFileName() + ".bitfield.tmp" );

		// TODO some more check, it can be memory leadk with current_bitfield
		return {
				FilteredDisk(std::move(disk), std::move(*current_bitfield.release()), entry_size)
        , BufferedDisk(&tmp_1_disks[7], new_table_sizes[7] * new_entry_size)
        , std::move(output_files)
        , std::move(new_table_sizes)
    };
}

#endif  // SRC_CPP_PHASE2_HPP
