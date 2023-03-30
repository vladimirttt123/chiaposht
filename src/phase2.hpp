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
#include "phases.hpp"
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
		ReadStreamToDisk table7;
    std::vector<std::unique_ptr<SortManager>> output_files;
    std::vector<uint64_t> table_sizes;
};


inline void ScanTable( IReadDiskStream *disk, int table_index, const int64_t &table_size, int16_t const &entry_size,
											 bitfield &current_bitfield, bitfield &next_bitfield, const uint32_t &num_threads,
											 uint8_t const pos_offset_size, uint8_t const k ){
	Timer scan_timer;
	std::cout << "\ttable " << table_index << ": scan " << std::flush;

	next_bitfield.clear();

	// read_index is the number of entries we've processed so far (in the
	// current table) i.e. the index to the current entry. This is not used
	// for table 7

	int64_t read_cursor = 0;
	const auto max_threads = std::max((uint32_t)1, num_threads);
	auto threads = std::make_unique<std::thread[]>( max_threads );
	std::mutex read_mutex[2], union_mutex;
	// ensure buffer size is even.
	const int64_t read_bufsize = (BUF_SIZE/entry_size)*entry_size; // allign size to entry length
	// Run the threads
	for( uint32_t i = 0; i < max_threads; i++ ){
		threads[i] = std::thread( [table_index, entry_size, pos_offset_size, k, read_bufsize, &read_mutex, &union_mutex]
															(IReadDiskStream *disk, int64_t *read_cursor, const bitfield * current_bitfield, bitfield *next_bitfield){
			auto buffer = std::make_unique<uint8_t[]>(read_bufsize);
			bitfieldReader cur_bitfield( *current_bitfield );
			auto processed = std::make_unique<uint64_t[]>( read_bufsize/entry_size*2 );
			int64_t buf_size = 0, buf_start = 0;

			while( true ){
				{	// Read next buffer
					const std::lock_guard<std::mutex> lk(read_mutex[0]);
					buf_start = *read_cursor;
					buf_size = disk->Read( buffer, read_bufsize );
					if( buf_size == 0 ) return;// nothing to read -> exit

					*read_cursor += buf_size;
				}
				{ // Setting limits could read data from file than need a mutex
					const std::lock_guard<std::mutex> lk(read_mutex[1]);
					cur_bitfield.setLimits( buf_start/entry_size, buf_size/entry_size );
				}

				uint32_t processed_count = 0;

				// Convert buffer to numbers in final bitfield
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

					uint64_t entry_pos = entry_pos_offset >> kOffsetSize;
					uint64_t entry_offset = entry_pos_offset & ((1U << kOffsetSize) - 1);

					// mark the two matching entries as used (pos and pos+offset)
					processed[processed_count++] = entry_pos;
					processed[processed_count++] = entry_pos + entry_offset;
				}

				if( processed_count > 0 ){
					const std::lock_guard<std::mutex> lk(union_mutex);
					next_bitfield->set( processed.get(), processed_count );
				}
			}
		}, disk, &read_cursor, &current_bitfield, &next_bitfield );
	}

	// Wait for job done
	for( uint32_t i = 0; i < max_threads; i++ )
		threads[i].join();

	scan_timer.PrintElapsed( "time =" );
}

inline void SortRegularTableThread( IReadDiskStream * disk, const uint64_t &table_size,
																				const int16_t &entry_size, const int16_t &new_entry_size,
																				uint64_t *read_position, uint64_t *global_write_counter,
																				const bitfield *current_bitfield, const bitfield_index &index,
																				SortManager * sort_manager, std::mutex *sync_mutex,
																				const uint8_t &pos_offset_size, const uint8_t &k )
{
	uint8_t const write_counter_shift = 128 - k;
	uint8_t const pos_offset_shift = write_counter_shift - pos_offset_size;

	uint64_t buf_size = (BUF_SIZE/entry_size)*entry_size;
	auto buffer = std::make_unique<uint8_t[]>(buf_size);
	SortManager::ThreadWriter writer = SortManager::ThreadWriter( *sort_manager );

	bitfieldReader cur_bitfield = bitfieldReader( *current_bitfield );
	uint64_t write_counter = 0;

	while( true ){
		{	// Read next buffer
			const std::lock_guard<std::mutex> lk(*sync_mutex);

			buf_size = disk->Read( buffer, buf_size );
			if( buf_size == 0 ) return;// nothing to read -> exit

			cur_bitfield.setLimits( *read_position/entry_size, buf_size/entry_size );
			*read_position += buf_size;
			write_counter = *global_write_counter;
			*global_write_counter += cur_bitfield.count( 0, buf_size/entry_size );

//			assert( *global_write_counter == (uint64_t)current_bitfield->count(0 , *read_position/entry_size) );
		}

		for( uint64_t buf_ptr = 0; buf_ptr < buf_size; buf_ptr += entry_size ){

			// skipping
			if( !cur_bitfield.get( buf_ptr/entry_size ) ) continue;

			uint8_t const* entry = buffer.get() + buf_ptr;

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
			uint128_t new_entry = (uint128_t)write_counter << write_counter_shift;
			new_entry |= (uint128_t)entry_pos_offset << pos_offset_shift;

			writer.Add( new_entry );

			++write_counter;
		}
	}
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

				int64_t const table_size = table_sizes[table_index];
				int16_t const entry_size = cdiv(k + kOffsetSize + (table_index == 7 ? k : 0), 8);

				std::cout << "Backpropagating on table " << table_index << "  size: " << table_size <<  std::endl;
        std::cout << "Progress update: " << progress_percent[7 - table_index] << std::endl;

				// 2 instances of bitfield can be larger than memory_size
				// for example recomended minimum for k32 is 900MiB and
				// it is acound 2^32 entries, that implice  2^29 bytes of bitfield, or around 512MiB
				// that twice is 1GiB or larger than 900KiB
				// in such case we can operate with one bitfield only and the current one write
				// to disk and use FilteredDisk class
				if( table_index < 7 && std::max(bitfield::memSize(table_size), current_bitfield->memSize())*2 > memory_size ){
					std::cout << "Warning: Cannot fit 2 bitfields in buffer, flushing one to disk. Need Ram:"
										<< (std::max(bitfield::memSize(table_size), current_bitfield->memSize())*2.0/1024/1024/1024) << "GiB" << std::endl;
					current_bitfield->FlushToDisk( tmp_1_disks[table_index].GetFileName() + ".bitfield.tmp" );
				}
				auto next_bitfield = std::make_unique<bitfield>( std::max( current_bitfield->size(), table_size ) );

				{ // Scope for reader
					auto table_reader = std::unique_ptr<IReadDiskStream>(  table_index == 7 ?
										CreateLastTableReader( &tmp_1_disks[table_index], k, entry_size, (flags&NO_COMPACTION)==0 ):
										//new AsyncStreamReader(
											new ReadFileStream( &tmp_1_disks[table_index], table_size * entry_size )
										//		, (BUF_SIZE/entry_size)*entry_size )
										);
					ScanTable( table_reader.get(), table_index, table_size, entry_size,
										 *current_bitfield.get(), *next_bitfield.get(), num_threads, pos_offset_size, k );
				}
				std::cout << "\ttable " << table_index << ": sort " << std::flush;
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

				if( table_index == 7 ){
					// Instead of rewritting table 7 we store the bitfield and read it rewtited in phase 3
					next_bitfield->FlushToDisk( tmp_1_disks[table_index].GetFileName() + ".bitfield.tmp" );
				}
				else{
					bitfield_index const index(*next_bitfield.get());

					auto sort_manager = std::make_unique<SortManager>(
							table_index == 2 ? ( memory_size /*it should be minus size of bitfield of table 1*/): memory_size / 2,
							num_buckets,
							log_num_buckets,
							new_entry_size,
							tmp_dirname,
							filename + ".p2.t" + std::to_string(table_index),
							uint32_t(k), // bits_begin
							0, // strip_size
							k,
							2, // Phase
							table_index,
							num_threads,
							(flags&NO_COMPACTION)==0 );

					uint64_t read_position = 0, write_counter = 0;
					std::mutex sort_mutext;
					auto threads = std::make_unique<std::thread[]>( num_threads - 1 );

					{	// scope for reader stream
						auto table_stream = AsyncStreamReader(
									new ReadFileStream( &tmp_1_disks[table_index], table_size*entry_size ),
									(BUF_SIZE/entry_size)*entry_size );
						for( uint64_t t = 0; t < num_threads - 1; t++ )
							threads[t] = std::thread(	SortRegularTableThread, &table_stream,
																				table_size, entry_size, new_entry_size,
																				&read_position, &write_counter, current_bitfield.get(),
																				index, sort_manager.get(), &sort_mutext, pos_offset_size, k );

						SortRegularTableThread( &table_stream, table_size, entry_size, new_entry_size,
																			 &read_position, &write_counter, current_bitfield.get(),
																			 index, sort_manager.get(), &sort_mutext, pos_offset_size, k );

						for( uint64_t t = 0; t < num_threads - 1; t++ )
							threads[t].join();
					}

					assert( (uint64_t)current_bitfield->count(0, table_size) == sort_manager->Count() );
					assert( write_counter == sort_manager->Count() );

					std::cout << write_counter << " entries ";

					// clear disk caches and memory
					// TODO real flush by flag.
					// Table 2 is used immediatly after in phase 3 than do not clean it caches fully.
					sort_manager->FlushCache( table_index != 2 );  // close all files
					//sort_manager->FreeMemory(); // should do nothing at this point

					output_files[table_index - 2] = std::move(sort_manager);
					new_table_sizes[table_index] = write_counter;
			}

			sort_timer.PrintElapsed( ", time =" );
			current_bitfield.swap(next_bitfield);
			next_bitfield->FreeMemory( table_index != 6 );// do not delete file of table 7 bitfield

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

    std::cout << "table " << table_index << " new size: " << new_table_sizes[table_index] << std::endl;

		// to free memory we flush the bitfield to disk
		//current_bitfield->flush_to_disk( tmp_1_disks[table_index].GetFileName() + ".bitfield.tmp" );

		// TODO some more check, it can be memory leaks with current_bitfield
		BufferedDisk disk_table1(&tmp_1_disks[1], table_size * entry_size);
		return {
				FilteredDisk(std::move(disk_table1), current_bitfield.release(), entry_size)
				, ReadStreamToDisk( CreateLastTableReader( &tmp_1_disks[7], k, new_entry_size,
														new_table_sizes[7], tmp_1_disks[7].GetFileName() + ".bitfield.tmp",
														(flags&NO_COMPACTION)==0 ), new_entry_size )
        , std::move(output_files)
        , std::move(new_table_sizes)
    };
}

#endif  // SRC_CPP_PHASE2_HPP
