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

#ifndef SRC_CPP_PHASE4_HPP_
#define SRC_CPP_PHASE4_HPP_

#include "disk.hpp"
#include "encoding.hpp"
#include "entry_sizes.hpp"
#include "phase3.hpp"
#include "pos_constants.hpp"
#include "util.hpp"
#include "progress.hpp"
#include "park_writer.hpp"


struct p4_syncer {
	p4_syncer( uint64_t first_bucket_end_pos ){
		bucket_end_pos = first_bucket_end_pos;
	}

	inline bool wait_for_bucket( const uint64_t pos, uint64_t &processed ){
		uint64_t b_end = bucket_end_pos.load(std::memory_order::relaxed);
		if( pos < b_end ) return false;

		processed_bytes.fetch_add( processed, std::memory_order::relaxed );
		processed_bytes.notify_all();
		processed = 0;

		do{
			bucket_end_pos.wait( b_end, std::memory_order::relaxed );
		}	while( pos >= (b_end = bucket_end_pos.load(std::memory_order::relaxed) ) ); // wait for bucket

		return true;
	}

	inline void wait_for_processed( uint64_t pos ){
		for( uint64_t old; pos > (old = processed_bytes.load(std::memory_order::relaxed)); )
			processed_bytes.wait( old, std::memory_order::relaxed );
	}

	inline void set_new_end( uint64_t pos ){
		bucket_end_pos.store( pos, std::memory_order::relaxed );
		bucket_end_pos.notify_all();
	}
private:
	std::atomic_uint64_t bucket_end_pos, processed_bytes = 0;
};

// Writes the checkpoint tables. The purpose of these tables, is to store a list of ~2^k values
// of size k (the proof of space outputs from table 7), in a way where they can be looked up for
// proofs, but also efficiently. To do this, we assume table 7 is sorted by f7, and we write the
// deltas between each f7 (which will be mostly 1s and 0s), with a variable encoding scheme
// (C3). Furthermore, we create C1 checkpoints along the way.  For example, every 10,000 f7
// entries, we can have a C1 checkpoint, and a C3 delta encoded entry with 10,000 deltas.

// Since we can't store all the checkpoints in
// memory for large plots, we create checkpoints for the checkpoints (C2), that are meant to be
// stored in memory during proving. For example, every 10,000 C1 entries, we can have a C2
// entry.

// The final table format for the checkpoints will be:
// C1 (checkpoint values)
// C2 (checkpoint values into)
// C3 (deltas of f7s between C1 checkpoints)
void RunPhase4(uint8_t k, uint8_t pos_size, FileDisk &tmp2_disk, Phase3Results &res,
							 const uint8_t flags, const int max_phase4_progress_updates, uint32_t num_threads)
{
		const uint32_t P7_park_size = Util::ByteAlign((k + 1) * kEntriesPerPark) / 8;
    uint64_t number_of_p7_parks =
				((res.final_entries_written == 0 ? 0 : res.final_entries_written - 1) / kEntriesPerPark) + 1;

		const uint64_t begin_byte_C1 = res.final_table_begin_pointers[7] + number_of_p7_parks * P7_park_size;

		const uint32_t C1_entry_buf_size_bytes = Util::ByteAlign(k) / 8;
		const uint64_t total_C1_entries = cdiv(res.final_entries_written, kCheckpoint1Interval);
		const uint64_t begin_byte_C2 = begin_byte_C1 + (total_C1_entries + 1) * C1_entry_buf_size_bytes;
		const uint64_t total_C2_entries = cdiv(total_C1_entries, kCheckpoint2Interval);
		const uint64_t begin_byte_C3 = begin_byte_C2 + (total_C2_entries + 1) * C1_entry_buf_size_bytes;

		const uint32_t size_C3 = EntrySizes::CalculateC3Size(k);
		const uint64_t end_byte = begin_byte_C3 + (total_C1_entries)*size_C3;

    res.final_table_begin_pointers[8] = begin_byte_C1;
    res.final_table_begin_pointers[9] = begin_byte_C2;
    res.final_table_begin_pointers[10] = begin_byte_C3;
    res.final_table_begin_pointers[11] = end_byte;

		uint64_t plot_file_reader = 0;
    uint64_t final_file_writer_1 = begin_byte_C1;


		StreamBuffer C2( ( cdiv( res.final_entries_written, kCheckpoint1Interval * kCheckpoint2Interval) + 1 ) * C1_entry_buf_size_bytes );
    uint64_t num_C1_entries = 0;
    uint32_t right_entry_size_bytes = res.right_entry_size_bits / 8;


		std::cout << "Progress update: 0.98" << std::endl;
    std::cout << "\tStarting to write C1 and C3 tables" << std::endl;
		res.table7_sm->EnsureSortingStarted();

		std::vector<uint8_t> deltas_to_write;
		if( num_threads <= 1 ){
			uint64_t prev_y = 0;
			uint64_t final_file_writer_2 = begin_byte_C3;
			uint64_t final_file_writer_3 = res.final_table_begin_pointers[7];

			uint8_t * P7_entry_buf = new uint8_t[P7_park_size];
			auto C1_entry_buf = new uint8_t[C1_entry_buf_size_bytes];
			auto C3_entry_buf = new uint8_t[size_C3];

			const int progress_update_increment = res.final_entries_written / max_phase4_progress_updates;
			ParkBits to_write_p7;

			// We read each table7 entry, which is sorted by f7, but we don't need f7 anymore. Instead,
			// we will just store pos6, and the deltas in table C3, and checkpoints in tables C1 and C2.
			for (uint64_t f7_position = 0; f7_position < res.final_entries_written; f7_position++) {
				auto right_entry_buf = res.table7_sm->ReadEntry(plot_file_reader);

				plot_file_reader += right_entry_size_bytes;
				uint64_t entry_y = Util::SliceInt64FromBytes(right_entry_buf, k);
				uint64_t entry_new_pos = Util::SliceInt64FromBytes(right_entry_buf, k, pos_size);


				if (f7_position % kEntriesPerPark == 0 && f7_position > 0) {
					memset(P7_entry_buf, 0, P7_park_size);
					to_write_p7.ToBytes(P7_entry_buf);
					tmp2_disk.Write(final_file_writer_3, (P7_entry_buf), P7_park_size);
					final_file_writer_3 += P7_park_size;
					to_write_p7 = ParkBits();
				}

				to_write_p7 += ParkBits(entry_new_pos, k + 1);

				if( f7_position % kCheckpoint1Interval == 0 ) {
					Bits entry_y_bits = Bits(entry_y, k);
					entry_y_bits.ToBytes(C1_entry_buf);
					tmp2_disk.Write(final_file_writer_1, C1_entry_buf, C1_entry_buf_size_bytes);
					final_file_writer_1 += C1_entry_buf_size_bytes;
					if (num_C1_entries > 0) {
						final_file_writer_2 = begin_byte_C3 + (num_C1_entries - 1) * size_C3;
						size_t num_bytes =
								Encoding::ANSEncodeDeltas(deltas_to_write, kC3R, C3_entry_buf + 2) + 2;

						// We need to be careful because deltas are variable sized, and they need to fit
						assert(size_C3 * 8/*why x8?*/ > num_bytes);

						// Write the size
						Util::IntToTwoBytes(C3_entry_buf, num_bytes - 2);

						tmp2_disk.Write(final_file_writer_2, (C3_entry_buf), num_bytes);
						final_file_writer_2 += num_bytes;
					}
					if (f7_position % (kCheckpoint1Interval * kCheckpoint2Interval) == 0) {
						assert( C2.size() > (C2.used()+C1_entry_buf_size_bytes) );
						C2.add( C1_entry_buf, C1_entry_buf_size_bytes );
					}
					deltas_to_write.clear();
					++num_C1_entries;
				} else {
					deltas_to_write.push_back(entry_y - prev_y);
				}
				prev_y = entry_y;

				if (flags & SHOW_PROGRESS && f7_position % progress_update_increment == 0) {
					progress(4, f7_position, res.final_entries_written);
				}
			}

			// Writes the final park to disk
			memset(P7_entry_buf, 0, P7_park_size);
			to_write_p7.ToBytes( P7_entry_buf );

			tmp2_disk.Write(final_file_writer_3, (P7_entry_buf), P7_park_size);
			// final_file_writer_3 += P7_park_size;

			if (!deltas_to_write.empty()) {
				size_t num_bytes = Encoding::ANSEncodeDeltas(deltas_to_write, kC3R, C3_entry_buf + 2);
				memset(C3_entry_buf + num_bytes + 2, 0, size_C3 - (num_bytes + 2));
				final_file_writer_2 = begin_byte_C3 + (num_C1_entries - 1) * size_C3;

				// Write the size
				Util::IntToTwoBytes(C3_entry_buf, num_bytes);

				tmp2_disk.Write(final_file_writer_2, (C3_entry_buf), size_C3);
				// final_file_writer_2 += size_C3;
				// Encoding::ANSFree(kC3R);
			}

			Bits( (uint128_t)0, Util::ByteAlign(k) ).ToBytes(C1_entry_buf);
			tmp2_disk.Write(final_file_writer_1, C1_entry_buf, C1_entry_buf_size_bytes);
			final_file_writer_1 += C1_entry_buf_size_bytes;
			assert( final_file_writer_1 == begin_byte_C1 + (cdiv( res.final_entries_written, kCheckpoint1Interval ) + 1)*C1_entry_buf_size_bytes );

			delete[] P7_entry_buf;
			delete[] C1_entry_buf;
			delete[] C3_entry_buf;
		} else { // num_threads > 1
			std::mutex write_mutex;
			p4_syncer thread_sync( res.table7_sm->CurrentBucketEnd() );

			ParkAggregator P7_parker( &tmp2_disk, res.final_table_begin_pointers[7], P7_park_size );
			const uint64_t final_size_written = res.final_entries_written*res.table7_sm->EntrySize();
			std::atomic_uint64_t P7_read_pos = 0;
			const uint32_t P7_read_chunk = kEntriesPerPark*right_entry_size_bytes;

			auto P7_thread_func = [&k, &pos_size, &res, &right_entry_size_bytes, &P7_parker, &P7_read_pos,
														 &final_size_written, &write_mutex, &thread_sync, &P7_read_chunk](){
				uint64_t processed = 0;

				while( true ){
					ParkBits to_write_p7;

					uint64_t start_to_read = P7_read_pos.fetch_add( P7_read_chunk, std::memory_order::relaxed );
					if( start_to_read >= final_size_written ) return; // thread is done

					const uint64_t park_idx = start_to_read / P7_read_chunk;
					uint8_t * P7_entry_buf = P7_parker.get_for( park_idx );

					thread_sync.wait_for_bucket( start_to_read, processed );

					uint64_t read_up_to= std::min( final_size_written, start_to_read + P7_read_chunk );
					const uint8_t * right_entry_buf = res.table7_sm->CurrentBucketBuffer( std::min( res.table7_sm->CurrentBucketEnd(), read_up_to ) )
																					 + start_to_read - res.table7_sm->CurrentBucketStart();

					for( ; start_to_read < read_up_to; start_to_read += right_entry_size_bytes,
																						 right_entry_buf += right_entry_size_bytes, processed += right_entry_size_bytes ){
						if( thread_sync.wait_for_bucket( start_to_read, processed ) )
							right_entry_buf = res.table7_sm->CurrentBucketBuffer( std::min( res.table7_sm->CurrentBucketEnd(), read_up_to ) );

						uint64_t entry_new_pos = Util::SliceInt64FromBytes( right_entry_buf, k, pos_size );
						to_write_p7 += ParkBits(entry_new_pos, k + 1);
					}
					memset( P7_entry_buf, 0, P7_parker.park_size_bytes );
					to_write_p7.ToBytes(P7_entry_buf);
					P7_parker.finish( park_idx, write_mutex );
				}
			};

			ParkAggregator C1_parker( &tmp2_disk, begin_byte_C1, C1_entry_buf_size_bytes );
			ParkAggregator C3_parker( &tmp2_disk, begin_byte_C3, size_C3 );
			std::atomic_uint64_t C_read_pos = 0;

			const uint32_t C_read_chunk = kCheckpoint1Interval*right_entry_size_bytes;
			auto C_thread_func = [&k, &res, &right_entry_size_bytes, &C1_parker, &C3_parker, &C_read_pos,
														&final_size_written, &C_read_chunk, &size_C3, &C2,
														&C1_entry_buf_size_bytes, &write_mutex, &thread_sync](){
				std::vector<uint8_t> deltas_to_write;
				uint64_t processed = 0;

				while( true ){
					uint64_t start_to_read = C_read_pos.fetch_add( C_read_chunk, std::memory_order::relaxed );
					if( start_to_read >= final_size_written ) return; // thread is done

					uint64_t park_idx = start_to_read / C_read_chunk;
					auto c3_buf = C3_parker.get_for( park_idx );
					auto c1_buf = C1_parker.get_for( park_idx );

					thread_sync.wait_for_bucket( start_to_read, processed );

					uint64_t read_up_to= std::min( final_size_written, start_to_read + C_read_chunk );
					const uint8_t * right_entry_buf = res.table7_sm->CurrentBucketBuffer( std::min( res.table7_sm->CurrentBucketEnd(), read_up_to ) )
																					 + start_to_read - res.table7_sm->CurrentBucketStart();

					uint64_t entry_y = Util::SliceInt64FromBytes(right_entry_buf, k);
					Bits entry_y_bits = Bits(entry_y, k);
					entry_y_bits.ToBytes( c1_buf );
					C1_parker.finish( park_idx, write_mutex );

					if( (park_idx%kCheckpoint2Interval) == 0 ){
						memcpy( C2.get() + (park_idx/kCheckpoint2Interval)*C1_entry_buf_size_bytes, c1_buf, C1_entry_buf_size_bytes );
						//C2.add( c1_buf, C1_entry_buf_size_bytes ); // this row can be problematic with 5000+ threadas :)
					}

					uint64_t prev_y = Util::SliceInt64FromBytes(right_entry_buf, k);
					for( start_to_read+=right_entry_size_bytes, right_entry_buf += right_entry_size_bytes, processed+=right_entry_size_bytes
							 ; start_to_read < read_up_to; start_to_read += right_entry_size_bytes,
																						 right_entry_buf += right_entry_size_bytes,
																						 processed += right_entry_size_bytes ){
						if( thread_sync.wait_for_bucket( start_to_read, processed ) )
							right_entry_buf = res.table7_sm->CurrentBucketBuffer( std::min( res.table7_sm->CurrentBucketEnd(), read_up_to ) );

						uint64_t entry_y = Util::SliceInt64FromBytes(right_entry_buf, k);

						deltas_to_write.push_back(entry_y - prev_y);
						prev_y = entry_y;
					};

					size_t num_bytes =
							Encoding::ANSEncodeDeltas(deltas_to_write, kC3R, c3_buf + 2) + 2;

					// We need to be careful because deltas are variable sized, and they need to fit
					assert(size_C3 * 8 > num_bytes);

					// Write the size
					Util::IntToTwoBytes( c3_buf, num_bytes - 2);
					memset( c3_buf + num_bytes, 0, C3_parker.park_size_bytes - num_bytes ); // clear park tail

					C3_parker.finish( park_idx, write_mutex );

					deltas_to_write.clear();
				};
			}; // end of C_thread_func

			// start threads
			num_threads = num_threads/2*2; // allign to even
			std::unique_ptr<std::thread> threads[num_threads];
			for( uint32_t t = 0; t < num_threads; t += 2 ){
				threads[t].reset( new std::thread(P7_thread_func) );
				threads[t+1].reset( new std::thread(C_thread_func) );
			}

			// loop to switch buckets
			while( !res.table7_sm->CurrentBucketIsLast() ){
				thread_sync.wait_for_processed( res.table7_sm->CurrentBucketEnd() * 2 );

				res.table7_sm->SwitchNextBucket();

				thread_sync.set_new_end( res.table7_sm->CurrentBucketEnd() );
			}

			for( uint32_t t = 0; t < num_threads; t++ )
				threads[t]->join();

			C2.setUsed( C2.size() - C1_entry_buf_size_bytes ); // add hoc to write correct C2
			final_file_writer_1 = begin_byte_C1 + (cdiv( res.final_entries_written, kCheckpoint1Interval ) + 1)*C1_entry_buf_size_bytes;
		}

    Encoding::ANSFree(kC3R);
    res.table7_sm.reset();

    std::cout << "\tFinished writing C1 and C3 tables" << std::endl;
    std::cout << "\tWriting C2 table" << std::endl;

		assert( C2.used()+C1_entry_buf_size_bytes == C2.size() );
		memset( C2.get() + C2.used(), 0, C1_entry_buf_size_bytes ); // last 0 entry?
		tmp2_disk.Write( final_file_writer_1, C2.get(), C2.size() );
    std::cout << "\tFinished writing C2 table" << std::endl;


    final_file_writer_1 = res.header_size - 8 * 3;
    uint8_t table_pointer_bytes[8];

    // Writes the pointers to the start of the tables, for proving
    for (int i = 8; i <= 10; i++) {
        Util::IntToEightBytes(table_pointer_bytes, res.final_table_begin_pointers[i]);
        tmp2_disk.Write(final_file_writer_1, table_pointer_bytes, 8);
        final_file_writer_1 += 8;
    }

    std::cout << "\tFinal table pointers:" << std::endl << std::hex;

    for (int i = 1; i <= 10; i++) {
        std::cout << "\t" << (i < 8 ? "P" : "C") << (i < 8 ? i : i - 7);
        std::cout << ": 0x" << res.final_table_begin_pointers[i] << std::endl;
    }
    std::cout << std::dec;
}
#endif  // SRC_CPP_PHASE4_HPP
