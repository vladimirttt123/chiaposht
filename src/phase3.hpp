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
#include "exceptions.hpp"
#include "phase2.hpp"
#include "pos_constants.hpp"
#include "sort_manager.hpp"
#include "progress.hpp"

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

// This writes a number of entries into a file, in the final, optimized format. The park
// contains a checkpoint value (which is a 2k bits line point), as well as EPP (entries per
// park) entries. These entries are each divided into stub and delta section. The stub bits are
// encoded as is, but the delta bits are optimized into a variable encoding scheme. Since we
// have many entries in each park, we can approximate how much space each park with take. Format
// is: [2k bits of first_line_point]  [EPP-1 stubs] [Deltas size] [EPP-1 deltas]....
// [first_line_point] ...
void PrepareParkBuffer(
		uint32_t park_size_bytes,
		uint128_t first_line_point,
		const std::vector<uint8_t> &park_deltas,
		const std::vector<uint64_t> &park_stubs,
		uint8_t k,
		uint8_t table_index,
		uint8_t *park_buffer,
		uint64_t const park_buffer_size)
{
		uint8_t *index = park_buffer;

		first_line_point <<= 128 - 2 * k;
		Util::IntTo16Bytes(index, first_line_point);
		index += EntrySizes::CalculateLinePointSize(k);

		// We use ParkBits instead of Bits since it allows storing more data
		ParkBits park_stubs_bits;
		for (uint64_t stub : park_stubs) {
				park_stubs_bits.AppendValue(stub, (k - kStubMinusBits));
		}
		uint32_t stubs_size = EntrySizes::CalculateStubsSize(k);
		uint32_t stubs_valid_size = cdiv(park_stubs_bits.GetSize(), 8);
		park_stubs_bits.ToBytes(index);
		memset(index + stubs_valid_size, 0, stubs_size - stubs_valid_size);
		index += stubs_size;

		// The stubs are random so they don't need encoding. But deltas are more likely to
		// be small, so we can compress them
		double R = kRValues[table_index - 1];
		uint8_t *deltas_start = index + 2;
		size_t deltas_size = Encoding::ANSEncodeDeltas(park_deltas, R, deltas_start);

		if (!deltas_size) {
				// Uncompressed
				deltas_size = park_deltas.size();
				Util::IntToTwoBytesLE(index, deltas_size | 0x8000);
				memcpy(deltas_start, park_deltas.data(), deltas_size);
		} else {
				// Compressed
				Util::IntToTwoBytesLE(index, deltas_size);
		}

		index += 2 + deltas_size;

		if ((uint32_t)(index - park_buffer) > park_buffer_size) {
				std::cout << "index-park_buffer " << index - park_buffer << " park_buffer_size "
									<< park_buffer_size << std::endl;
				throw InvalidStateException(
						"Overflowed park buffer, writing " + std::to_string(index - park_buffer) +
						" bytes. Space: " + std::to_string(park_buffer_size));
		}
		memset(index, 0x00, park_size_bytes - (index - park_buffer));
}

void WriteParkToFile(
		FileDisk &final_disk,
		uint64_t table_start,
		uint64_t park_index,
		uint32_t park_size_bytes,
		uint128_t first_line_point,
		const std::vector<uint8_t> &park_deltas,
		const std::vector<uint64_t> &park_stubs,
		uint8_t k,
		uint8_t table_index,
		uint8_t *park_buffer,
		uint64_t const park_buffer_size)
{
	PrepareParkBuffer( park_size_bytes, first_line_point, park_deltas,
										 park_stubs, k, table_index, park_buffer, park_buffer_size );
	// Parks are fixed size, so we know where to start writing. The deltas will not go over
	// into the next park.
	uint64_t writer = table_start + park_index * park_size_bytes;
	final_disk.Write(writer, (uint8_t *)park_buffer, park_size_bytes);
}

struct ParkWriterTS{

	ParkWriterTS( FileDisk *final_disk, std::mutex *write_mutex, const uint8_t &k,
								const uint8_t &table_index, const uint64_t &table_start, const uint32_t &park_size_bytes )
		: final_disk(final_disk), write_mutex(write_mutex), k(k)
		, park_buffer_size( EntrySizes::CalculateLinePointSize(k)
												+ EntrySizes::CalculateStubsSize(k) + 2
												+ EntrySizes::CalculateMaxDeltasSize(k, 1) )
		, park_buffer( new uint8_t[park_buffer_size] )
		, table_start(table_start)
		, table_index(table_index)
		, park_size_bytes(park_size_bytes )
	{}

	void Write( uint64_t park_index, uint128_t first_line_point,
					std::unique_ptr<std::vector<uint8_t>> &park_deltas,
					std::unique_ptr<std::vector<uint64_t>> &park_stubs ){

		PrepareParkBuffer( park_size_bytes,
										 first_line_point, *park_deltas.get(), *park_stubs.get(),
										 k, table_index,
										 park_buffer.get(), park_buffer_size );

		uint64_t writer = table_start + park_index * park_size_bytes;
		{
			std::lock_guard<std::mutex> lk(*write_mutex);
			final_disk->Write( writer, park_buffer.get(), park_size_bytes );
		}

		park_deltas->clear();
		park_stubs->clear();
	}

private:
	FileDisk *final_disk;
	std::mutex *write_mutex;
	const uint8_t k;
	uint64_t const park_buffer_size;
	std::unique_ptr<uint8_t[]> park_buffer;
	const uint64_t table_start;
	const uint8_t table_index;

	// The park size must be constant, for simplicity, but must be big enough to store EPP
	// entries. entry deltas are encoded with variable length, and thus there is no
	// guarantee that they won't override into the next park. It is only different (larger)
	// for table 1
	const uint32_t park_size_bytes;
};



struct OldData{
	uint64_t sort_keys[kMaxMatchesSingleEntry];
	uint64_t offsets[kMaxMatchesSingleEntry];
	uint16_t count = 0;

	inline void add( uint64_t sort_key, uint64_t offset ){
		sort_keys[count] = sort_key;
		offsets[count] = offset;
		count++;
	}
};


struct EntryAsynRewriter{
	EntryAsynRewriter( uint8_t k, SortManager *R_sort_manager, uint32_t num_threads, const uint32_t right_entry_size_bytes,
										 uint64_t (&left_new_pos)[kCachedPositionsSize], OldData (&old_data)[kReadMinusWrite] )
		: k(k), num_threads( std::max( 1U, num_threads ) ), POSITION_LIMIT((uint64_t)1 << k)
		, line_point_size(2 * k - 1), right_sort_key_size(k)
		, right_entry_size_bytes(right_entry_size_bytes), sm(R_sort_manager)
		, writer(*R_sort_manager), BATCH_SIZE( this->num_threads <= 1 ? 0 : (BUF_SIZE/sizeof(BatchEntry)) )
		, left_new_pos(left_new_pos), old_data(old_data) {

		if( this->num_threads > 1 ){
			next_batch.reset( new Batch( BATCH_SIZE ) );
			full_batch.store( new Batch( BATCH_SIZE ) );
		}
	}

	inline void next( const uint64_t &current_pos ){
		if( current_pos + 1 >= kReadMinusWrite ) {
			uint64_t const write_pointer_pos = current_pos - kReadMinusWrite + 1;
			uint64_t left_new_pos_1 = left_new_pos[write_pointer_pos % kCachedPositionsSize];
			OldData &old_datum = old_data[write_pointer_pos % kReadMinusWrite];

			if( num_threads <= 1 ){ // no threads
				for (uint32_t counter = 0; counter < old_datum.count; counter++ ){
					writeNext(	writer,	left_new_pos_1,
											left_new_pos[old_datum.offsets[counter] % kCachedPositionsSize],
											old_datum.sort_keys[counter] );
				}
			} else { // processing with thread
				for (uint32_t counter = 0; counter < old_datum.count; counter++ ){
					assert( next_batch->size < BATCH_SIZE );
					next_batch->add( left_new_pos_1,
								left_new_pos[old_datum.offsets[counter] % kCachedPositionsSize],
								old_datum.sort_keys[counter] );

					while( next_batch->size >= BATCH_SIZE ){ // queue is full
						next_batch.reset( full_batch.exchange( next_batch.release(), std::memory_order_relaxed ) );
						if( processing_threads.size() == 0
								|| ( next_batch->size > 0 && processing_threads.size() < num_threads ) ){
							processing_threads.emplace_back( [this](){threaded_porcessor();} );
						}
						if( next_batch->size > 0 ) // wait for threads will process
							std::this_thread::sleep_for( 5us );
					}
				}
			}
		}
	}

	~EntryAsynRewriter(){
		if( num_threads > 1 ){
			finished = true;

			for( uint32_t i = 0; i < next_batch->size; i++ )
				processBatchEntry( writer, next_batch->entries[i]);

			for( auto &t : processing_threads )
				t.join();

			delete full_batch.load( std::memory_order_relaxed );
		}
	}
private:
	const uint8_t k;
	const uint32_t num_threads;
	const uint64_t POSITION_LIMIT;
	const uint8_t  line_point_size;
	const uint32_t right_sort_key_size;
	const uint32_t right_entry_size_bytes;


	struct BatchEntry{
		uint64_t left_new_pos_1;
		uint64_t left_new_pos_2;
		uint64_t old_sort_key;
	};
	struct Batch{
		Batch( uint32_t batch_size ){ entries.reset( new BatchEntry[batch_size] ); }
		std::unique_ptr<BatchEntry[]>entries;
		uint32_t size = 0;
		inline void add( uint64_t pos_1, uint64_t pos_2, uint64_t key ){
			entries[size].left_new_pos_1 = pos_1;
			entries[size].left_new_pos_2 = pos_2;
			entries[size++].old_sort_key = key;
		}
	};

	SortManager * sm;
	SortManager::ThreadWriter writer;
	const uint32_t BATCH_SIZE;
	std::unique_ptr<Batch> next_batch;
	std::atomic<Batch *> full_batch;

	bool finished = false;
	std::vector<std::thread> processing_threads;

	uint64_t (&left_new_pos)[kCachedPositionsSize];
	OldData (&old_data)[kReadMinusWrite];

	// ================ FUNCS =================
	inline void writeNext( SortManager::ThreadWriter & tw, uint64_t left_new_pos_1,	uint64_t left_new_pos_2, uint64_t old_sort_key ){
		// A line point is an encoding of two k bit values into one 2k bit value.
		uint128_t line_point =
				Encoding::SquareToLinePoint( left_new_pos_1, left_new_pos_2 );

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
		to_write.AppendValue( old_sort_key, right_sort_key_size );

		uint8_t entry_buf[right_entry_size_bytes+8];
		to_write.ToBytes( entry_buf );
		tw.Add( entry_buf );
	}
	inline void processBatchEntry( SortManager::ThreadWriter & tw, BatchEntry &entry ){
		writeNext( tw, entry.left_new_pos_1, entry.left_new_pos_2, entry.old_sort_key );
	}

	void threaded_porcessor(){
		SortManager::ThreadWriter thread_writer(*sm);
		std::unique_ptr<Batch> batch( new Batch(BATCH_SIZE) );

		while( !finished || full_batch.load(std::memory_order_relaxed)->size > 0 ){
			batch.reset( full_batch.exchange( batch.release(), std::memory_order_relaxed ) );
			if( batch->size == 0 ){
				if( finished ) return;
				std::this_thread::sleep_for( 5us );
			}
			else {
				assert( batch->size == BATCH_SIZE ); // only full batches should be here
				for( uint32_t i = 0; i < batch->size; i++ )
					processBatchEntry( thread_writer, batch->entries[i] );
				batch->size = 0; // clear the batch
			}
		}
	}

};


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
		uint32_t num_buckets,
		uint32_t log_num_buckets,
		const uint8_t flags,
		uint32_t num_threads)
{
		uint8_t const pos_size = k;
		uint8_t const line_point_size = 2 * k - 1;

		std::vector<uint64_t> final_table_begin_pointers(12, 0);
		final_table_begin_pointers[1] = header_size;

		uint8_t table_pointer_bytes[8];
		Util::IntToEightBytes(table_pointer_bytes, final_table_begin_pointers[1]);
		tmp2_disk.Write(header_size - 10 * 8, table_pointer_bytes, 8);

		uint32_t new_pos_entry_size_bytes = 0;

		std::unique_ptr<SortManager> L_sort_manager;
		std::unique_ptr<SortManager> R_sort_manager;

		// These variables are used in the WriteParkToFile method. They are preallocatted here
		// to save time.
//		ParkWriter parker( &tmp2_disk, k );

		// Iterates through all tables, starting at 1, with L and R pointers.
		// For each table, R entries are rewritten with line points. Then, the right table is
		// sorted by line_point. After this, the right table entries are rewritten as (sort_key,
		// new_pos), where new_pos is the position in the table, where it's sorted by line_point,
		// and the line_points are written to disk to a final table. Finally, table_i is sorted by
		// sort_key. This allows us to compare to the next table.
		double progress_percent[] = {0.66, 0.73, 0.79, 0.85, 0.92, 0.98};
		for (int table_index = 1; table_index < 7; table_index++) {
				Timer table_timer;
				Timer computation_pass_1_timer;
				std::cout << "Compressing tables " << table_index << " and " << (table_index + 1)
									<< std::endl;
				std::cout << "Progress update: " << std::setprecision(2) << progress_percent[table_index - 1] << std::endl;

				Disk& right_disk = res2.disk_for_table(table_index + 1);

				// Sort key is k bits for all tables. For table 7 it is just y, which
				// is k bits, and for all other tables the number of entries does not
				// exceed 0.865 * 2^k on average.
				const uint32_t right_sort_key_size = k;

				const uint32_t p2_entry_size_bytes = EntrySizes::GetKeyPosOffsetSize(k);
				const uint32_t right_entry_size_bytes = EntrySizes::GetMaxEntrySize(k, table_index + 1, false);

				uint64_t left_reader = 0;
				uint64_t right_reader = 0;
				uint64_t right_reader_count = 0;

				if (table_index > 1) {
						L_sort_manager->FreeMemory();
				}

				// We read only from this SortManager during the second pass, so all
				// memory is available
				R_sort_manager = std::make_unique<SortManager>(
						memory_manager,
						num_buckets,
						log_num_buckets,
						right_entry_size_bytes,
						tmp_dirname,
						filename + ".p3.t" + std::to_string(table_index + 1),
						0, // begin_bits
						0, // stripe_size
						k, // plot size
						3, // Phase
						table_index,
						num_threads,
						(flags&NO_COMPACTION)==0 );
				StreamBuffer entry_buffer( right_entry_size_bytes );

				bool should_read_entry = true;
				uint64_t left_new_pos[kCachedPositionsSize];

				OldData old_data[kReadMinusWrite];

				bool end_of_right_table = false;

				uint64_t end_of_table_pos = 0;
				uint64_t greatest_pos = 0;

				uint8_t const* left_entry_disk_buf = nullptr;

				uint64_t entry_sort_key, entry_pos, entry_offset;
				uint64_t cached_entry_sort_key = 0;
				uint64_t cached_entry_pos = 0;
				uint64_t cached_entry_offset = 0;

				{ // scope for async rewriter
					EntryAsynRewriter rewriter( k, R_sort_manager.get(), num_threads, right_entry_size_bytes, left_new_pos, old_data );

					// Similar algorithm as Backprop, to read both L and R tables simultaneously
					for( uint64_t current_pos = 0;
							 !end_of_right_table || (current_pos - end_of_table_pos <= kReadMinusWrite); current_pos++ ) {

						old_data[current_pos % kReadMinusWrite].count = 0;

						if (end_of_right_table || current_pos <= greatest_pos) {
								while (!end_of_right_table) {
										if (should_read_entry) {
												if (right_reader_count == res2.table_sizes[table_index + 1]) {
														end_of_right_table = true;
														end_of_table_pos = current_pos;
														right_disk.FreeMemory();
														break;
												}
												// The right entries are in the format from backprop, (sort_key, pos, offset)
												uint8_t const* right_entry_buf = right_disk.Read(right_reader, p2_entry_size_bytes);
												right_reader += p2_entry_size_bytes;
												right_reader_count++;

												entry_sort_key =
														Util::SliceInt64FromBytes(right_entry_buf, right_sort_key_size);
												entry_pos = Util::SliceInt64FromBytes(
														right_entry_buf, right_sort_key_size, pos_size);
												entry_offset = Util::SliceInt64FromBytes(
														right_entry_buf, right_sort_key_size + pos_size, kOffsetSize);
										} else if (cached_entry_pos == current_pos) {
												entry_sort_key = cached_entry_sort_key;
												entry_pos = cached_entry_pos;
												entry_offset = cached_entry_offset;
										} else {
												break;
										}

										should_read_entry = true;

										if (entry_pos + entry_offset > greatest_pos) {
												greatest_pos = entry_pos + entry_offset;
										}
										if (entry_pos == current_pos) {
											old_data[entry_pos % kReadMinusWrite].add( entry_sort_key, entry_pos + entry_offset );
										} else {
												should_read_entry = false;
												cached_entry_sort_key = entry_sort_key;
												cached_entry_pos = entry_pos;
												cached_entry_offset = entry_offset;
												break;
										}
								}
						}

						if( current_pos < res2.table_sizes[table_index] ) {
							// The left entries are in the new format: (sort_key, new_pos), except for table
							// 1: (y, x).

							// We read the "new_pos" from the L table, which for table 1 is just x. For
							// other tables, the new_pos
							if (table_index == 1) {
									left_entry_disk_buf = res2.table1.ReadNext();

									// Only k bits, since this is x
									left_new_pos[current_pos % kCachedPositionsSize] =
											Util::SliceInt64FromBytes(left_entry_disk_buf, k);
							} else {
									left_entry_disk_buf = L_sort_manager->ReadEntry(left_reader);
									left_reader += new_pos_entry_size_bytes;

									// k+1 bits in case it overflows
									left_new_pos[current_pos % kCachedPositionsSize] =
											Util::SliceInt64FromBytes(left_entry_disk_buf, right_sort_key_size, k);
							}
						}

						// Rewrites each right entry as (line_point, sort_key)
						rewriter.next( current_pos );
					} // end of loop of first computation pass
				} // end srope for async rewriter

				// Remove no longer needed file
				res2.table1.FreeMemory();

				// Flush cache so all entries are written to buckets
				R_sort_manager->FlushCache();
				R_sort_manager->FreeMemory();

				computation_pass_1_timer.PrintElapsed("\tFirst computation pass time:");

				Timer computation_pass_2_timer;

				if (table_index > 1) {
						// Make sure all files are removed
						L_sort_manager.reset();
				}

				// In the second pass we read from R sort manager and write to L sort
				// manager, and they both handle table (table_index + 1)'s data. The
				// newly written table consists of (sort_key, new_pos). Add one extra
				// bit for 'new_pos' to the 7-th table as it may have more than 2^k
				// entries.
				new_pos_entry_size_bytes = cdiv(2 * k + (table_index == 6 ? 1 : 0), 8);

				// For tables below 6 we can only use a half of memory_size since it
				// will be sorted in the first pass of the next iteration together with
				// the next table, which will use the other half of memory_size.
				// Tables 6 and 7 will be sorted alone, so we use all memory for them.
				L_sort_manager = std::make_unique<SortManager>(
						//(table_index >= 5) ? memory_size : (memory_size / 2),
						memory_manager,
						num_buckets,
						log_num_buckets,
						new_pos_entry_size_bytes,
						tmp_dirname,
						filename + ".p3s.t" + std::to_string(table_index + 1),
						0, // bits_begin
						0, // strip_size
						k, // plot size
						4, // Phase
						table_index + 1,
						num_threads,
						(flags&NO_COMPACTION)==0 );

				// Now we will write on of the final tables, since we have a table sorted by line point.
				// The final table will simply store the deltas between each line_point, in fixed space
				// groups(parks), with a checkpoint in each group.
				uint8_t const sort_key_shift = 128 - right_sort_key_size;
				uint8_t const index_shift = sort_key_shift - (k + (table_index == 6 ? 1 : 0));
				const uint32_t park_size_bytes = EntrySizes::CalculateParkSize( k, table_index );

				// At this point we know how many park will be written,
				// than we can evalueate next table pointer
				uint64_t to_be_written = R_sort_manager->Count();
				final_table_begin_pointers[table_index + 1] = final_table_begin_pointers[table_index] +
						( to_be_written/kEntriesPerPark + ( (to_be_written%kEntriesPerPark) ? 1 : 0) ) * park_size_bytes;

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
								sort_buf_size = R_sort_manager->Read( sort_buf, sort_buf_size );
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

				// Start threads
				if( num_threads <= 1 )
					parking_thread();
				else {
					std::unique_ptr<std::thread> threads[num_threads];
					for( uint32_t t = 0; t < num_threads; t++ )
						threads[t].reset( new std::thread(parking_thread) );

					for( uint32_t t = 0; t < num_threads; t++ )
						threads[t]->join();
				}

				R_sort_manager.reset();
				L_sort_manager->FlushCache();

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
