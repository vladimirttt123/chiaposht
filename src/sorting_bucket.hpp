// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//    http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef SRC_CPP_SORTING_BUCKET_HPP_
#define SRC_CPP_SORTING_BUCKET_HPP_

#include <atomic>
#include <mutex>
#include "util.hpp"
#include "quicksort.hpp"
#include "disk_streams.hpp"

#define STATS_UINT_TYPE uint16_t

// ==================================================================================================
struct SortingBucket{
	const EntryCompactionData &compaction_data;
	const bool parallel_read;
	const uint8_t subbucket_bits;

	SortingBucket( const std::string &fileName, MemoryManager &memory_manager, STATS_UINT_TYPE* stats_mem,
								uint16_t bucket_no, uint8_t subbucket_bits, const EntryCompactionData &compaction_data, bool parallel_read = true )
			: compaction_data(compaction_data), parallel_read(parallel_read)
			, subbucket_bits(subbucket_bits)
			,	disk( new BucketStream( fileName, memory_manager, bucket_no, compaction_data ) )

#ifndef NDEBUG
			, bucket_no_( bucket_no )
#endif
			, statistics( stats_mem )
			, memory_manager(memory_manager)
	{
		// clear statistics
		memset( statistics, 0, sizeof(STATS_UINT_TYPE)<<subbucket_bits );
	}

	inline uint64_t Size() const { return compaction_data.entry_size_bytes*entries_count; }
	inline uint64_t Count() const { return entries_count; }
	inline uint16_t EntrySize() const { return compaction_data.entry_size_bytes; }

	std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
	uint64_t prepare_time = 0, read_time = 0, sort_time = 0;

	// Adds Entry to bucket - NOT THREADSAFE!!!
	inline void AddEntry( StreamBuffer &entry, const uint32_t &statistics_bits ){
		assert( disk ); // Check not closed
		assert( entry.size() >= compaction_data.entry_size_bytes );
		assert( entry.used() == compaction_data.entry_size_bytes  );
		assert( statistics_bits < ((uint32_t)1<<subbucket_bits) );
		assert( Util::ExtractNum( entry.get(), compaction_data.entry_size_bytes,
														(compaction_data.begin_bits + compaction_data.log_num_buckets), subbucket_bits ) == statistics_bits );

		statistics[statistics_bits]++;
		entries_count++;
		disk->Write( entry );
	}

	// Flush current buffer to disk
	// It should be at end of writting.
	void Flush( bool free_memory = false ){

		if( free_memory ){
			if( disk ) disk->FlusToDisk();

			if( entries_count > 0 ){
				// Save statistics to file
				statistics_file.reset( CreateFileStream( disk->getFileName() + ".statistics.tmp", memory_manager ) );
				StreamBuffer buf( (uint8_t*)statistics, sizeof(STATS_UINT_TYPE)<<subbucket_bits, sizeof(STATS_UINT_TYPE)<<subbucket_bits );
				statistics_file->Write( buf );
				((IBlockWriter*)statistics_file.get())->Close();
				buf.release(); // to not delete bufer that is external for the buffer
			}
			statistics = NULL; // no reset stats is external now
		} else {
			if( disk ) disk->Close(); // in case of cached enabled allow to clean for next buffer.
		}

	}

	void FreeMemory(){} // TODO: remove the proc

	/* Like destructor totaly removes the bucket including underlying file */
	void Remove(){
		if( statistics_file ) {
			statistics_file->Remove();
			statistics_file.reset();
		}
		if( !disk ) return; // already removed
		disk->EndToRead();
		disk.reset();
		statistics = NULL; // no reset stats is external now
	}


	void SortToMemory( uint8_t * memory, uint32_t num_threads = 2, uint32_t num_read_threads = 2, std::mutex * read_mutex = NULL ){
		assert( disk );

		if( entries_count == 0 ){
			if( read_mutex != NULL ) read_mutex->unlock();
			return; // nothing to sort
		}

		start_time = std::chrono::high_resolution_clock::now();
		const uint32_t subbuckets_count = 1<<subbucket_bits;
		STATS_UINT_TYPE *stats = statistics;
		auto local_stats( Util::allocate<STATS_UINT_TYPE>(statistics?0:subbuckets_count) ); // empty local stats

		if( statistics == NULL ) {
			assert( statistics_file );
			// Read statistics from file
			uint8_t * stats_pntr = (uint8_t*)( stats = local_stats.get() );
			for( StreamBuffer buf; statistics_file->Read( buf ) > 0; stats_pntr += buf.used() )
				memcpy( stats_pntr, buf.get(), buf.used() );

			assert( (STATS_UINT_TYPE*)stats_pntr == stats + subbuckets_count ); // read all stats

			statistics_file->Remove();
		}


		auto subbucket_positions_ptr = Util::allocate<uint64_t>( subbuckets_count );
		uint64_t* subbucket_positions = subbucket_positions_ptr.get();

		// Calculate initial subbuckets positions.
		subbucket_positions[0] = 0;
		for( uint32_t i = 1; i < subbuckets_count; i++ )
			subbucket_positions[i] = subbucket_positions[i-1] + (((uint64_t)stats[i-1])*compaction_data.entry_size_bytes);

		// Last position should end at count.
		if( subbucket_positions[subbuckets_count-1]/compaction_data.entry_size_bytes + stats[subbuckets_count-1] != Count() )
			throw InvalidValueException( "Sort statistics failed - may be some overflow happened. Try decrease subbucket bits number (-S parameter)." );


		auto end_prepare_time = std::chrono::high_resolution_clock::now();
		prepare_time = (end_prepare_time - start_time)/std::chrono::milliseconds(1);

		// Read from file to buckets, do not run threads if less than 1024 entries to read
		if( num_read_threads <= 1 || /*read_size*/ Count() < 1024 ){

			for( StreamBuffer buf; disk->Read( buf ) > 0; )
				FillSubBuckets( memory, buf.get(), buf.used(), subbucket_positions, compaction_data.entry_size_bytes );

			assert( subbucket_positions[0] == stats[0]*compaction_data.entry_size_bytes ); // check first bucket is full
			assert( subbucket_positions[subbuckets_count-1]/compaction_data.entry_size_bytes == Count() ); // check last bucket is full
		} else { // Read by threads
			if( num_read_threads <= 6 ){
				 // read in 2 directions one fills from start and second from back without synchronizations.
				auto subbucket_end_positions = std::make_unique<uint64_t[]>( subbuckets_count );
				// Calculate initial end subbuckets positions.
				for( uint32_t i = 1; i < subbuckets_count; i++ )
					subbucket_end_positions[i-1] = subbucket_positions[i] - compaction_data.entry_size_bytes;
				subbucket_end_positions[subbuckets_count-1] = compaction_data.entry_size_bytes*entries_count - compaction_data.entry_size_bytes;

				std::vector<std::thread> threads;
				std::mutex forward_mutex, backward_mutex;

				if( num_read_threads == 2 ){ // in case of 2 threads it is possible to do without locks on memory
					auto thread_function = [this, &memory]( uint64_t * positions, int16_t const direction){
						StreamBuffer buf( BUF_SIZE/compaction_data.entry_size_bytes*compaction_data.entry_size_bytes );
						std::unique_ptr<IBlockReader> reader( disk->CreateReader( true, parallel_read ) );

						while(  reader->Read( buf ) > 0 )
							FillSubBuckets( memory, buf.get(), buf.used(), positions, direction );
					};

					// start threads
					threads.emplace_back( thread_function, subbucket_positions, compaction_data.entry_size_bytes );
					threads.emplace_back( thread_function, subbucket_end_positions.get(), -compaction_data.entry_size_bytes );
				 } else {
					auto thread_function = [this, &memory]( uint64_t * positions, int16_t direction, std::mutex *mut ){
						StreamBuffer buf( BUF_SIZE/compaction_data.entry_size_bytes*compaction_data.entry_size_bytes );
						std::unique_ptr<IBlockReader> reader( disk->CreateReader() );

						while(  reader->Read( buf ) > 0 ){
							std::lock_guard<std::mutex> lk(*mut);
							FillSubBuckets( memory, buf.get(), buf.used(), positions, direction );
						}
					};

					for( uint t = 0; t < num_read_threads; t++ )
						if( (t&1) == 0 )
							threads.emplace_back( thread_function, subbucket_positions, compaction_data.entry_size_bytes, &forward_mutex );
						else
							threads.emplace_back( thread_function, subbucket_end_positions.get(), -compaction_data.entry_size_bytes, &backward_mutex );
				 }

				 for (auto& t : threads) t.join(); // wait for threads

#ifndef NDEBUG
				 // Chcek everything read as written.
				 for( uint32_t i = 0; i < subbuckets_count; i++ )
					 assert( subbucket_positions[i] == subbucket_end_positions[i] + compaction_data.entry_size_bytes );
#endif
			} else {
				auto a_subbucket_positions_ptr = Util::allocate<std::atomic_uint64_t>( subbuckets_count );
				auto a_subbucket_positions = a_subbucket_positions_ptr.get();
				for( uint32_t i = 0; i < subbuckets_count; i++ )
					a_subbucket_positions[i].store( subbucket_positions[i], std::memory_order_relaxed );

				// Define thread function
				auto thread_func = [this, &memory, &a_subbucket_positions](){
					StreamBuffer buf( BUF_SIZE/compaction_data.entry_size_bytes*compaction_data.entry_size_bytes );
					std::unique_ptr<IBlockReader> reader( disk->CreateReader() );

					while(  reader->Read( buf ) > 0 ){
						for( auto buf_ptr = buf.get() + buf.used() - compaction_data.entry_size_bytes; buf_ptr >= buf.get(); buf_ptr -= compaction_data.entry_size_bytes ){
							auto b_bits = (uint32_t)Util::ExtractNum64( buf_ptr, compaction_data.begin_bits_bucket, subbucket_bits );
							// Put next entry to its bucket and move pointer inside bucket to next entry position
							memcpy( memory + a_subbucket_positions[b_bits].fetch_add( compaction_data.entry_size_bytes, std::memory_order_relaxed), buf_ptr, compaction_data.entry_size_bytes );
						}
					}
				};

				{	// Start threads
					std::vector<std::thread> threads;
					for( uint32_t t = 0; t < num_read_threads; t++ )
						threads.emplace_back( thread_func );

					for (auto& t : threads) t.join();
				}

#ifndef NDEBUG
				// Chcek everything read as written.
				for( uint32_t i = 0; i < subbuckets_count-1; i++ )
					assert( a_subbucket_positions[i].load( std::memory_order_relaxed ) == subbucket_positions[i+1] );
				// last bucket at the count
				assert( a_subbucket_positions[subbuckets_count-1].load( std::memory_order_relaxed ) == Count()*compaction_data.entry_size_bytes );
#endif
			}

			// Fix bucket_positions for sorting step in such way that bucket_positions[i] reprsents end position of subbucket.
			subbucket_positions[0] = stats[0]*compaction_data.entry_size_bytes;
			for( uint32_t i = 1; i < subbuckets_count; i++ )
				subbucket_positions[i] = subbucket_positions[i-1] + stats[i]*compaction_data.entry_size_bytes;
		} // end read by threads

		// Clean underling resources
		if( read_mutex != NULL ) read_mutex->unlock(); // for external purposes mark end of reading.
		disk->EndToRead();
		total_reads = disk->getTotalReads(); // save reads statistics
		instant_reads = disk->getInstantReads();// save reads statistics
		disk.reset();
		read_time = (std::chrono::high_resolution_clock::now() - end_prepare_time)/std::chrono::milliseconds(1);

		if( num_threads <= 1 ){
			// Sort first bucket
			QuickSort::Sort( memory, compaction_data.entry_size_bytes, stats[0], compaction_data.begin_bits_bucket + subbucket_bits );

			for( uint32_t i = 1; i < subbuckets_count; i++ ){
				assert( subbucket_positions[i] == (subbucket_positions[i-1] + stats[i]*compaction_data.entry_size_bytes) ); // check any bucket is full
				QuickSort::Sort( memory + subbucket_positions[i-1], compaction_data.entry_size_bytes, stats[i], compaction_data.begin_bits_bucket + subbucket_bits );
			}
		}else {
			// Sort in threads
			auto threads = std::make_unique<std::thread[]>(num_threads);

			for( uint32_t i = 0; i < num_threads; i++ )
				threads[i] = std::thread( [&num_threads, &subbuckets_count, this, &subbucket_positions, &memory, &stats ]( uint32_t thread_no ){
						for( uint32_t t = thread_no; t < subbuckets_count; t += num_threads ){
							assert( t == 0 || subbucket_positions[t] == (subbucket_positions[t-1] + stats[t]*compaction_data.entry_size_bytes) ); // check any bucket is full

							uint64_t start_pos = t == 0 ? 0 : subbucket_positions[t-1];
							QuickSort::Sort( memory + start_pos, compaction_data.entry_size_bytes, stats[t], compaction_data.begin_bits_bucket + subbucket_bits );
						}
					}, i );

			for( uint32_t i = 0; i < num_threads; i++ )
				threads[i].join();
		}

		statistics = NULL; // clear ram - statistics is external than no reset

		sort_time = (std::chrono::high_resolution_clock::now() - start_time)/std::chrono::milliseconds(1);

		assert( Util::CheckSort( memory, compaction_data.entry_size_bytes, compaction_data.begin_bits, entries_count ) == 0 );
	}

	inline uint64_t getTotalReads() const { return total_reads; }
	inline uint64_t getInstantReads() const { return instant_reads; }

private:
	std::unique_ptr<std::mutex> addMutex = std::make_unique<std::mutex>();
	std::unique_ptr<BucketStream> disk;
#ifndef NDEBUG
	const uint16_t bucket_no_;
#endif
	// this is the number of entries for each subbucket
	STATS_UINT_TYPE* statistics;
	std::unique_ptr<IBlockWriterReader> statistics_file;
	uint64_t entries_count = 0;
	MemoryManager &memory_manager;

	uint64_t total_reads = 0, instant_reads = 0;


	/* Used for reading from disk */
	inline void FillSubBuckets( uint8_t* memory, const uint8_t* disk_buffer, const uint64_t &buf_size, uint64_t * subbucket_positions, int64_t direction ){
		for( uint32_t buf_ptr = 0; buf_ptr < buf_size; buf_ptr += compaction_data.entry_size_bytes ){
			// Define bucket to put into
			uint64_t subbucket = Util::ExtractNum64( disk_buffer + buf_ptr, compaction_data.begin_bits_bucket, subbucket_bits );

			// Check overflow
			assert( subbucket_positions[subbucket] < Size() );
			assert( Util::ExtractNum64( disk_buffer + buf_ptr, compaction_data.begin_bits, compaction_data.log_num_buckets ) == bucket_no_ );

			// Put next entry to its bucket
			memcpy( memory + subbucket_positions[subbucket], disk_buffer + buf_ptr, compaction_data.entry_size_bytes );
			// move pointer inside bucket to next entry position
			subbucket_positions[subbucket] += direction;
		}
	}

	#pragma region CacheBucket support {
	// This function in particular for CachBucket
	inline void BulkAdd( uint8_t * entries, const uint32_t * stats, uint32_t count ){
		assert( disk ); // Check not closed

		// Adding to statistics
		for( uint32_t i = 0; i < count; i++ )	{
			assert( Util::ExtractNum64( entries + i*compaction_data.entry_size_bytes
																, compaction_data.begin_bits + compaction_data.log_num_buckets, subbucket_bits ) == stats[i] );
			statistics[stats[i]]++;
		}

		StreamBuffer buf( 0 );
		// try to do without copy... may be need to change that write recieves const stream buffer
		buf.reset( entries, count * compaction_data.entry_size_bytes, count * compaction_data.entry_size_bytes );
		disk->Write( buf );
		buf.release();
		entries_count += count;
	}

	// This function in particular for CachBucket
	inline bool TryAddEntriesTS( uint8_t * entries, const uint32_t * stats, uint32_t count ){

		if( !addMutex->try_lock() ) return false;

		BulkAdd( entries, stats, count );

		addMutex->unlock();
		return true;
	}

	inline void AddEntriesTS( uint8_t * entries, const uint32_t * stats, uint32_t count ){
		std::lock_guard<std::mutex> lk( *addMutex.get() );
		BulkAdd( entries, stats, count );
	}
	friend struct CacheBucket;
	#pragma endregion CacheBucket support }

};

#endif // SRC_CPP_SORTING_BUCKET_HPP_
