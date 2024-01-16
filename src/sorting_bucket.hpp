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
const inline uint32_t CacheBucketSize = 64; // mesured in number of entries
const inline uint32_t CacheBucketSizeLimit = 60; // mesured in number of entries

struct EntriesCache{
	uint8_t * memory;
	const uint16_t entry_size, start_bit;
	const uint32_t entries_count;
	uint32_t used_count = 0;


	EntriesCache( uint8_t* memory, const uint16_t entry_size, const uint32_t entries_count, const uint16_t start_bit )
			: memory( memory ), entry_size(entry_size), start_bit(start_bit), entries_count(entries_count)
	{	}

	// returns is cache full
	inline bool merge( uint8_t *entries, uint16_t count ){
		assert( used_count + count <= entries_count );
		assert( BucketSort::CheckSort( entries, entry_size, start_bit, count ) == 0 );
		assert( BucketSort::CheckSort( get(), entry_size, start_bit, used_count ) == 0 );

		uint8_t *dst = memory + (entries_count-count-used_count)*entry_size;

		if( used_count == 0 ){
			memcpy( dst, entries, count*entry_size );
		} else {
			for( uint8_t *src = get(), *src_end = memory + entries_count*entry_size , *end = entries + count*entry_size;
					 entries < end; dst += entry_size ){
				assert( dst + entry_size <= src_end );

				if( Util::MemCmpBits( src, entries, entry_size, start_bit ) < 0 ){
					memcpy( dst, src, entry_size );
					src += entry_size;
					if( src == src_end ){
						assert( dst + entry_size + (end - entries) == src_end );
						memcpy( dst + entry_size, entries, end - entries );
						entries = end; // to halt the loop
					}
				} else {
					memcpy( dst, entries, entry_size );
					entries += entry_size;
				}
			}
		}

		used_count += count;

		assert( BucketSort::CheckSort( get(), entry_size, start_bit, used_count ) == 0 );

		return ( entries_count - used_count ) < CacheBucketSize;
	}

	inline uint8_t * shift( uint32_t count ){
		assert( BucketSort::CheckSort( get(), entry_size, start_bit, entries_count ) == 0 );
		assert( count <= used_count );

		used_count -= count;
		return memory + (entries_count-used_count - count)*entry_size;
	}

	inline uint8_t * get(){
		assert( used_count <= entries_count );
		return memory + (entries_count-used_count)*entry_size;
	}
};

// ==================================================================================================
struct SortingBucket{
	bool parallel_read;

	SortingBucket( const std::string &fileName, MemoryManager &memory_manager, STATS_UINT_TYPE* stats_mem,
								uint16_t bucket_no, uint8_t log_num_buckets, uint16_t entry_size,
								uint32_t begin_bits, uint8_t bucket_bits_count,
								bool enable_compression = true, int16_t sequence_start_bit = -1, bool parallel_read = true,
								uint8_t *presort_cache = nullptr, uint32_t presort_cache_entries_no = 0 )
			: parallel_read(parallel_read),
				disk( new BucketStream( fileName, memory_manager, bucket_no, log_num_buckets, entry_size, begin_bits, enable_compression, sequence_start_bit ) ),
#ifndef NDEBUG
				bucket_no_( bucket_no ),
				log_num_buckets_( log_num_buckets ),
#endif
				bucket_bits_count_(bucket_bits_count),
				entry_size_(entry_size),
				begin_bits_(begin_bits),
				statistics( stats_mem ),
				memory_manager(memory_manager),
				entries_cache( presort_cache, entry_size, presort_cache_entries_no, sequence_start_bit )
	{
		// clear statistics
		memset( statistics, 0, sizeof(STATS_UINT_TYPE)<<bucket_bits_count );
	}
	inline uint64_t SortedPosision() const { return sorted_pos; }
	inline uint32_t SortedPositionWait( uint64_t min_position ){
		uint32_t wait_times = 0;
		for( uint64_t current ; (current = sorted_pos.load(std::memory_order::relaxed)) < min_position; ){
			wait_times++;
			sorted_pos.wait( current, std::memory_order::relaxed );
		};
		return wait_times;
	}
	inline uint64_t Size() const { return entry_size_*entries_count; }
	inline uint64_t Count() const { return entries_count; }
	inline uint16_t EntrySize() const { return entry_size_; }

	std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
	uint64_t prepare_time = 0, read_time = 0, sort_time = 0;

	// Adds Entry to bucket - NOT THREADSAFE!!!
	inline void AddEntry( StreamBuffer &entry, const uint32_t &statistics_bits ){
		assert( disk ); // Check not closed
		assert( entry.size() >= entry_size_ );
		assert( entry.used() == entry_size_  );
		assert( statistics_bits < ((uint32_t)1<<bucket_bits_count_) );
		assert( Util::ExtractNum( entry.get(), entry_size_, begin_bits_, bucket_bits_count_ ) == statistics_bits );

		statistics[statistics_bits]++;
		entries_count++;
		disk->Write( entry );
	}



	// Flush current buffer to disk
	// It should be at end of writting.
	void Flush( bool free_memory = false ){

		if( entries_cache.used_count > 0 ){
			StreamBuffer buf( 0 );
			buf.reset( entries_cache.get(), entries_cache.used_count * entry_size_, entries_cache.used_count * entry_size_ );
			disk->Write( buf );
			buf.release();
			entries_cache.used_count = 0;
		}

		if( free_memory ){
			if( disk ) disk->FlusToDisk();

			if( entries_count > 0 ){
				// Save statistics to file
				statistics_file.reset( CreateFileStream( disk->getFileName() + ".statistics.tmp", memory_manager ) );
				StreamBuffer buf( (uint8_t*)statistics, sizeof(STATS_UINT_TYPE)<<bucket_bits_count_, sizeof(STATS_UINT_TYPE)<<bucket_bits_count_ );
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
		const uint32_t subbuckets_count = 1<<bucket_bits_count_;
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
			subbucket_positions[i] = subbucket_positions[i-1] + (((uint64_t)stats[i-1])*entry_size_);

		// Last position should end at count.
		if( subbucket_positions[subbuckets_count-1]/entry_size_ + stats[subbuckets_count-1] != Count() )
			throw InvalidValueException( "Sort statistics failed - may be some overflow happened. Try decrease subbucket bits number (-S parameter)." );


		auto end_prepare_time = std::chrono::high_resolution_clock::now();
		prepare_time = (end_prepare_time - start_time)/std::chrono::milliseconds(1);

		// Read from file to buckets, do not run threads if less than 1024 entries to read
		if( num_read_threads <= 1 || /*read_size*/ Size() < 1024*entry_size_ ){

			for( StreamBuffer buf; disk->Read( buf ) > 0; )
				 FillBuckets( memory, buf.get(), buf.used(), subbucket_positions, entry_size_ );

			assert( subbucket_positions[0] == stats[0]*entry_size_ ); // check first bucket is full
			assert( subbucket_positions[subbuckets_count-1]/entry_size_ == Count() ); // check last bucket is full
		} else { // Read by threads
			if( num_read_threads <= 6 ){
				 // read in 2 directions one fills from start and second from back without synchronizations.
				auto subbucket_end_positions = std::make_unique<uint64_t[]>( subbuckets_count );
				// Calculate initial end subbuckets positions.
				for( uint32_t i = 1; i < subbuckets_count; i++ )
					subbucket_end_positions[i-1] = subbucket_positions[i] - entry_size_;
				subbucket_end_positions[subbuckets_count-1] = entry_size_*entries_count - entry_size_;

				std::vector<std::thread> threads;
				std::mutex forward_mutex, backward_mutex;

				if( num_read_threads == 2 ){ // in case of 2 threads it is possible to do without locks on memory
					auto thread_function = [this, &memory]( uint64_t * positions, int16_t const direction){
						StreamBuffer buf( BUF_SIZE/entry_size_*entry_size_ );
						std::unique_ptr<IBlockReader> reader( disk->CreateReader( true, parallel_read ) );

						while(  reader->Read( buf ) > 0 )
							FillBuckets( memory, buf.get(), buf.used(), positions, direction );
					};

					// start threads
					threads.emplace_back( thread_function, subbucket_positions, entry_size_ );
					threads.emplace_back( thread_function, subbucket_end_positions.get(), -entry_size_ );
				 } else {
					auto thread_function = [this, &memory]( uint64_t * positions, int16_t direction, std::mutex *mut ){
						StreamBuffer buf( BUF_SIZE/entry_size_*entry_size_ );
						std::unique_ptr<IBlockReader> reader( disk->CreateReader() );

						while(  reader->Read( buf ) > 0 ){
							std::lock_guard<std::mutex> lk(*mut);
							FillBuckets( memory, buf.get(), buf.used(), positions, direction );
						}
					};

					for( uint t = 0; t < num_read_threads; t++ )
						if( (t&1) == 0 )
							threads.emplace_back( thread_function, subbucket_positions, entry_size_, &forward_mutex );
						else
							threads.emplace_back( thread_function, subbucket_end_positions.get(), -entry_size_, &backward_mutex );
				 }

				 for (auto& t : threads) t.join(); // wait for threads

#ifndef NDEBUG
				 // Chcek everything read as written.
				 for( uint32_t i = 0; i < subbuckets_count; i++ )
					 assert( subbucket_positions[i] == subbucket_end_positions[i] + entry_size_ );
#endif
			} else {
				auto a_subbucket_positions_ptr = Util::allocate<std::atomic_uint64_t>( subbuckets_count );
				auto a_subbucket_positions = a_subbucket_positions_ptr.get();
				for( uint32_t i = 0; i < subbuckets_count; i++ )
					a_subbucket_positions[i].store( subbucket_positions[i], std::memory_order_relaxed );

				// Define thread function
				auto thread_func = [this, &memory, &a_subbucket_positions](){
					StreamBuffer buf( BUF_SIZE/entry_size_*entry_size_ );
					std::unique_ptr<IBlockReader> reader( disk->CreateReader() );

					while(  reader->Read( buf ) > 0 ){
						for( auto buf_ptr = buf.get() + buf.used() - entry_size_; buf_ptr >= buf.get(); buf_ptr -= entry_size_ ){
							auto b_bits = (uint32_t)Util::ExtractNum64( buf_ptr, begin_bits_, bucket_bits_count_ );
							// Put next entry to its bucket and move pointer inside bucket to next entry position
							memcpy( memory + a_subbucket_positions[b_bits].fetch_add( entry_size_, std::memory_order_relaxed), buf_ptr, entry_size_ );
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
				assert( a_subbucket_positions[subbuckets_count-1].load( std::memory_order_relaxed ) == Count()*entry_size_ );
#endif
			}

			// Fix bucket_positions for sorting step in such way that bucket_positions[i] reprsents end position of subbucket.
			subbucket_positions[0] = stats[0]*entry_size_;
			for( uint32_t i = 1; i < subbuckets_count; i++ )
				subbucket_positions[i] = subbucket_positions[i-1] + stats[i]*entry_size_;
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
			QuickSort::Sort( memory, entry_size_, stats[0], begin_bits_ + bucket_bits_count_ );
			for( uint32_t i = 1; i < subbuckets_count; i++ ){
				assert( subbucket_positions[i] == (subbucket_positions[i-1] + stats[i]*entry_size_) ); // check any bucket is full
				QuickSort::Sort( memory + subbucket_positions[i-1], entry_size_, stats[i], begin_bits_ + bucket_bits_count_ );
			}
		}else {
			// Sort in threads
			auto threads = std::make_unique<std::thread[]>(num_threads);
			uint64_t sorted_positions[num_threads];
			memset( sorted_positions, 0, num_threads*8 ); // init for zeros...

			for( uint32_t i = 0; i < num_threads; i++ )
				threads[i] = std::thread( [&num_threads, &subbuckets_count, this, &subbucket_positions, &memory, &stats, &sorted_positions ]( uint32_t thread_no ){
						for( uint32_t t = thread_no; t < subbuckets_count; t += num_threads ){
							assert( t == 0 || subbucket_positions[t] == (subbucket_positions[t-1] + stats[t]*entry_size_) ); // check any bucket is full

							uint64_t start_pos = t == 0 ? 0 : subbucket_positions[t-1];

							QuickSort::Sort( memory + start_pos, entry_size_, stats[t], begin_bits_ + bucket_bits_count_ );

							// now update global sorted position
							sorted_positions[thread_no] = subbucket_positions[t];

							if( thread_no == 0 ){ // only 1st thread updates

								uint64_t minVal = sorted_positions[0];
								for( uint j = 1; j < num_threads; j++ )
									if( minVal > sorted_positions[j] )
										minVal = sorted_positions[j];

								bool isMin = true;
								for( uint j = 0; isMin && j < num_threads; j++ )
									isMin = minVal <= sorted_positions[j];

								if( isMin ){
									assert( minVal >= sorted_pos );

									sorted_pos.store( minVal, std::memory_order::relaxed );
									sorted_pos.notify_all();
								}
							}
						}
					}, i );

			for( uint32_t i = 0; i < num_threads; i++ )
				threads[i].join();
		}

		sorted_pos = entries_count*entry_size_; // mark finished sort
		sorted_pos.notify_all();

		statistics = NULL; // clear ram - statistics is external no than no reset

		sort_time = (std::chrono::high_resolution_clock::now() - start_time)/std::chrono::milliseconds(1);
#ifndef NDEBUG
		// Check sort
		for( uint64_t i = 1; i < entries_count; i++ ){
			auto cmp = Util::MemCmpBits( memory + (i-1) * entry_size_, memory + i*entry_size_, entry_size_, begin_bits_ );
			assert( cmp < 0 );
		}
#endif
	}

	inline uint64_t getTotalReads() const { return total_reads; }
	inline uint64_t getInstantReads() const { return instant_reads; }

private:
	std::unique_ptr<std::mutex> addMutex = std::make_unique<std::mutex>();
	std::unique_ptr<BucketStream> disk;
#ifndef NDEBUG
	const uint16_t bucket_no_;
	const uint8_t log_num_buckets_;
#endif
	const uint8_t bucket_bits_count_;
	const uint16_t entry_size_;
	const uint16_t begin_bits_;
	// this is the number of entries for each subbucket
	STATS_UINT_TYPE* statistics;
	std::unique_ptr<IBlockWriterReader> statistics_file;
	uint64_t entries_count = 0;
	std::atomic_uint64_t sorted_pos = 0;
	MemoryManager &memory_manager;

	EntriesCache entries_cache;

	uint64_t total_reads = 0, instant_reads = 0;



	/* Used for reading from disk */
	inline void FillBuckets( uint8_t* memory, const uint8_t* disk_buffer, const uint64_t &buf_size, uint64_t * bucket_positions, int64_t direction ){
		for( uint32_t buf_ptr = 0; buf_ptr < buf_size; buf_ptr += entry_size_ ){
			// Define bucket to put into
			uint64_t bucket_bits = Util::ExtractNum64( disk_buffer + buf_ptr, begin_bits_, bucket_bits_count_ );

			// Check overflow
			assert( bucket_positions[bucket_bits] < Size() );
			assert( Util::ExtractNum64( disk_buffer + buf_ptr, begin_bits_ - log_num_buckets_, log_num_buckets_ ) == bucket_no_ );

			// Put next entry to its bucket
			memcpy( memory + bucket_positions[bucket_bits], disk_buffer + buf_ptr, entry_size_ );
			// move pointer inside bucket to next entry position
			bucket_positions[bucket_bits] += direction;
		}
	}

	#pragma region CacheBucket support {
	// This function in particular for CachBucket
	inline void BulkAdd( uint8_t * entries, const uint32_t * stats, uint32_t count ){
		assert( disk ); // Check not closed

		// Adding to statistics
		for( uint32_t i = 0; i < count; i++ )	{
			assert( Util::ExtractNum( entries + i*entry_size_, entry_size_, begin_bits_, bucket_bits_count_ ) == stats[i] );
			statistics[stats[i]]++;
		}

		StreamBuffer buf( 0 );
		if( entries_cache.memory != nullptr ){
			if( entries_cache.merge( entries, count ) ){
				buf.reset( entries_cache.shift( count ), count * entry_size_, count * entry_size_ );
				disk->Write( buf );
			}
		} else {
			// try to do without copy... may be need to change that write recieves const stream buffer
			buf.reset( entries, count * entry_size_, count * entry_size_ );
			disk->Write( buf );
		}
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
