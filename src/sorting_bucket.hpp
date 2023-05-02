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

#include <mutex>
#include "util.hpp"
#include "quicksort.hpp"
#include "disk_streams.hpp"


// ==================================================================================================
struct SortingBucket{
	SortingBucket( const std::string &fileName, MemoryManager &memory_manager, uint16_t bucket_no, uint8_t log_num_buckets, uint16_t entry_size,
								 uint32_t begin_bits, uint8_t bucket_bits_count, bool enable_compression = true, uint16_t sequence_start_bit = -1 )
		:	disk( new BucketStream( fileName, memory_manager, bucket_no, log_num_buckets, entry_size, begin_bits, enable_compression, sequence_start_bit ) )
		, bucket_no_( bucket_no )
		, log_num_buckets_( log_num_buckets )
		, bucket_bits_count_(bucket_bits_count)
		, entry_size_(entry_size)
		, begin_bits_(begin_bits)
		, statistics( new uint32_t[1<<bucket_bits_count] )
		, memory_manager(memory_manager)
	{
		// clear statistics
		memset( statistics.get(), 0, sizeof(uint32_t)*(1<<bucket_bits_count) );
	}

	inline uint64_t Size() const { return entry_size_*entries_count; }
	inline uint64_t Count() const { return entries_count; }
	inline uint16_t EntrySize() const { return entry_size_; }

	uint64_t prepare_time = 0, read_time = 0, sort_time = 0;

	// Adds Entry to bucket
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

	inline void AddBulkTS( StreamBuffer & entries, const uint32_t * stats ){

		std::lock_guard<std::mutex> lk( *addMutex.get() );

		auto count = entries.used() / entry_size_;

		// Adding to statistics
		for( uint32_t i = 0; i < count; i++ )	{
			assert( Util::ExtractNum( entries.get() + i*entry_size_, entry_size_, begin_bits_, bucket_bits_count_ ) == stats[i] );
			statistics[stats[i]]++;
		}

		disk->Write( entries );
		entries_count += count;
	}

	// Flush current buffer to disk
	// It should be at end of writting.
	void Flush( bool free_memory = false ){

		if( free_memory ){
			disk->FlusToDisk();

			if( entries_count > 0 ){
				// Save statistics to file
				statistics_file.reset( CreateFileStream( disk->getFileName() + ".statistics.tmp", memory_manager ) );
				StreamBuffer buf(  (uint8_t*)statistics.release(), sizeof(uint32_t)<<bucket_bits_count_, sizeof(uint32_t)<<bucket_bits_count_ );
				statistics_file->Write( buf );
				((IBlockWriter*)statistics_file.get())->Close();
			}
			statistics.reset();
		}
	}

	void FreeMemory(){ memory_.reset(); }

	/* Like destructor totaly removes the bucket including underlying file */
	void Remove(){
		if( statistics_file ) {
			statistics_file->Remove();
			statistics_file.reset();
		}
		if( !disk ) return; // already removed
		disk->EndToRead();
		disk.reset();
		statistics.reset();
	}

	const inline uint8_t* get() const{
		assert(memory_);
		return memory_.get();
	}

	void SortToMemory( uint32_t num_threads = 2 ){
		assert( disk );

		if( memory_ || entries_count == 0 ) return; // already sorted or nothing to sort

		auto start_time = std::chrono::high_resolution_clock::now();
		if( !statistics ) {
			assert( statistics_file );
			// Read statistics from file
			StreamBuffer buf;
			statistics.reset( new uint32_t[1<<bucket_bits_count_] );// allocate ram for statistics
			for( uint64_t read = 0; read < sizeof(uint32_t)<<bucket_bits_count_; read += buf.used() ){
				statistics_file->Read( buf );
				memcpy( statistics.get() + read, buf.get(), buf.used() );
			}
			statistics_file->Remove();
			statistics_file.reset();
		}

		// Init memory to sort into
		memory_.reset( Util::NewSafeBuffer(Size()) );
		uint8_t* memory = memory_.get();

		uint32_t buckets_count = 1<<bucket_bits_count_;
		auto bucket_positions = std::make_unique<uint64_t[]>( buckets_count );

		// Calculate initial buckets positions.
		bucket_positions[0] = 0;
		for( uint32_t i = 1; i < buckets_count; i++ )
			bucket_positions[i] = bucket_positions[i-1] + (statistics[i-1]*entry_size_);

		// Last position should end at count.
		assert( bucket_positions[buckets_count-1]/entry_size_ + statistics[buckets_count-1] == Count() );


		auto end_prepare_time = std::chrono::high_resolution_clock::now();
		prepare_time = (end_prepare_time - start_time)/std::chrono::milliseconds(1);

		// Read from file to buckets, do not run threads if less than 1024 entries to read
		if( num_threads <= 1 || /*read_size*/ Size() < 1024*entry_size_ ){

			for( StreamBuffer buf; disk->Read( buf ) > 0; )
				 FillBuckets( memory, buf.get(), buf.used(), bucket_positions.get(), entry_size_ );

			assert( bucket_positions[0] == statistics[0]*entry_size_ ); // check first bucket is full
			assert( bucket_positions[buckets_count-1]/entry_size_ == Count() ); // check last bucket is full
		} else { // Read by threads
			// read in 2 directions one fills from start and second from back.
			auto back_bucket_positions = std::make_unique<uint64_t[]>( buckets_count );
			for( uint32_t i = 0; i < buckets_count - 1; i++ )
				back_bucket_positions[i] = bucket_positions[i+1]-entry_size_;
			back_bucket_positions[buckets_count-1] = back_bucket_positions[buckets_count-2] + statistics[buckets_count-1] * entry_size_;

			uint32_t max_threads = std::max( (uint32_t)1, num_threads/2 );
			const uint32_t num_sub_locks_bits = std::min( 8U, std::min( (uint32_t)bucket_bits_count_, (uint32_t)std::log2(max_threads)+1 ) );
			const uint32_t num_sub_locks = 1U << num_sub_locks_bits;
			const uint32_t sub_locks_mask = num_sub_locks - 1;

			auto mutForward = std::make_unique<std::mutex[]>(num_sub_locks);
			auto mutBackward = std::make_unique<std::mutex[]>(num_sub_locks);
			// Define thread function
			auto thread_func = [this, &memory , &num_sub_locks, &sub_locks_mask]( std::mutex *mutWrite,  uint64_t* bucket_positions, const int64_t direction ){
				StreamBuffer buf( BUF_SIZE/entry_size_*entry_size_ );
				std::unique_ptr<IBlockReader> reader( disk->CreateReader() );

				while(  reader->Read( buf ) > 0 ){
					uint32_t bucket_bits[buf.size()/entry_size_ + 1];

					// Extract all bucket_bits
					for( uint64_t buf_ptr = 0; buf_ptr < buf.used(); buf_ptr += entry_size_ )
						bucket_bits[buf_ptr/entry_size_] = (uint32_t)Util::ExtractNum64( buf.get() + buf_ptr, begin_bits_, bucket_bits_count_ );

					for( uint32_t l = 0; l < num_sub_locks; l++ ){
						std::lock_guard lk( mutWrite[l] );
						for( uint64_t buf_ptr = 0; buf_ptr < buf.used(); buf_ptr += entry_size_ ){
							auto b_bits = bucket_bits[buf_ptr/entry_size_];
							if( (b_bits&sub_locks_mask) == l ){
								// Put next entry to its bucket
								memcpy( memory + bucket_positions[b_bits], buf.get() + buf_ptr, entry_size_ );
								// move pointer inside bucket to next entry position
								bucket_positions[b_bits] += direction;
							}
						}
					}
				}

			};

			// Start threads
			{
				std::vector<std::thread> threads;
				for( uint32_t t = 0; t < max_threads; t++ ){
					threads.emplace_back( thread_func, mutForward.get(), bucket_positions.get(), entry_size_ ); // thread forward
					threads.emplace_back( thread_func, mutBackward.get(), back_bucket_positions.get(), -(int64_t)entry_size_ ); // thread backward
				}

				for (auto& t : threads)
						t.join();
			}

#ifndef NDEBUG
			// Chcek everything read as written.
			for( uint32_t i = 0; i < buckets_count; i++ )
				assert( (int64_t)bucket_positions[i] == (int64_t)back_bucket_positions[i] + entry_size_ );
#endif
			// Clean underling resources
			disk->EndToRead();


			// Fix bucket_positions for sorting step in such way that bucket_positions[i] reprsents end position of subbucket.
			bucket_positions[0] = statistics[0]*entry_size_;
			for( uint32_t i = 1; i < buckets_count; i++ )
				bucket_positions[i] = bucket_positions[i-1] + statistics[i]*entry_size_;
		}

		read_time = (std::chrono::high_resolution_clock::now() - end_prepare_time)/std::chrono::milliseconds(1);

		if( num_threads <= 1 ){
			// Sort first bucket
			QuickSort::Sort( memory, entry_size_, statistics[0], begin_bits_ + bucket_bits_count_ );
			for( uint32_t i = 1; i < buckets_count; i++ ){
				assert( bucket_positions[i] == (bucket_positions[i-1] + statistics[i]*entry_size_) ); // check any bucket is full
				QuickSort::Sort( memory + bucket_positions[i-1], entry_size_, statistics[i], begin_bits_ + bucket_bits_count_ );
			}
		}else {
			// Sort in threads
			auto threads = std::make_unique<std::thread[]>(num_threads);
			for( uint32_t i = 0; i < num_threads; i++ )
				threads[i] = std::thread( [&num_threads, &buckets_count, this, &bucket_positions, &memory]( uint32_t thread_no ){
						for( uint32_t i = thread_no; i < buckets_count; i += num_threads ){
							assert( i == 0 || bucket_positions[i] == (bucket_positions[i-1] + statistics[i]*entry_size_) ); // check any bucket is full
							QuickSort::Sort( memory + (i == 0 ? 0 : bucket_positions[i-1]), entry_size_, statistics[i], begin_bits_ + bucket_bits_count_ );
						}
					}, i );

			for( uint32_t i = 0; i < num_threads; i++ )
				threads[i].join();
		}

		sort_time = (std::chrono::high_resolution_clock::now() - start_time)/std::chrono::milliseconds(1);
#ifndef NDEBUG
		// Check sort
		for( uint64_t i = 1; i < entries_count; i++ )
			assert( Util::MemCmpBits( memory + (i-1) * entry_size_, memory + i*entry_size_, entry_size_, begin_bits_ ) < 0 );
#endif
	}

private:
	std::unique_ptr<std::mutex> addMutex = std::make_unique<std::mutex>();
	std::unique_ptr<BucketStream> disk;
	const uint16_t bucket_no_;
	const uint8_t log_num_buckets_;
	const uint8_t bucket_bits_count_;

	const uint16_t entry_size_;
	const uint16_t begin_bits_;
	// this is the number of entries for each subbucket
	std::unique_ptr<uint32_t[]> statistics;
	std::unique_ptr<IBlockWriterReader> statistics_file;
	uint64_t entries_count = 0;
	std::unique_ptr<uint8_t[]> memory_;
	MemoryManager &memory_manager;


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

};

#endif // SRC_CPP_SORTING_BUCKET_HPP_
