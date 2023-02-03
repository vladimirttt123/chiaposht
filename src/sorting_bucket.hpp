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

#ifndef SRC_CPP_SORTING_BUCKET_HPP_
#define SRC_CPP_SORTING_BUCKET_HPP_

#include <mutex>
#include "disk.hpp"
#include "util.hpp"
#include "quicksort.hpp"



struct ABuffer{
	ABuffer( uint32_t entry_size )
		: size(0), allocated_length(BUF_SIZE - (BUF_SIZE%entry_size))
		, buffer(new uint8_t[allocated_length]){}

	inline bool Add( const uint8_t* &entry, const uint32_t &entry_size ){
		assert( (size+entry_size) <= allocated_length );
		memcpy(buffer.get() + size, entry, entry_size );
		return ( size += entry_size ) == allocated_length;
	}

	// returns number of NOT inserted records.
	inline uint32_t AddBulk( const uint8_t * &entries, const uint32_t &entry_size,
													 const uint32_t &first_entry, const uint32_t &num_entries ) {
		uint32_t insert_amount = std::min( allocated_length - size, entry_size*num_entries );
		memcpy( buffer.get() + size, entries + (first_entry*entry_size), insert_amount );
		size += insert_amount;
		return num_entries - insert_amount/entry_size;
	}

	inline bool IsEmpty() const { return size == 0;}
	inline bool IsFull() const { return size == allocated_length; }
	inline uint32_t Size() const { return size; }
	inline void EnsureAllocated() {
		if( !buffer )
			buffer.reset( new uint8_t[allocated_length] );
	}

	uint8_t * ReleaseBuffer( bool allocate_another = true ){
		auto res = buffer.release();
		size = 0;
		if( allocate_another )
			buffer.reset( new uint8_t[allocated_length] );
		return res;
	}


	// size of used buffer
	uint32_t size;
	// the buffer length
	const uint32_t allocated_length;
	std::unique_ptr<uint8_t[]> buffer;
};


// ==================================================================================================
struct SortingBucket{
	SortingBucket( const std::string &fileName, uint32_t bucket_no, uint32_t log_num_buckets, uint16_t entry_size, uint32_t bits_begin, uint8_t bucket_bits_count )
		:	disk( new FileDisk(fileName) )
		, bucket_no_( bucket_no )
		, log_num_buckets_( log_num_buckets )
		, bucket_bits_count_(bucket_bits_count)
		, entry_size_(entry_size)
		, bits_begin_(bits_begin)
		, statistics( new uint32_t[1<<bucket_bits_count] )
		, disk_buffer( new ABuffer( entry_size ) )
	{
		// clear statistics
		memset( statistics.get(), 0, sizeof(uint32_t)*(1<<bucket_bits_count) );
	}

	inline uint64_t Size() const { return entry_size_*entries_count; }
	inline uint64_t Count() const { return entries_count; }
	inline uint8_t	SubBucketBits() const { return bucket_bits_count_; }
	inline uint16_t EntrySize() const { return entry_size_; }

	uint64_t read_time = 0;
	uint64_t sort_time = 0;

	// Adds Entry to bucket
	inline void AddEntry( const uint8_t * entry, const uint32_t &statistics_bits ){
		assert(disk); // Check not closed
		assert( statistics_bits < ((uint32_t)1<<bucket_bits_count_) );
		assert( Util::ExtractNum( entry, entry_size_, bits_begin_, bucket_bits_count_ ) == statistics_bits );

		statistics[statistics_bits]++;
		entries_count++;
		if( disk_buffer->Add(entry, entry_size_ ) )
			Flush();
	}

	inline void AddEntryTS( const uint8_t *entry, const uint32_t &statistics_bits ){
		std::lock_guard<std::mutex> lk( *addMutex.get() );
		AddEntry( entry, statistics_bits );
	}

	inline void AddBulkTS( const uint8_t * entries, const uint32_t * stats, const uint32_t &count ){

		std::lock_guard<std::mutex> lk( *addMutex.get() );

		// Adding to buffer
		auto to_add = count;
		while( to_add > 0 ){
			 to_add = disk_buffer->AddBulk( entries, entry_size_, count-to_add, to_add );
			 if( disk_buffer->IsFull() ) Flush();
		}

		// Adding to statistics
		for( uint32_t i = 0; i < count; i++ )	{
			assert( Util::ExtractNum( entries + i*entry_size_, entry_size_, bits_begin_, bucket_bits_count_ ) == stats[i] );
			statistics[stats[i]]++;
		}
		entries_count += count;
	}

	// Flush current buffer to disk
	void Flush( bool free_memory = false ){
		if( !disk || disk_buffer->IsEmpty() ) return;

		// If more than 4 ( it should be number of threads ) buffers are waiting to write than wait to finish;
		if( disk_write_position + BUF_SIZE*4 < entry_size_*entries_count )
			WaitLastDiskWrite();

		int64_t size_to_write = disk_buffer->Size();
		disk_output_thread.reset( new std::thread( [this, size_to_write]( uint8_t* buffer, std::thread * prevThread ){
								auto buf_size = CompactBuffer( buffer, size_to_write );
								// Wait for previous disk write operation.
								if( prevThread ){
									prevThread->join();
									delete prevThread;
								}
								disk->Write( disk_write_position, buffer, buf_size );
								disk->Flush();
								delete[] buffer;
								disk_write_position += buf_size;
							}, disk_buffer->ReleaseBuffer( !free_memory ), disk_output_thread.release() )
					);
	}

	void CloseFile(){
		if( !disk ) return; // already removed
		WaitLastDiskWrite();
		disk->Close();
	}

	void FreeMemory(){
		memory_.reset();
	}

	/* Like destructor totaly removes the bucket including underlying file */
	void Remove(){
		if( !disk ) return; // already removed
		WaitLastDiskWrite();
		disk->Remove( disk_write_position == 0 );
		disk.reset();
		statistics.reset();
		disk_buffer.reset();
	}

	const inline uint8_t* get() const{
		assert(memory_);
		return memory_.get();
	}

	void SortToMemory( uint32_t num_threads = 2 ){
		if( memory_ ) return; // already sorted;

		// Init memory to sort into
		memory_.reset( new uint8_t[Size()] );
		uint8_t* memory = memory_.get();

		uint32_t buckets_count = 1<<bucket_bits_count_;
		auto bucket_positions = std::make_unique<uint64_t[]>( buckets_count );

		// Calculate initial buckets positions.
		bucket_positions[0] = 0;
		for( uint32_t i = 1; i < buckets_count; i++ )
			bucket_positions[i] = bucket_positions[i-1] + (statistics[i-1]*entry_size_);

		// Last position should end at count.
		assert( bucket_positions[buckets_count-1]/entry_size_ + statistics[buckets_count-1] == Count() );

		auto start_time = std::chrono::high_resolution_clock::now();
		uint64_t read_pos = 0, read_size = Size() - disk_buffer->Size();


		// Read from file to buckets, do not run threads if less than 1024 entries to read
		if( num_threads <= 1 || read_size < 1024*entry_size_ ){
			// if last buffer not flashed do not flash it, use it as already read.
			if( !disk_buffer->IsEmpty() )
				FillBuckets( memory, disk_buffer->buffer.get(), disk_buffer->Size(), bucket_positions.get(), entry_size_ );

			disk_buffer->EnsureAllocated();
			// Wait for last disk write finish
			WaitLastDiskWrite();

			while( read_pos < read_size ){
				auto buf_size = ReadCompactedBuffer( disk, disk_buffer->buffer.get(), disk_buffer->allocated_length, read_pos, read_size );
				GrowBuffer( disk_buffer->buffer.get(), buf_size );

				FillBuckets( memory, disk_buffer->buffer.get(), buf_size, bucket_positions.get(), entry_size_ );
			}
			// Clean Memory
			disk_buffer->buffer.reset();

			assert( bucket_positions[0] == statistics[0]*entry_size_ ); // check first bucket is full
			assert( bucket_positions[buckets_count-1]/entry_size_ == Count() ); // check last bucket is full
		} else {
			// read in 2 directions one fills from start and second from back.
			auto back_bucket_positions = std::make_unique<uint64_t[]>( buckets_count );
			for( uint32_t i = 0; i < buckets_count - 1; i++ )
				back_bucket_positions[i] = bucket_positions[i+1]-entry_size_;
			back_bucket_positions[buckets_count-1] = back_bucket_positions[buckets_count-2] + statistics[buckets_count-1] * entry_size_;

			// Wait for last disk write finish
			WaitLastDiskWrite();

			uint32_t max_threads = std::max( (uint32_t)1, num_threads/2 );
			const uint32_t num_sub_locks = ((uint32_t)1 ) << std::min( (uint32_t)bucket_bits_count_, (uint32_t)std::log2(max_threads)+2 );
			const uint32_t sub_locks_mask = num_sub_locks - 1;

			std::mutex mutRead;
			auto mutForward = std::make_unique<std::mutex[]>(num_sub_locks);
			auto mutBackward = std::make_unique<std::mutex[]>(num_sub_locks);
			// Define thread function
			auto thread_func = [&read_pos, &read_size, this, &memory , &num_sub_locks, &sub_locks_mask]( std::mutex *mutRead, std::mutex *mutWrite,  uint64_t* bucket_positions, const uint64_t direction ){
				while( true ){
					uint64_t buf_size;
					auto buf = std::make_unique<uint8_t[]>( disk_buffer->allocated_length );
					{
						std::lock_guard<std::mutex> lk(*mutRead);
						buf_size = ReadCompactedBuffer( disk, buf.get(), disk_buffer->allocated_length, read_pos, read_size );
						if( buf_size == 0 ) return;
					}
					GrowBuffer( buf.get(), buf_size );

					// Extract all bucket_bits
					uint32_t bucket_bits[buf_size/entry_size_];
					for( uint64_t buf_ptr = 0; buf_ptr < buf_size; buf_ptr += entry_size_ )
						bucket_bits[buf_ptr/entry_size_] = (uint32_t)Util::ExtractNum( buf.get() + buf_ptr, entry_size_, bits_begin_, bucket_bits_count_ );

					for( uint32_t l = 0; l < num_sub_locks; l++ ){
						std::lock_guard lk(mutWrite[l]);
						for( uint64_t buf_ptr = 0; buf_ptr < buf_size; buf_ptr += entry_size_ ){
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

			// Do it in parallel with forward reading
			if( !disk_buffer->IsEmpty() )
				FillBuckets( memory, disk_buffer->buffer.get(), disk_buffer->size, back_bucket_positions.get(), -(int64_t)entry_size_ );
			// Disk buffer can be cleaned we do not need it anymore.
			disk_buffer->buffer.reset();

			// Start threads
			{
				std::vector<std::thread> threads;
				for( uint32_t t = 0; t < max_threads; t++ ){
					threads.emplace_back( thread_func, &mutRead, mutForward.get(), bucket_positions.get(), entry_size_ );
					threads.emplace_back( thread_func, &mutRead, mutBackward.get(), back_bucket_positions.get(), -(int64_t)entry_size_ );
				}

				for (auto& t : threads)
						t.join();
			}

#ifndef NDEBUG
			// Chcek everything read as written.
			for( uint32_t i = 0; i < buckets_count; i++ )
				assert( (int64_t)bucket_positions[i] == (int64_t)back_bucket_positions[i] + entry_size_ );
#endif
			// Fix bucket_positions for sorting step
			bucket_positions[0] = statistics[0]*entry_size_;
			for( uint32_t i = 1; i < buckets_count; i++ )
				bucket_positions[i] = bucket_positions[i-1] + statistics[i]*entry_size_;
		}

		auto end_time = std::chrono::high_resolution_clock::now();
		read_time = (end_time - start_time)/std::chrono::milliseconds(1);

		if( num_threads <= 1 ){
			// Sort first bucket
			QuickSort::Sort( memory, entry_size_, statistics[0], bits_begin_ + bucket_bits_count_ );
			for( uint32_t i = 1; i < buckets_count; i++ ){
				assert( bucket_positions[i] == (bucket_positions[i-1] + statistics[i]*entry_size_) ); // check any bucket is full
				QuickSort::Sort( memory + bucket_positions[i-1], entry_size_, statistics[i], bits_begin_ + bucket_bits_count_ );
			}
		}else {
			// Sort in threads
			auto threads = std::make_unique<std::thread[]>(num_threads);
			for( uint32_t i = 0; i < num_threads; i++ )
				threads[i] = std::thread( [&num_threads, &buckets_count, this, &bucket_positions, &memory]( uint32_t thread_no ){
						for( uint32_t i = thread_no; i < buckets_count; i += num_threads ){
							assert( i == 0 || bucket_positions[i] == (bucket_positions[i-1] + statistics[i]*entry_size_) ); // check any bucket is full
							QuickSort::Sort( memory + (i == 0 ? 0 : bucket_positions[i-1]), entry_size_, statistics[i], bits_begin_ + bucket_bits_count_ );
						}
					}, i );

			for( uint32_t i = 0; i < num_threads; i++ )
				threads[i].join();
		}

		sort_time = (std::chrono::high_resolution_clock::now() - start_time)/std::chrono::milliseconds(1);
#ifndef NDEBUG
		// Check sort
		for( uint64_t i = 1; i < entries_count; i++ )
			assert( Util::MemCmpBits( memory + (i-1) * entry_size_, memory + i*entry_size_, entry_size_, bits_begin_ ) < 0 );
#endif
	}

	std::unique_ptr<std::mutex> addMutex = std::make_unique<std::mutex>();

private:
	std::unique_ptr<FileDisk> disk;
	const uint32_t bucket_no_;
	const uint32_t log_num_buckets_;
	const uint8_t bucket_bits_count_;
	std::unique_ptr<std::thread> disk_output_thread;
	uint64_t disk_write_position = 0;

	const uint16_t entry_size_;
	uint32_t bits_begin_;
	std::unique_ptr<uint32_t[]> statistics;
	std::unique_ptr<ABuffer> disk_buffer;
	uint64_t entries_count = 0;
	std::unique_ptr<uint8_t[]> memory_;


	inline void WaitLastDiskWrite(){
		if( disk_output_thread ){
			disk_output_thread->join();
			disk_output_thread.reset();
		}
	}

	/* Used for reading from disk */
	inline void FillBuckets( uint8_t* memory, const uint8_t* disk_buffer, const uint64_t &buf_size, uint64_t * bucket_positions, int64_t direction ){
		for( uint32_t buf_ptr = 0; buf_ptr < buf_size; buf_ptr += entry_size_ ){
			// Define bucket to put into
			uint64_t bucket_bits = Util::ExtractNum64( disk_buffer + buf_ptr, bits_begin_, bucket_bits_count_ );

			// Check overflow
			assert( bucket_positions[bucket_bits] < Size() );

			// Put next entry to its bucket
			memcpy( memory + bucket_positions[bucket_bits], disk_buffer + buf_ptr, entry_size_ );
			// move pointer inside bucket to next entry position
			bucket_positions[bucket_bits] += direction;
		}
	}


	inline uint32_t CompactBuffer( uint8_t * buffer, uint32_t buf_size ){
		if( log_num_buckets_ < 8 ) return buf_size;
		uint32_t res = buf_size;
		uint64_t new_entry_size = entry_size_;
		// Stage I: remove bucket number if possible by bytes alignings

		uint32_t start_bit = bits_begin_ - log_num_buckets_;
		uint64_t bytes_before = start_bit/8;
		new_entry_size--;
		if( (start_bit&7)==0 ){
			// if remove alligned by byte
			for( uint32_t i = bytes_before, j = i; i < buf_size; i++ )
				if( (i%entry_size_) != bytes_before )
					buffer[j++] = buffer[i];
				else {
					assert( buffer[i] == (bucket_no_&255) );
				}
		} else {
			uint8_t mask = (((uint8_t)0xff)>>(start_bit&7));
			for( uint32_t i = bytes_before, j = i; i < buf_size; i++ ){
				switch( (i%entry_size_) - bytes_before ){
					case 0:
						assert( (buffer[i]&mask) == (((bucket_no_>>(start_bit&7))&mask) ) );
						buffer[j] = buffer[i]&~mask;
						break;
					case 1:

						buffer[j++] |= buffer[i]&mask;
					break;
					default: buffer[j++] = buffer[i]; break;
				}
			}
		}
		res = buf_size - buf_size/entry_size_;

		return res;
	}

	inline uint32_t ReadCompactedBuffer( std::unique_ptr<FileDisk> &f, uint8_t *buffer, const uint32_t &buf_size, uint64_t &read_pos, const uint64_t &read_size ){
		auto to_read = std::min( (uint64_t)buf_size, read_size - read_pos );
		if( !to_read || log_num_buckets_ < 8 ) {
			if( to_read ) f->Read( read_pos, buffer, to_read );
			read_pos += to_read;
			return to_read;
		}

		int64_t compacted_size = to_read/entry_size_ *(entry_size_-1);
		f->Read( read_pos/entry_size_*(entry_size_-1), buffer, compacted_size );
		read_pos += to_read;

		return to_read;
	}

	inline void GrowBuffer( uint8_t* buffer, const uint32_t &buf_size ){
		if( log_num_buckets_ < 8 ) return; // buffer is not compacted
		const int64_t compacted_size = buf_size/entry_size_ *(entry_size_-1);
		const uint32_t start_bit = bits_begin_ - log_num_buckets_;
		uint64_t bytes_before = start_bit/8;
		uint8_t excluded = bucket_no_ & 255;
		if( (start_bit&7)==0 ){
			// if remove alligned by byte
			for( int64_t i = buf_size-1, j = compacted_size-1; i >= 0 ; i-- ){
				assert( j >= -1 );
				assert( i > 0 || j == -1 ); // lastx
				if( (i%entry_size_) == (int64_t)bytes_before  )
					buffer[i] = excluded;
				else buffer[i] = buffer[j--];
			}
		} else {
			uint8_t mask = (((uint8_t)0xff)>>(start_bit&7));
			uint8_t excluded_lo = (excluded>>(start_bit&7)), excluded_hi = excluded<<(8-(start_bit&7));
			for( int64_t i = buf_size-1, j = compacted_size-1; i >= (int64_t)bytes_before ; i-- ){
				switch( (i%entry_size_) - bytes_before ){
				case 0: buffer[i] = (buffer[j--]&(~mask)) | excluded_lo; break;
				case 1: buffer[i] = (buffer[j]&mask) | excluded_hi; break;
				default: 	buffer[i] =  buffer[j--]; break;
				}
			}
		}
#ifndef NDEBUG
		for( uint32_t i = 0; i < buf_size/entry_size_; i++ )
			assert( Util::ExtractNum( buffer + i*entry_size_, entry_size_, bits_begin_-log_num_buckets_, log_num_buckets_ ) == bucket_no_ );
#endif

	}
};

#endif // SRC_CPP_SORTING_BUCKET_HPP_
