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

#ifndef SRC_CPP_FAST_SORT_ON_DISK_HPP_
#define SRC_CPP_FAST_SORT_ON_DISK_HPP_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <chrono>
using namespace std::chrono_literals; // for operator""ns;

#include "./bits.hpp"
#include "./calculate_bucket.hpp"
#include "disk.hpp"
#include "exceptions.hpp"
#include "sorting_bucket.hpp"

const uint32_t CacheBucketSize = 256; // mesured in number of entries
// Small bucket used in thread writings
struct CacheBucket{
	explicit CacheBucket( SortingBucket &cacheFor )
		: entries( CacheBucketSize*cacheFor.EntrySize() )
		, entry_size( cacheFor.EntrySize() )
		, parent( cacheFor )
	{	}

	inline void Add( const uint8_t *entry, const uint32_t &stats ){
		if( entries.isFull() ) Flush();
		statistics[entries.used()/entry_size] = stats;
		entries.add( entry, entry_size );
	}

	inline void Flush(){
		if( entries.used() > 0 ){
			parent.AddBulkTS( entries, statistics );
			entries.setUsed( 0 );
		}
	}

	~CacheBucket(){Flush();}
private:
	uint32_t statistics[CacheBucketSize];
	StreamBuffer entries;
	const uint16_t entry_size;
	SortingBucket &parent;
};

struct SortedBucketBuffer{
	const uint64_t buffer_size;
	const bool in_background;

	SortedBucketBuffer( uint64_t buffer_size, bool in_background )
			: buffer_size( buffer_size )
			, in_background( in_background )
			, bucket_buffer( Util::NewSafeBuffer( buffer_size ) )
	{
	}

	~SortedBucketBuffer(){ FinishSort(); }

	void GetFrom( SortedBucketBuffer & other ){
		assert( other.buffer_size == this->buffer_size );

		FinishSort();


		bucket_buffer.swap( other.bucket_buffer );
		start_position = other.start_position;
		end_position = other.end_position;

		bucket = other.bucket;
		other.bucket = NULL;
		sorting_thread = other.sorting_thread;
		other.sorting_thread = NULL;

		bucket_no.store(other.bucket_no.load(std::memory_order_relaxed), std::memory_order_relaxed );
		other.bucket_no.store( -1, std::memory_order_relaxed );
	}

	// returns true if current sort just finished.
	void WaitForSortedTo( uint64_t position ){
		while( sorting_thread != NULL && position >= start_position && position <= end_position
					 && position > ( bucket->SortedCount() * bucket->EntrySize() + start_position ) )
				std::this_thread::sleep_for( 100ns );
	}

	bool isReading() const {
		return sorting_thread != NULL && bucket->SortedCount() == 0;
	}

	inline const uint8_t* buffer() const { return bucket_buffer.get(); }
	inline SortingBucket * Bucket() const { return bucket; }
	inline int32_t BucketNo() const { return bucket_no.load( std::memory_order_relaxed ); }
	inline uint64_t StartPosition() const { return start_position; }
	inline uint64_t EndPosition() const { return end_position; }

	// used when no enough memory for background presort
	void StartSortingNext( SortingBucket * next_bucket, int num_background_threads, int num_read_threads ){
		FinishSort();
		bucket_no++;
		start_position += bucket->Size();
		bucket = next_bucket;
		end_position = start_position + next_bucket->Size();
		RunSort( num_background_threads, num_read_threads );
	}

	void TryStartSorting( uint bucket_no, uint64_t start_position, SortingBucket * bucket, int num_background_threads, int num_read_threads ){
		int_fast32_t expected = -1;
		if( this->bucket_no.compare_exchange_weak( expected, bucket_no, std::memory_order_relaxed, std::memory_order_relaxed ) ){
			this->bucket = bucket;
			this->start_position = start_position;
			this->end_position = start_position + bucket->Size();

			RunSort( num_background_threads, num_read_threads );
		}
	}

	void StartSorting( uint bucket_no, uint64_t start_position, SortingBucket * bucket, int num_background_threads, int num_read_threads ){
		assert( bucket_no >= 0 );
		FinishSort();

		this->bucket_no.store( bucket_no, std::memory_order_relaxed );
		this->bucket = bucket;
		this->start_position = start_position;
		this->end_position = start_position + bucket->Size();

		RunSort( num_background_threads, num_read_threads );
	}
private:
	std::unique_ptr<uint8_t[]> bucket_buffer;
	std::atomic_int_fast32_t bucket_no = -1;
	uint64_t start_position, end_position;
	SortingBucket *bucket = NULL;
	std::thread *sorting_thread = NULL;

	inline void FinishSort(){
		if( sorting_thread != NULL ){
			sorting_thread->join();
			delete sorting_thread;
			sorting_thread = NULL;
		}
	}
	void RunSort( int num_background_threads, int num_read_threads ){
		assert( this->bucket_no >= 0 );
		if( in_background )
			sorting_thread = new std::thread( [num_background_threads, num_read_threads]( SortingBucket* bucket, uint8_t* buf ){
				bucket->SortToMemory( buf, num_background_threads, num_read_threads );}, this->bucket, this->bucket_buffer.get() );
		else
			bucket->SortToMemory( bucket_buffer.get(), num_background_threads, num_read_threads );
	}
};

class SortManager : public Disk, IReadDiskStream {
public:
	SortManager(
			MemoryManager &memory_manager,
			uint32_t const num_buckets,
			uint32_t const log_num_buckets,
			uint16_t const entry_size,
			const std::string &tmp_dirname,
			const std::string &filename,
			uint32_t const begin_bits,
			uint64_t const stripe_size,
			uint8_t k, uint8_t phase, uint8_t table_index,
			uint32_t num_threads = 2,
			bool enable_compaction = true )

			: memory_manager(memory_manager)
			, entry_size_(entry_size)
			, begin_bits_(begin_bits)
			, log_num_buckets_(log_num_buckets)
			, prev_bucket_buf_size(
					2 * (stripe_size + 10 * (kBC / pow(2, kExtraBits))) * entry_size)
			, num_threads( num_threads )
			, num_background_threads( std::max( num_threads>1?2U:1U, num_threads/2 ) )
			, num_read_threads( std::max( num_threads>1?2U:1U, num_threads/3 ) )
			, subbucket_bits( std::min( (uint8_t)(32-log_num_buckets), std::max( (uint8_t)2, (uint8_t)(k - log_num_buckets - kSubBucketBits) ) ) )
			, k_(k), phase_(phase), table_index_(table_index)
			, stats_mask( ( (uint64_t)1<<subbucket_bits)-1 )
	{
		assert( subbucket_bits > 0 );
		// Cross platform way to concatenate paths, gulrak library.
		std::vector<fs::path> bucket_filenames = std::vector<fs::path>();

		buckets_.reserve(num_buckets);
		for (size_t bucket_i = 0; bucket_i < num_buckets; bucket_i++) {
				std::ostringstream bucket_number_padded;
				bucket_number_padded << std::internal << std::setfill('0') << std::setw(3) << bucket_i;

				fs::path const bucket_filename =
						fs::path(tmp_dirname) /
						fs::path(filename + ".sort_bucket_" + bucket_number_padded.str() + ".tmp");
				uint16_t sequence_start = -1;
				if( k >= 32 ){ // for k < 32 sqeunce compaction not working. it need long enough entry size
					switch (phase) {
						case 1: sequence_start = table_index == 1 ? k : (k+kExtraBits); break;
						case 2: sequence_start = 0; break;
						case 4: sequence_start = k; break;
					}
				}
				buckets_.emplace_back( SortingBucket( bucket_filename.string(), memory_manager, bucket_i, log_num_buckets_,
																							entry_size, begin_bits_ + log_num_buckets, subbucket_bits,
																							enable_compaction, sequence_start ) );
		}
	}

		// Class to support writing to sort cache by threads in safe way
		struct ThreadWriter{
			explicit ThreadWriter( SortManager &parent ) : parent_(parent)
						, buckets_cache(new std::unique_ptr<CacheBucket>[parent.buckets_.size()])
//						, begin_bytes( parent.begin_bits_/8 ), begin_bits( parent.begin_bits_&7 )
//						, bits_shift( 32 - parent.log_num_buckets_ - parent.subbucket_bits )
			{
				for( uint32_t i = 0; i < parent.buckets_.size(); i++ )
					buckets_cache[i].reset( new CacheBucket( parent.buckets_[i] ) );
			}

			inline void Add( const uint8_t *entry ){
//				uint64_t const bucket_index = Util::ExtractNum32( entry + begin_bytes, begin_bits ) >> bits_shift;
				uint64_t const bucket_index =
						Util::ExtractNum64(entry, parent_.begin_bits_, parent_.log_num_buckets_ + parent_.subbucket_bits );

				buckets_cache[bucket_index>>parent_.subbucket_bits]->Add( entry, bucket_index & parent_.stats_mask );
			}

			inline void Add( uint128_t entry ){
				uint8_t bytes[16];
				Util::IntTo16Bytes(bytes, entry);
				Add( bytes );
			}

			inline void Flush(){
				for( uint32_t i = 0; i < parent_.buckets_.size(); i++ )
					buckets_cache[i]->Flush();
			}

			~ThreadWriter(){Flush();}
		private:
			SortManager &parent_;
			std::unique_ptr<std::unique_ptr<CacheBucket>[]> buckets_cache;
//			const uint8_t begin_bytes;
//			const uint8_t begin_bits;
//			const uint8_t bits_shift;
		}; // ====== END OF ThreadWriter

		// returned number of entries
		inline uint64_t Count() const {
			uint64_t res = 0;
			for( auto &b : buckets_ )
				res += b.Count();
			return res;
		}


		inline void AddToCache( StreamBuffer &entry )
    {
        uint64_t const bucket_index =
						Util::ExtractNum64( entry.get(), begin_bits_, log_num_buckets_ + subbucket_bits );
				buckets_[bucket_index>>subbucket_bits].AddEntry( entry, bucket_index & stats_mask );
		}


		uint8_t const* Read(uint64_t begin, uint64_t length) override
    {
        assert(length <= entry_size_);
        return ReadEntry(begin);
    }

    void Write(uint64_t, uint8_t const*, uint64_t) override
    {
        assert(false);
        throw InvalidStateException("Invalid Write() called on SortManager");
    }

    void Truncate(uint64_t new_size) override
    {
        if (new_size != 0) {
            assert(false);
            throw InvalidStateException("Invalid Truncate() called on SortManager");
        }

        FlushCache();
        FreeMemory();
    }

		std::string GetFileName() override { return "<SortManager>"; }

    void FreeMemory() override
    {
			for (auto& b : buckets_)
				b.FreeMemory();

			prev_bucket_buf_.reset();

			if( sorted_current ){
				memory_manager.release( sorted_current->buffer_size * (sorted_next?2:1) ); // release current and next bucket's memory
				sorted_current.reset();
				sorted_next.reset();
			}
    }

		//#Start IReadStream implementation
		uint32_t Read( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override {
			if( stream_read_position == 0 ) StartFirstBucket(); // initiate if needed

			for( uint32_t buf_start = 0; buf_start < buf_size; ){
				if( hasMoreBuckets && stream_read_position == sorted_current->EndPosition() )
					SwitchNextBucket();
				auto to_copy = std::min( sorted_current->EndPosition() - stream_read_position, (uint64_t)(buf_size-buf_start) );
				if( to_copy == 0 ) return buf_start; // is it possible case?
				EnsureSortedTo( stream_read_position + to_copy );
				memcpy( buf.get() + buf_start, memory_start() + stream_read_position - sorted_current->StartPosition(), to_copy );
				buf_start += to_copy;
				stream_read_position += to_copy;
			}

			return buf_size;
		};
		bool atEnd() const override{
			return !hasMoreBuckets && sorted_current && (uint32_t)sorted_current->BucketNo() >= buckets_.size()
						 && stream_read_position >= sorted_current->EndPosition();
		};
		void Close() override { FreeMemory(); };

		inline uint64_t GetReadPosition() const { return stream_read_position; }
		//#End IReadStream implementation

		const uint8_t *ReadEntry( uint64_t position )
    {
			if( !sorted_current ) StartFirstBucket(); // ensure sorting started

			if( position < sorted_current->StartPosition() ) {
				assert( sorted_current->BucketNo() > 0 );
				if( position < prev_bucket_position_start )
					throw InvalidStateException("Invalid prev bucket start");

				// this is allocated lazily, make sure it's here
				assert(prev_bucket_buf_);
				return prev_bucket_buf_.get() + (position - prev_bucket_position_start);
			}

			while( position >= sorted_current->EndPosition() )
				SwitchNextBucket(); // this will throw exception if beyond last bucket

//			if( !(this->final_position_end > position) ) // this should throw before on sorting not exists bucket...
//					throw InvalidValueException("Position too large");

//			if( position < sorted_current->StartPosition() ) // this seems impossible becasue we check this case before
//				throw InvalidValueException( "Position too small" );
			EnsureSortedTo( position + 1 /* at least byte of requested entry :) */ );
			return memory_start() + ( position - sorted_current->StartPosition() );
    }

		inline bool CloseToNewBucket(uint64_t position) const
    {
			if( !hasMoreBuckets ) return false;// no more backets

			return position > sorted_current->EndPosition()
						 || position + prev_bucket_buf_size / 2 >= sorted_current->EndPosition(); // TODO simplify to less work by saving proximity position in sorted_current
    }

		void TriggerNewBucket( uint64_t position )
		{
			if( position < prev_bucket_position_start )
				throw InvalidValueException("Triggering bucket too early");


			if( !sorted_current )
				SwitchNextBucket();
			else {
				if( position > sorted_current->EndPosition() )
					throw InvalidValueException("Triggering bucket too late");

				// save some of the current bucket, to allow some reverse-tracking
				// in the reading pattern,
				// position is the first position that we need in the new array
				sorted_current->WaitForSortedTo( sorted_current->EndPosition() ); // wait for everithing is sorted
				uint64_t const cache_size = sorted_current->EndPosition() - position;
				if( !prev_bucket_buf_ ) // init on first use
					prev_bucket_buf_.reset( Util::NewSafeBuffer( prev_bucket_buf_size ) );

				memcpy( prev_bucket_buf_.get(),
							 memory_start() + position - sorted_current->StartPosition(),
							 cache_size );
				memset( prev_bucket_buf_.get() + cache_size, 0x00, this->prev_bucket_buf_size - cache_size );


				prev_bucket_position_start = position;

				if( hasMoreBuckets ) SwitchNextBucket();
			}

    }

		void FlushCache( bool isFull = false )
		{
			for (auto& b : buckets_)
				 b.Flush( isFull );
    }

		inline uint64_t BiggestBucketSize() const {
			uint64_t res = 0;
			for( auto &b : buckets_ )
				if( b.Size() > res ) res = b.Size();
			return res;
		}

    ~SortManager()
    {
			if( sorted_current )
				memory_manager.release( sorted_current->buffer_size * (sorted_next? 2 : 1) ); // release current and next bucket's memory

			// Close and delete files in case we exit without doing the sort
			for (auto& b : buckets_)
				b.Remove( ); // remove buckets files
    }

private:

    // Size of the whole memory array
		MemoryManager &memory_manager;
    // Size of each entry
		const uint16_t entry_size_;
    // Bucket determined by the first "log_num_buckets" bits starting at "begin_bits"
		const uint32_t begin_bits_;
    // Log of the number of buckets; num bits to use to determine bucket
		const uint32_t log_num_buckets_;

		std::vector<SortingBucket> buckets_;
		bool hasMoreBuckets = true;

		const uint64_t prev_bucket_buf_size;
    std::unique_ptr<uint8_t[]> prev_bucket_buf_;
    uint64_t prev_bucket_position_start = 0;

		std::unique_ptr<SortedBucketBuffer> sorted_current, sorted_next;

		uint32_t num_threads, num_background_threads, num_read_threads;
		const uint8_t subbucket_bits;

		const uint8_t k_, phase_, table_index_;
		const uint32_t stats_mask;
		uint64_t time_total_wait = 0;

		uint64_t stream_read_position = 0;

		const inline uint8_t* memory_start() const {return sorted_current->buffer();}

		bool StartFirstBucket(){
			if( sorted_current ) return false;

			if( buckets_.size() == 0 )
				throw InvalidStateException( "No buckets to sort" );

			// find biggerst bucket and reserv ram to sort it
			uint64_t reserved_buffer_size  = BiggestBucketSize();

			if( !memory_manager.request(  reserved_buffer_size, true ) ) {
				std::cout<< "!!! Not enough memory for sort in RAM. Need to sort " <<
						(reserved_buffer_size / (1024.0 * 1024.0 * 1024.0)) << "GiB!!!" << std::endl;
				std::cout<< "!!!Going to continue using bigger buffer than in parameters!!!" << std::endl;

				memory_manager.requier( reserved_buffer_size );
			}

			if( num_threads > 1 ){
				if( memory_manager.request( reserved_buffer_size, true ) )
					sorted_next.reset( new SortedBucketBuffer( reserved_buffer_size, true ) ); // reserve for next bucket background sorting
				else
					std::cout << "Warning buffer is too small to allow background sorting, the performance may be slower." << std::endl
										<< "It is possible to increase buffer or to add number of buckets to improve multithreaded performance." << std::endl;
			}

			sorted_current.reset( new SortedBucketBuffer( reserved_buffer_size, !!sorted_next ) ); // reserve for current bucket
			sorted_current->StartSorting( 0, 0, &buckets_[0], num_background_threads, num_read_threads );

			return true;
		}
		// this done without switching to next bucket
		void EnsureSortedTo( uint64_t position = 0 ){
			StartFirstBucket(); // to ensure sorting started

			sorted_current->WaitForSortedTo( position );

			TrySortNextBucket();
		}

		inline void TrySortNextBucket() {
			if( hasMoreBuckets && sorted_next && sorted_next->BucketNo() == -1 && !sorted_current->isReading() ){
				uint next_bucket_no = sorted_current->BucketNo() + 1;
				if( next_bucket_no < buckets_.size() )
					sorted_next->TryStartSorting( next_bucket_no, sorted_current->EndPosition(),
																		&buckets_[next_bucket_no], num_background_threads, num_read_threads );
				else hasMoreBuckets = false;
			}
		}

		void ShowStatistics(){
			double const total_ram = memory_manager.getTotalSize() / (1024.0 * 1024.0 * 1024.0);
			double const cache_ram = memory_manager.getAccessibleRam() / (1024.0 * 1024.0 * 1024.0);
			double const free_ram = memory_manager.getFreeRam() / (1024.0 * 1024.0 * 1024.0);
			double const qs_ram = sorted_current->Bucket()->Size() / (1024.0 * 1024.0 * 1024.0);

			std::cout << "\r\tk" << (uint32_t)k_ << " p" << (uint32_t)phase_ << " t" << (uint32_t)table_index_
								<< " Bucket " << sorted_current->BucketNo() << " size: "
								<< std::setprecision( qs_ram > 10 ? 1:( qs_ram>1? 2 : 3) ) <<  qs_ram << "GiB,"
								<< " ram: " << std::fixed
								<< std::setprecision( total_ram > 10 ? 1:( total_ram>1? 2 : 3) ) << total_ram << "GiB";
			if( memory_manager.CacheEnabled )
				std::cout << std::fixed << ", cache: "
								<< std::setprecision( cache_ram > 10 ? 1:( cache_ram>1? 2 : 3) ) << cache_ram << "GiB";
			std::cout << std::fixed << ", free: "
								<< std::setprecision( free_ram > 10 ? 1:( free_ram>1? 2 : 3) ) << free_ram << "GiB"
								<< std::flush;

			std::cout << ", times: ( read:" << (sorted_current->Bucket()->read_time)/1000.0 << "s, total: "
								<< (sorted_current->Bucket()->sort_time)/1000.0 << "s )" << std::flush;

			//if( !hasMoreBuckets )
				std::cout << std::endl;
		}

		void SwitchNextBucket()
		{
			if( StartFirstBucket() ) return; // for first bucket

			if( !hasMoreBuckets )
				throw InvalidValueException( "Trying to sort bucket which does not exist." );

			ShowStatistics();

			if( sorted_next && sorted_next->BucketNo() == sorted_current->BucketNo() + 1){
				sorted_current->GetFrom( *sorted_next );
				TrySortNextBucket();
			} else {
				uint next_bucket_no = sorted_current->BucketNo() + 1;
				sorted_current->StartSortingNext( &buckets_[next_bucket_no], num_background_threads, num_read_threads );
				hasMoreBuckets = next_bucket_no + 1 < buckets_.size();
			}
//			if( next_bucket_to_sort >= buckets_.size() )
//				throw InvalidValueException( "Trying to sort bucket which does not exist." );

//			if( next_bucket_to_sort == 0 ) PrepareToSort();


//			SortingBucket& b = buckets_[next_bucket_to_sort];


//			double const total_ram = memory_manager.getTotalSize() / (1024.0 * 1024.0 * 1024.0);
//			double const cache_ram = memory_manager.getAccessibleRam() / (1024.0 * 1024.0 * 1024.0);
//			double const free_ram = memory_manager.getFreeRam() / (1024.0 * 1024.0 * 1024.0);
//			double const qs_ram = b.Size() / (1024.0 * 1024.0 * 1024.0);

//			std::cout << "\r\tk" << (uint32_t)k_ << " p" << (uint32_t)phase_ << " t" << (uint32_t)table_index_
//								<< " Bucket " << next_bucket_to_sort << " size: "
//								<< std::setprecision( qs_ram > 10 ? 1:( qs_ram>1? 2 : 3) ) <<  qs_ram << "GiB,"
//								<< " ram: " << std::fixed
//								<< std::setprecision( total_ram > 10 ? 1:( total_ram>1? 2 : 3) ) << total_ram << "GiB";
//			if( memory_manager.CacheEnabled )
//				std::cout << std::fixed << ", cache: "
//								<< std::setprecision( cache_ram > 10 ? 1:( cache_ram>1? 2 : 3) ) << cache_ram << "GiB";
//			std::cout << std::fixed << ", free: "
//								<< std::setprecision( free_ram > 10 ? 1:( free_ram>1? 2 : 3) ) << free_ram << "GiB"
//								<< std::flush;


//			int64_t wait_time = 0;
//			if( next_bucket_sorting_thread ){
//				std::cout << ", thr: (bg: " << num_background_treads << ", read: " << num_read_threads << ")" << std::flush;
//				auto start_time = std::chrono::high_resolution_clock::now();
//				next_bucket_sorting_thread->join(); // wait for sorting finish
//				auto end_time = std::chrono::high_resolution_clock::now();
//				next_bucket_sorting_thread.reset(); // free resources

//				wait_time = (end_time - start_time)/std::chrono::milliseconds(1);
//				time_total_wait += wait_time;
//				if( wait_time == 0 )
//					wait_time = b.sort_time - (start_time-b.start_time)/std::chrono::milliseconds(1); // negative value

//				int64_t time_threashold = b.sort_time/20; // 5%

//				// if we are waiting for sorts than we can add threads to it
//				if( wait_time < 0 && (-wait_time) > time_threashold ){
//					if( num_background_treads > 2 && (num_background_treads*2) > num_read_threads ){
//							num_background_treads--; // if no wait time than adjust number of computations threads down
//					}
//					else if( num_read_threads > 2 )
//						num_read_threads--; // if number of computational threads cant be down down read thredds
//				}
//				else if( wait_time > time_threashold ){
//					if( num_background_treads < num_threads )
//						num_background_treads++;

//					if( num_read_threads < (num_threads<<2) && (b.getInstantReads()<<3/*12.5%*/) > b.getTotalReads() )
//						num_read_threads++; // have space to improve reading time.
//				}
//				sorted_current_bucket.swap( sorted_next_bucket ); // put sorted result in action
//			}
//			else{
//				b.SortToMemory( sorted_current_bucket.get(), num_threads, num_read_threads );
//			}

//			if( b.Size() > 0 )
//				std::cout << ", times: ( wait:" << wait_time/1000.0 << "s, read:" << b.read_time/1000.0 << "s, total: "
//									<< b.sort_time/1000.0 << "s )" << std::flush;

//			this->final_position_start = this->final_position_end;
//			this->final_position_end += b.Size();
//			this->next_bucket_to_sort += 1;

//			if( this->next_bucket_to_sort >= buckets_.size()
//					|| ( buckets_[this->next_bucket_to_sort].Size() == 0 && b.Size() > 0 ) ){
//				// the final bucket - show some stats
//				uint64_t read_time = 0, sort_time = 0;
//				for(uint32_t i = 0; i < next_bucket_to_sort; i++ ){
//					read_time += buckets_[i].read_time;
//					sort_time += buckets_[i].sort_time;
//				}
//				std::cout << std::endl << "\tk" << (uint32_t)k_ << " p" << (uint32_t)phase_ << " t" << (uint32_t)table_index_
//									<< " avg times:( wait: " << time_total_wait/1000.0/next_bucket_to_sort << "s, read: "
//									<< read_time/1000.0/next_bucket_to_sort << "s, total sort:"
//									<< sort_time/1000.0/next_bucket_to_sort << "s )" << std::endl;
//			}

//			if( sorted_next_bucket ){ // reserved sort for next bucket in case it is threaded
//				if( this->next_bucket_to_sort < buckets_.size() )
//					next_bucket_sorting_thread.reset(
//								new std::thread( [this](){ buckets_[this->next_bucket_to_sort].SortToMemory(
//																					sorted_next_bucket.get(), num_background_treads, num_read_threads );} ) );
//				else { // no next sort -> we can release the next sort buffer
//					sorted_next_bucket.reset();
//					memory_manager.release( reserved_buffer_size );
//				}
//			}

//			b.Remove(); // Deletes the bucket file
		}
};

#endif  // SRC_CPP_FAST_SORT_ON_DISK_HPP_
