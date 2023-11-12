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

const inline uint32_t CacheBucketSize = 256; // mesured in number of entries
const inline uint32_t CacheBucketSizeLimit = 250; // mesured in number of entries

// Small bucket used in thread writings
struct CacheBucket{
	explicit CacheBucket( SortingBucket &cacheFor )
			: entries( CacheBucketSize*cacheFor.EntrySize() )
			, entry_size( cacheFor.EntrySize() )
			, limit_size( CacheBucketSizeLimit*cacheFor.EntrySize() )
			, parent( cacheFor )
	{	}

	inline void Add( const uint8_t *entry, const uint32_t &stats ){
		statistics[entries.used()/entry_size] = stats;
		entries.add( entry, entry_size );

		if( entries.used() > limit_size ){
			if( parent.TryAddEntriesTS( entries, statistics ) )
				entries.setUsed( 0 );
			else if( entries.isFull() ) {
				parent.AddEntriesTS( entries, statistics );
				entries.setUsed( 0 );
			}
		}
	}

	~CacheBucket(){
		if( entries.used() > 0 )
			parent.AddEntriesTS( entries, statistics );
	}

private:
	uint32_t statistics[CacheBucketSize];
	StreamBuffer entries;
	const uint16_t entry_size;
	const uint32_t limit_size;
	SortingBucket &parent;
};

class SortManager;
struct SortedBucketBuffer;
// this function used to show statistics from thread when sorting is finished
void ShowStats( const SortManager* sort_mngr, const SortedBucketBuffer *sbuf );

// ------------------------------------------------------------------
struct SortedBucketBuffer{
	const uint64_t buffer_size;
	const SortManager *sort_manager;
	int num_background_threads, num_read_threads;

	SortedBucketBuffer( const SortManager *sort_manager, uint64_t buffer_size, std::mutex *read_mutex, int num_background_threads, int num_read_threads )
			: buffer_size( buffer_size )
			, sort_manager( sort_manager )
			, num_background_threads(num_background_threads)
			, num_read_threads( num_read_threads )
			, bucket_buffer( Util::NewSafeBuffer( buffer_size ) )
			, read_mutex( read_mutex )
	{
	}

	~SortedBucketBuffer(){ FinishSort(); }

	inline uint64_t WiatsCount() const { return waits_count; }
	inline uint64_t WiatsReadCount() const { return read_waits_count; }
	inline double WiatsTimeSec() const { return wait_time; }

	// returns true if current sort just finished.
	inline void WaitForSortedTo( uint64_t position ){
		assert( position >= start_position && position <= end_position );

		if( sorting_thread.load( std::memory_order_relaxed ) != NULL && position >= start_position && position <= end_position ){
			waits_count++;
			if( bucket->SortedPosision() <= 0 ) read_waits_count++;

			auto start_time = std::chrono::high_resolution_clock::now();

			if( position == end_position ) FinishSort();
			else {
				while( sorting_thread.load( std::memory_order_relaxed ) != NULL && position >= start_position
							 && (position - start_position) > bucket->SortedPosision() ){
					std::this_thread::sleep_for( 100us );
				}

				if( bucket->SortedPosision() >= bucket->Size() )
					FinishSort();
			}

			wait_time += (std::chrono::high_resolution_clock::now() - start_time)/std::chrono::milliseconds(1)/1000.0;
		}
	}

	inline bool WaitForSorted(){ return FinishSort(); }

	inline const uint8_t* buffer() const { return bucket_buffer.get(); }
	inline SortingBucket * Bucket() const { return bucket; }
	inline int32_t BucketNo() const { return bucket_no; }
	inline uint64_t StartPosition() const { return start_position; }
	inline uint64_t EndPosition() const { return end_position; }

	// Sort next bucket synchronously
	void SortBucket( SortingBucket * next_bucket ){
		bucket_no++;
		bucket = next_bucket;

		start_position = end_position;
		end_position = start_position + next_bucket->Size();

		bucket->SortToMemory( bucket_buffer.get(), num_background_threads, num_read_threads );
		ShowStats( this->sort_manager, this );
	}

	void StartSorting( uint bucket_no, uint64_t start_position, SortingBucket * bucket ){
		assert( bucket_no >= 0 );

		this->bucket_no = (int)bucket_no;
		this->bucket = bucket;

		this->start_position = start_position;
		this->end_position = start_position + bucket->Size();

		waits_count = 0;
		FinishSort( new std::thread(
										[this]( SortingBucket* bucket, uint8_t* buf, std::mutex *r_mutex ){
												r_mutex->lock();
												bucket->SortToMemory( buf, num_background_threads, num_read_threads, r_mutex );
												ShowStats( this->sort_manager, this );
												// adjust number of threads for next sorting?
												if( read_waits_count > 0 && num_read_threads < num_background_threads )
													num_read_threads++;
											}, this->bucket, this->bucket_buffer.get(), this->read_mutex ) );
	}
private:
	std::unique_ptr<uint8_t[]> bucket_buffer;
	int bucket_no = -1;
	uint64_t start_position = 0, end_position = 0;
	SortingBucket *bucket = NULL;
	std::atomic<std::thread*> sorting_thread = NULL;
	std::mutex *read_mutex;
	uint64_t waits_count = 0, read_waits_count = 0;
	double wait_time = 0;


	inline bool FinishSort( std::thread * newThread = NULL){
		auto prev = sorting_thread.exchange( newThread, std::memory_order_relaxed );

		if( prev == NULL ) return false;

		prev->join();
		delete prev;
		return true;
	}
};

// ====================================================
class SortManager : public Disk, IReadDiskStream {
public:
	SortManager(
			MemoryManager &memory_manager,
			uint32_t  num_buckets,
			uint32_t log_num_buckets_not_used,
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
			, prev_bucket_buf_size(
					2 * (stripe_size + 10 * (kBC / pow(2, kExtraBits))) * entry_size)
			, num_threads( num_threads )
			, k_(k), phase_(phase), table_index_(table_index)
	{

		uint64_t sorting_size = (1ULL<<k)*entry_size;
		double expected_buckets_no = num_buckets*( phase_==1 ?1.0:( phase_==3?0.64:0.8 ) );
		const auto expected_bucket_size = sorting_size/(double)expected_buckets_no;
		const auto memory_size = memory_manager.getTotalSize()/(isSingleSort()?1:2);
		const double buckets_in_ram = num_threads > 1 ? 2.1 : 1.05;

//		std::cout << "Estimation by k: " << (int)k << ", entry_size: " << entry_size
//							<< ", size: " << sorting_size << ", buckets_no: " << expected_buckets_no
//							<< ", bucket_size: " << (expected_bucket_size/1024/1024) << "MiB, memory size: " << memory_size
//							<< ", num buckets in ram: " << buckets_in_ram << std::endl;

		if( memory_size/expected_bucket_size < buckets_in_ram ){
			uint32_t need_buckets = num_buckets;
			for( auto size = expected_bucket_size; memory_size/size < buckets_in_ram; size /= 2 )
				need_buckets *= 2;
			std::cout << std::setprecision(2)
								<< "Warning! Expected bucket size " << (expected_bucket_size/1024/1024) << "MiB for table "
								<< (uint32_t)table_index_ << " with " << num_buckets
								<< " buckets. Available buffer " << (memory_size>>20) << "MiB is too small for this size."
								<< " Increase number of buckets to " << need_buckets << std::endl;
			num_buckets = need_buckets;
		}

		log_num_buckets_ = log2(num_buckets);
		subbucket_bits = std::min( (uint8_t)(32-log_num_buckets_), std::max( (uint8_t)2, (uint8_t)(k - log_num_buckets_ - kSubBucketBits) ) );
		stats_mask = ( (uint64_t)1<<subbucket_bits)-1;

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
																							entry_size, begin_bits_ + log_num_buckets_, subbucket_bits,
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

		// region Disk inheritance implementaion {
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
//			for (auto& b : buckets_)
//				b.FreeMemory();

			prev_bucket_buf_.reset();

			if( sorted_current ){
				memory_manager.release( sorted_current->buffer_size * (sorted_next?2:1) ); // release current and next bucket's memory
				sorted_current.reset();
				sorted_next.reset();
			}
    }
		// endregion }

		// region IReadStream implementation {
		uint32_t Read( std::unique_ptr<uint8_t[]> &buf, const uint32_t &buf_size ) override {
			EnsureSortingStarted();

			for( uint32_t buf_start = 0; buf_start < buf_size; ){
				assert( stream_read_position <= sorted_current->EndPosition() );

				if( hasMoreBuckets && stream_read_position == sorted_current->EndPosition() )
					SwitchNextBucket();

				assert( stream_read_position <= sorted_current->EndPosition() );

				int64_t to_copy = std::min( sorted_current->EndPosition() - stream_read_position, (uint64_t)(buf_size-buf_start) );
				assert( to_copy >= 0 );

				if( to_copy == 0 ) {
					assert( checkSort( buf.get(), buf_start ) );
					return buf_start; // is it possible case?
				}

				assert( (stream_read_position + to_copy) <= sorted_current->EndPosition() );
				assert( (stream_read_position + to_copy) >= sorted_current->StartPosition() );

				sorted_current->WaitForSortedTo( stream_read_position + to_copy );

				assert(  (stream_read_position + to_copy) <= ( sorted_current->Bucket()->SortedPosision() +  sorted_current->StartPosition() ) );
				assert( sorted_current->EndPosition() - sorted_current->StartPosition() <= sorted_current->buffer_size );

				memcpy( buf.get() + buf_start, memory_start() + stream_read_position - sorted_current->StartPosition(), to_copy );
				buf_start += to_copy;
				stream_read_position += to_copy;
			}

			assert( checkSort( buf.get(), buf_size ) );

			return buf_size;
		};

		bool checkSort( uint8_t const *buf, uint64_t buf_size ) const {
			for( uint64_t i = entry_size_; i < buf_size; i += entry_size_ )
				if( Util::MemCmpBits( buf + i - entry_size_, buf + i, entry_size_, begin_bits_ ) > 0 ){
					return false;
				}
			return true;
		}

		bool atEnd() const override{
			return !hasMoreBuckets && sorted_current && (uint32_t)sorted_current->BucketNo() >= buckets_.size()
						 && stream_read_position >= sorted_current->EndPosition();
		};
		void Close() override { FreeMemory(); };

		inline uint64_t GetReadPosition() const { return stream_read_position; }
		// endregion IReadStream implementation }

		const uint8_t *ReadEntry( uint64_t position )
    {
			EnsureSortingStarted();
			assert( sorted_current );

			if( position < sorted_current->StartPosition() ) {
				assert( sorted_current->BucketNo() > 0 );
				if( position < prev_bucket_position_start )
					throw InvalidStateException("Invalid prev bucket start");

				// this is allocated lazily, make sure it's here
				assert(prev_bucket_buf_);
				return prev_bucket_buf_.get() + (position - prev_bucket_position_start);
			}

			while( position >= sorted_current->EndPosition() ){
				SwitchNextBucket(); // this will throw exception if beyond last bucket
				assert( position < sorted_current->EndPosition() );
			}

			sorted_current->WaitForSortedTo( position + 1 );

			assert( position >= sorted_current->StartPosition() && position < sorted_current->EndPosition() );

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
				StartSorting();
			else {
				if( position > sorted_current->EndPosition() )
					throw InvalidValueException("Triggering bucket too late");

				// save some of the current bucket, to allow some reverse-tracking
				// in the reading pattern,
				// position is the first position that we need in the new array
				sorted_current->WaitForSorted(); // wait for everithing is sorted
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
				b.Remove(); // remove buckets files
    }

private:

	// Size of the whole memory array
	MemoryManager &memory_manager;
	// Size of each entry
	const uint16_t entry_size_;
	// Bucket determined by the first "log_num_buckets" bits starting at "begin_bits"
	const uint32_t begin_bits_;
	// Log of the number of buckets; num bits to use to determine bucket
	uint32_t log_num_buckets_;

	std::vector<SortingBucket> buckets_;
	bool hasMoreBuckets = true;

	const uint64_t prev_bucket_buf_size;
	std::unique_ptr<uint8_t[]> prev_bucket_buf_;
	uint64_t prev_bucket_position_start = 0;

	std::unique_ptr<SortedBucketBuffer> sorted_current, sorted_next;
	std::mutex bucket_read_mutex;

	uint32_t num_threads;
	uint8_t subbucket_bits;

	const uint8_t k_, phase_, table_index_;
	uint32_t stats_mask;
	//uint64_t time_total_wait = 0;
	std::chrono::time_point<std::chrono::high_resolution_clock> start_sorting_time;

	uint64_t stream_read_position = 0;

	const inline uint8_t* memory_start() const {return sorted_current->buffer();}

	void inline EnsureSortingStarted() { if( !sorted_current ) StartSorting(); }

	bool StartSorting(){
		if( sorted_current ) return false;

		assert( buckets_.size() > 0 );

		start_sorting_time = std::chrono::high_resolution_clock::now();

		// find biggerst bucket and reserv ram to sort it
		uint64_t reserved_buffer_size  = BiggestBucketSize();

		if( !memory_manager.request(  reserved_buffer_size, true ) ) {
			std::cout<< "!!! Not enough memory for sort in RAM. Need to sort " <<
					(reserved_buffer_size / (1024.0 * 1024.0 * 1024.0)) << "GiB!!!" << std::endl;
			std::cout<< "!!!Going to continue using bigger buffer than in parameters!!!" << std::endl;

			memory_manager.requier( reserved_buffer_size );
		}

		uint32_t num_background_threads = std::max( num_threads>1?2U:1U, num_threads/( isSingleSort() ? 1 : 2) )
				, num_read_threads = std::max( num_threads>1?2U:1U, num_threads/( isSingleSort() ? 2 : 3)  );

		sorted_current.reset( new SortedBucketBuffer( this, reserved_buffer_size, &bucket_read_mutex, num_background_threads, num_read_threads ) ); // reserve for current bucket
		sorted_current->SortBucket( &buckets_[0] ); // first one always sync sort


		if( num_threads > 1 && buckets_.size() > 1 ){
			if( memory_manager.request( reserved_buffer_size, true ) ){
				sorted_next.reset( new SortedBucketBuffer( this, reserved_buffer_size, &bucket_read_mutex, num_background_threads, 1 /* num_read_threads */) ); // reserve for next bucket background sorting
				sorted_next->StartSorting( 1, sorted_current->EndPosition(), &buckets_[1] ); // backgound sorting for next bucket
			}
			else
				std::cout << "Warning buffer is too small to allow background sorting, the performance may be slower." << std::endl
									<< "It is possible to increase buffer or to add number of buckets to improve multithreaded performance." << std::endl;
		}

		return true;
	}

	// checks if this sort is single in total process
	inline bool isSingleSort() const { return phase_ == 1 || phase_ == 3 || (phase_ == 4 && table_index_ >= 6); }

	void ShowStatistics( const SortedBucketBuffer *stats_of ) const {
		double const total_ram = memory_manager.getTotalSize() / (1024.0 * 1024.0 * 1024.0);
		// double const cache_ram = memory_manager.getAccessibleRam() / (1024.0 * 1024.0 * 1024.0);
		double const free_ram = memory_manager.getFreeRam() / (1024.0 * 1024.0 * 1024.0);
		double const qs_ram = stats_of->Bucket()->Size() / (1024.0 * 1024.0 * 1024.0);
		int max_bucket = buckets_.size();
		while( max_bucket > 0 && buckets_[max_bucket-1].Count() == 0 ) max_bucket--;
		std::cout << "\r\tk" << (uint32_t)k_ << " p" << (uint32_t)phase_ << " t" << (uint32_t)table_index_
							<< " Bucket " << (stats_of->BucketNo()+1) << "/" << max_bucket << " size: "
							<< std::setprecision( qs_ram > 10 ? 1:( qs_ram>1? 2 : 3) ) <<  qs_ram << "GiB,"
							<< " ram: " << std::fixed
							<< std::setprecision( total_ram > 10 ? 1:( total_ram>1? 2 : 3) ) << total_ram << "GiB";
//			if( memory_manager.CacheEnabled )
//				std::cout << std::fixed << ", cache: "
//								<< std::setprecision( cache_ram > 10 ? 1:( cache_ram>1? 2 : 3) ) << cache_ram << "GiB";
		std::cout << std::fixed << ", free: "
							<< std::setprecision( free_ram > 10 ? 1:( free_ram>1? 2 : 3) ) << free_ram << "GiB"
							<< std::flush;

		auto showTime = [](double time){
			bool isMin = time > 600;
			std::cout << std::setprecision( isMin ? 0 : 2 ) << ( time /(isMin ? 60 : 1) )
								<< (isMin ? "m":"s" );
		};

		auto read_time = (stats_of->Bucket()->read_time)/1000.0;
		auto total_time = (stats_of->Bucket()->sort_time)/1000.0;
		auto passed_time = (std::chrono::high_resolution_clock::now() - start_sorting_time)/std::chrono::milliseconds(1)/1000.0;
		auto estimated_time = ( stats_of->BucketNo() == 0 ? passed_time : (passed_time-total_time)/stats_of->BucketNo() )
																													 * max_bucket - passed_time;

		std::cout << ", thr: (r: " << stats_of->num_read_threads << ", srt: " << stats_of->num_background_threads << ")";
		std::cout << std::setprecision( std::min( read_time, total_time ) < 10 ? 2 : 1 )
							<< ", times: ( read: " << read_time << "s, bucket: "
							<< total_time << "s, total: ";
		showTime( passed_time );
		std::cout <<", est. left: ";
		showTime( estimated_time );
		std::cout << " )";
		if( stats_of->WiatsCount() > 0 )
			std::cout << ", waits: (cnt: " << stats_of->WiatsCount() << ", r_cnt: " << stats_of->WiatsReadCount() << ", time: "
								<< std::setprecision(2) << stats_of->WiatsTimeSec() << "s)";
		std::cout << std::flush;

#ifdef NDEBUG
		if( !hasMoreBuckets || buckets_[stats_of->BucketNo()+1].Size() == 0 )
#endif
			std::cout << std::endl;
	}

	void SwitchNextBucket()
	{
		if( StartSorting() ) return; // for first bucket

		if( !hasMoreBuckets )
			throw InvalidValueException( "Trying to sort bucket which does not exist." );


		if( sorted_next ){
			assert( sorted_next->BucketNo() == sorted_current->BucketNo() + 1 );
			assert( !sorted_current->WaitForSorted() ); // sort on current should be finished at this point

			sorted_current.swap( sorted_next ); // bring sorted to front


			uint next_bucket_no = sorted_current->BucketNo() + 1;
			hasMoreBuckets = next_bucket_no < buckets_.size() && buckets_[next_bucket_no].Count() > 0;
			if( hasMoreBuckets )
				sorted_next->StartSorting( next_bucket_no, sorted_current->EndPosition(), &buckets_[next_bucket_no] );
		} else {
			uint next_bucket_no = sorted_current->BucketNo() + 1;
			sorted_current->SortBucket( &buckets_[next_bucket_no] );
			hasMoreBuckets = next_bucket_no + 1 < buckets_.size() && buckets_[next_bucket_no+1].Count() > 0;
		}

		if( !hasMoreBuckets ) std::cout << std::endl;
	}

	friend void ShowStats( const SortManager* sort_mngr, const SortedBucketBuffer *sbuf );
};

void ShowStats( const SortManager* sort_mngr, const SortedBucketBuffer *sbuf ){
	sort_mngr->ShowStatistics( sbuf );
}

#endif  // SRC_CPP_FAST_SORT_ON_DISK_HPP_
