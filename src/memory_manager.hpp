// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//    http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef SRC_CPP_CACHE_MANAGER_HPP_
#define SRC_CPP_CACHE_MANAGER_HPP_

#include <mutex>
#include <assert.h>
#include <thread>
#include <chrono>
#include <vector>
#include "util.hpp"

extern uint64_t BUF_SIZE;

using namespace std::chrono_literals; // for operator""ms;

struct ICacheConsumer{
	// this function should return true if everithing is freed
	virtual bool FreeCache( int64_t size_to_free ) = 0;
};

struct MemoryManager{
	const bool CacheEnabled;
	const int64_t total_size;

	MemoryManager( uint64_t size, bool withCache = true )
		: CacheEnabled(withCache), total_size(size) {}

	inline uint64_t getAccessibleRam() const {
		return  total_size - used_ram + cleanable_ram;
	}

	inline int64_t getFreeRam() const {
		return total_size - used_ram;
	}

	inline int64_t getNotAllocatedRam() const {
		return total_size - used_ram - regular_buffers.size() * BUF_SIZE + BUF_SIZE;
	}

	inline int64_t getInUseRam() const {
		return used_ram;
	}

	inline int64_t getNotWritten() const {
		return not_written;
	}

	void SetMode( bool isForceClean, bool isFIFO ){
		this->isFIFO = isFIFO;
		if( isForceClean != isBackgroundClean ){

			if( background_cleaner ){
				isBackgroundClean = false;
				background_cleaner->join(); // wait for previous thread to finish
				background_cleaner.reset();
			}

			isBackgroundClean = isForceClean;

			if( CacheEnabled && isBackgroundClean ){
				// start background cleaner thread
				background_cleaner.reset( new std::thread( [this]{
					while( this->isBackgroundClean ){
						std::this_thread::sleep_for( 100ms );
						if( isNeedClean )
							isNeedClean = !CleanOne();
					}
				}) );
			}
		}
	}

	inline bool request( const uint64_t &size, bool forced = false ){

		if(getFreeRam() >= (int64_t)size || (forced && CleanCache( size )) ){
			{
				std::scoped_lock lk (sync_size);
				used_ram += size;
			}
			FreeBuffers( maxStoredBuffers() );
			return true;
		}
		if( isBackgroundClean ) isNeedClean = true;

		return false;
	}

	inline void requier( const uint64_t & size ){
		CleanCache( size );
		std::scoped_lock lk(sync_size);
		used_ram += size;
	}

	inline void release( const uint32_t &size ){
		std::scoped_lock lk(sync_size);

		assert( (int64_t)size <= used_ram );
		used_ram -= size;
	}

	inline int32_t registerConsumer( ICacheConsumer * consumer ){
		if( !CacheEnabled ) return -1; // disabled caching
		std::scoped_lock lk ( sync_consumers );
		if( min_consumer_idx > 0 ){
			consumers[--min_consumer_idx] = consumer;
			return min_consumer_idx;
		}
		consumers.push_back( consumer );
		return consumers.size()-1;
	}

	inline void unregisterConsumer( ICacheConsumer * consumer, uint32_t idx ){
		if( idx >= consumers.size() || consumers[idx] != consumer ) return; // check before lock to prevent deadlocks
		std::scoped_lock lk ( sync_consumers );
		if( idx >= min_consumer_idx && idx < consumers.size() && consumers[idx] == consumer ){
			consumers[idx] = nullptr;
			if( idx == min_consumer_idx ) min_consumer_idx++;
		}
	}

	inline uint8_t* consumerRequest(){
		uint8_t* res = nullptr;
		{
			std::lock_guard<std::mutex> lk(sync_buffers);
			if( regular_buffers.size() > 0 ){
				res = regular_buffers[regular_buffers.size()-1];
				regular_buffers.pop_back();
			}
		}

		if( res == nullptr && getFreeRam() >= (int64_t)BUF_SIZE )
			res = Util::NewSafeBuffer( BUF_SIZE );

		if( res != nullptr ){
			std::scoped_lock lk (sync_size);
			used_ram += BUF_SIZE;
			cleanable_ram += BUF_SIZE;
		} else	if( isBackgroundClean ) isNeedClean = true;

		return res;
	}

	inline void consumerRelease( uint8_t* buffer, uint32_t not_written = 0 ){
		{
			std::scoped_lock lk (sync_size);
			used_ram -= BUF_SIZE;
			cleanable_ram -= BUF_SIZE;
			not_written += not_written;
		}
		if( buffer != nullptr )
		{
			std::lock_guard<std::mutex> lk(sync_buffers);
			if( regular_buffers.size() < maxStoredBuffers() )
				regular_buffers.push_back( buffer );
			else delete[] buffer;
		}
	}


	~MemoryManager(){
		if( background_cleaner ){
			isBackgroundClean = false;
			background_cleaner->join();
		}
		FreeBuffers( 0 );
	}
private:

	int64_t used_ram = 0, cleanable_ram = 0, not_written = 0;
	std::mutex sync_size, sync_consumers, sync_buffers;
	std::vector<ICacheConsumer*> consumers;
	uint32_t min_consumer_idx = 0;
	bool isBackgroundClean = false;
	bool isFIFO = false;
	bool isNeedClean = false;
	std::unique_ptr<std::thread> background_cleaner;

	std::vector<uint8_t*> regular_buffers;

	inline uint32_t maxStoredBuffers() const {
		return std::max( 1UL, (total_size-used_ram)/BUF_SIZE );
	}

	inline void FreeBuffers( uint32_t max_to_leave ){
		std::lock_guard<std::mutex> lk(sync_buffers);
		for( ; regular_buffers.size() > max_to_leave; regular_buffers.pop_back() )
			delete[] regular_buffers[regular_buffers.size()-1];
	}

	inline bool CleanCache( int64_t need_size ){

		while( CleanOne() && getFreeRam() < need_size )
			;

		return getFreeRam() >= need_size;
	}

	inline bool CleanOne(){
		ICacheConsumer * cur = nullptr;
		uint32_t idx;

		{	// find consumer to clean
			std::lock_guard<std::mutex> lk(sync_consumers);
			if( isFIFO ){ // first in first out -> remove oldest
				while( cur == nullptr && min_consumer_idx < consumers.size() ){
					cur = consumers[ idx = min_consumer_idx];
					consumers[min_consumer_idx++] = nullptr;
				}
			}
			else{
				for( int32_t i = consumers.size()-1; cur == nullptr && i >= (int32_t)min_consumer_idx; i-- ){
					cur = consumers[ idx = i ];
					consumers[i] = nullptr;
				}
			}
		}
		if( cur != nullptr ) {
			cur->FreeCache( total_size );
			return true;
		}

		return false;
	}

};


#endif // SRC_CPP_CACHE_MANAGER_HPP_
