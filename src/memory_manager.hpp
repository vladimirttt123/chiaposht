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

#include <atomic>
#include <mutex>
#include <assert.h>
#include <thread>
#include <chrono>
#include <vector>
#include "util.hpp"

extern uint64_t BUF_SIZE;

using namespace std::chrono_literals; // for operator""ms;

enum cache_free_status : uint8_t{
	NO_CLEAN = 0,
	PARTIAL_CLEAN = 1,
	FULL_CLEAN = 2
};

struct ICacheConsumer{
	// this function should return true if everithing is freed
	virtual cache_free_status FreeCache( int64_t size_to_free ) = 0;
};


struct BuffersStack{
	const uint64_t max_buffers;

	BuffersStack( uint64_t max_buffers )
			: max_buffers(max_buffers), buffers_size( max_buffers + 250UL /*it could be threads count or should be more than treads count*/ )
			, buffers( new std::atomic<uint8_t*>[buffers_size] ),	min_idx(0), max_idx(0) {
		for( uint64_t i = 0; i < buffers_size; i++ )
			buffers[i].store( nullptr, std::memory_order_relaxed );
	}

	inline uint64_t size() const { return std::max(0L, (int64_t)max_idx.load(std::memory_order_relaxed) - (int64_t)min_idx.load(std::memory_order_relaxed)); }

	inline void put( uint8_t* buf ){
		if( buf == nullptr ) return;
		if( size() <= max_buffers ) {
			buf = buffers[ max_idx.fetch_add( 1, std::memory_order_relaxed) % buffers_size ].exchange( buf, std::memory_order_relaxed );
			assert( buf == nullptr ); // the place should be empty
		}

		if( buf != nullptr )
			delete [] buf;
	}

	inline uint8_t* get(){
		auto idx = min_idx.load( std::memory_order_relaxed );
		while( idx < max_idx ){
			if( min_idx.compare_exchange_strong( idx, idx+1, std::memory_order_relaxed) )
				return buffers[idx%buffers_size].exchange( nullptr, std::memory_order_relaxed );

			idx = min_idx.load( std::memory_order_relaxed );
		}
		return nullptr;
	}

	~BuffersStack(){
		for( uint32_t i = 0; i < buffers_size; i++ )
			if( buffers[i] != nullptr)
				delete [] buffers[i];
	}
private:
	const uint64_t buffers_size;
	std::unique_ptr<std::atomic<uint8_t*>[]> buffers;
	std::atomic_uint64_t min_idx, max_idx;
};

struct MemoryManager{
	const bool CacheEnabled;

	MemoryManager( uint64_t size = 0, int64_t max_cache_size = 0 )
		: CacheEnabled(max_cache_size > 0), total_size(size), max_cache_size(max_cache_size)
		, regular_buffers( std::max( 0L, max_cache_size/(int64_t)BUF_SIZE ) ){}

	inline uint64_t getTotalSize() const { return total_size; }

	inline uint64_t getAccessibleRam() const {
		return  total_size - used_ram + cleanable_ram;
	}

	inline int64_t getFreeRam() const {
		return total_size - used_ram;
	}
	inline int64_t getFreeCache() const {
		return max_cache_size - cleanable_ram;
	}

	inline int64_t getNotAllocatedRam() const {
		return total_size - used_ram - regular_buffers.size() * BUF_SIZE + BUF_SIZE;
	}

	inline int64_t getInUseRam() const { return used_ram;	}
	inline int64_t getNotWritten() const { return not_written; }

	void SetMode( bool isForceClean, bool isFIFO ){
		this->isFIFO = isFIFO;
		this->isForcedClean = isForceClean;
	}

	inline bool request( const uint64_t &size, bool forced = false ){

		if( getFreeRam() >= (int64_t)size || ( ( forced || isForcedClean ) && CleanCache( size )) ){
			used_ram += size;
			FreeBuffers( maxStoredBuffers() );
			return true;
		}

		return false;
	}

	inline void requier( const uint64_t & size ){
		CleanCache( size );
		used_ram += size;
	}

	inline void release( const uint64_t &size ){
		assert( size <= used_ram );
		used_ram -= size;
	}

	inline int32_t registerConsumer( ICacheConsumer * consumer ){
		if( !CacheEnabled ) return -1; // disabled caching
		std::lock_guard lk ( sync_consumers );
//	REUSING of old indexes is problematic...
//		if( min_consumer_idx > 0 ){
//			consumers[--min_consumer_idx] = consumer;
//			return min_consumer_idx;
//		}
		consumers.push_back( consumer );
		return consumers.size()-1;
	}

	inline void unregisterConsumer( ICacheConsumer * consumer, uint32_t idx ){
		if( idx >= consumers.size() || consumers[idx] != consumer ) return; // check before lock to prevent deadlocks
		std::lock_guard lk ( sync_consumers );
		if( idx >= min_consumer_idx && idx < consumers.size() && consumers[idx] == consumer ){
			consumers[idx] = nullptr;
			if( idx == min_consumer_idx ) min_consumer_idx++;
		}
	}

	inline uint8_t* consumerRequest(){
		uint8_t* res = regular_buffers.get();

		if( res == nullptr ){

			if( getFreeCache() >= (int64_t)BUF_SIZE ){
				res = Util::NewSafeBuffer( BUF_SIZE ); // have space for new buffer
			} else if( isForcedClean && CleanCache( BUF_SIZE ) ){
				res = regular_buffers.get();
			}
		}

		if( res != nullptr ){
			used_ram += BUF_SIZE;
			cleanable_ram += BUF_SIZE;
		}

		return res;
	}

	inline void consumerRelease( uint8_t* buffer, uint64_t cache_hit_size_size = 0 ){

		used_ram -= BUF_SIZE;
		cleanable_ram -= BUF_SIZE;
		not_written += cache_hit_size_size;

		if( buffer != nullptr )
			regular_buffers.put( buffer );
	}


	~MemoryManager(){
		CleanCache( total_size );
		FreeBuffers( 0 );
	}
private:
	const int64_t total_size, max_cache_size;


	std::atomic_ullong used_ram = 0, cleanable_ram = 0, not_written = 0;
	std::mutex sync_consumers;
	std::vector<ICacheConsumer*> consumers;
	uint32_t min_consumer_idx = 0;
	bool isFIFO = false;
	bool isForcedClean = false;

	BuffersStack regular_buffers;

	inline uint32_t maxStoredBuffers() const {
		return std::max( 1ULL, (total_size-used_ram)/BUF_SIZE );
	}

	inline void FreeBuffers( uint32_t max_to_leave ){
		while( regular_buffers.size() > max_to_leave ){
			auto buf = regular_buffers.get();
			assert( buf != nullptr );
			delete[]buf;
		}
	}

	inline bool CleanCache( int64_t need_size ){

		while( CleanOne( need_size - getFreeCache() ) && getFreeCache() < need_size )
			;

		return getFreeCache() >= need_size;
	}

	inline bool CleanOne( int64_t size_to_clean ){
		ICacheConsumer * cur = nullptr;
		uint32_t idx;

			// find consumer to clean
		if( isFIFO ){ // first in first out -> remove oldest
			if( min_consumer_idx < consumers.size() ) { // check to less locks
				std::lock_guard<std::mutex> lk(sync_consumers);
				while( cur == nullptr && min_consumer_idx < consumers.size() ){
					cur = consumers[ idx = min_consumer_idx];
					consumers[min_consumer_idx++] = nullptr;
				}
			}
		}
		else{
			std::lock_guard<std::mutex> lk(sync_consumers);
			for( int32_t i = consumers.size()-1; cur == nullptr && i >= (int32_t)min_consumer_idx; i-- ){
				cur = consumers[ idx = i ];
				consumers[i] = nullptr;
			}
		}

		if( cur != nullptr ) {
			auto clean_res = cur->FreeCache( size_to_clean );
			if( clean_res == NO_CLEAN )
				CleanOne( size_to_clean ); // clean another one

			if( clean_res != FULL_CLEAN ){ // not fully cleaned need to return the consumer back
				// no need of lock because there is no consumer index reuse
				consumers[idx] = cur;
				if( idx < min_consumer_idx ){
					std::lock_guard<std::mutex> lk(sync_consumers);
					if( idx < min_consumer_idx ) min_consumer_idx = idx;
				}
			}

			return true;
		}

		return false;
	}

};


#endif // SRC_CPP_CACHE_MANAGER_HPP_
