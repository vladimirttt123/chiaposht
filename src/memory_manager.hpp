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
#include <vector>

//enum cache_type : uint8_t {
//	BITFIELD = 1,
//	SORT_BUFFER = 2,
//	BUCKET_STREAM = 3,
//	TABLE_FILE = 4
//};

enum cache_mode : uint8_t {
	SKIP_NEW = 1,
	FLUSH_OLD = 2
};

struct ICacheConsumer{
	virtual uint64_t getUsedCache() const = 0;
	virtual void FreeCache() = 0;
};


struct MemoryManager{
	bool isForced = false;
	bool isFIFO = false;

	MemoryManager( uint64_t size ) : total_size(size){}

	inline uint64_t getAccessibleRam() const {
		uint64_t res = total_size - consumed;

		for( int32_t i = consumers.size() - 1; i >= (int32_t)consumers_start; i-- )
			if( consumers[i] != nullptr )
				res += consumers[i]->getUsedCache();

		return res;
	}

	inline uint64_t getFreeRam() const {
		return total_size - consumed;
	}

	inline bool request( const uint64_t &size, ICacheConsumer *consumer ){
		return request( size, false, consumer );
	}

	inline bool request( const uint64_t &size, bool forced = false, ICacheConsumer *consumer = nullptr ){

		{
			std::scoped_lock lk (syncMutex );

			if( total_size - consumed >= size ){
				consumed += size;
				return true;
			}

			// check if forced and we can clean enough
			if( !(forced || isForced )
					|| ( getAccessibleRam() - ( consumer == nullptr ? 0 : consumer->getUsedCache() ) ) < size  )
				return false;
		}

		// no lock here to allow releasing
		while( releaseConsumer( consumer ) && total_size - consumed < size )
			/* empty body */;

		{ // second check for space
			std::scoped_lock lk (syncMutex );

			if( total_size - consumed >= size ){
				consumed += size;
				return true;
			}
		}
		return false;
	}

	inline void release( const uint32_t &size ){
		std::scoped_lock lk (syncMutex );
		assert( (int64_t)size <= consumed );
		consumed -= size;
	}

	inline uint32_t registerConsumer( ICacheConsumer * consumer ){
		std::scoped_lock lk (syncMutex );
		consumers.push_back( consumer );
		return consumers.size()-1;
	}

	inline void unregisterConsumer( ICacheConsumer * cleaner, uint32_t idx ){
		std::scoped_lock lk (syncMutex );
		if( idx < consumers.size() && consumers[idx] == cleaner )
			consumers[idx] = nullptr;
		else{ // TODO update idx in consumer...
			for( uint32_t i = consumers_start; i < consumers.size(); i++ )
				if( consumers[i] == cleaner ){
					consumers[i] = nullptr;
					return;
				}
		}
	}

	void reorginizeConsumers(){
		std::scoped_lock lk (syncMutex );
		uint32_t cur = 0;
		for( uint32_t i = consumers_start; i < consumers.size(); i++ ){
			if( consumers[i] != nullptr )
				consumers[cur++] = consumers[i];
		}
		consumers_start = 0;
		consumers.resize(cur);
	}


private:
	const uint64_t total_size;
	int64_t consumed = 0;
	std::mutex syncMutex;
	std::vector<ICacheConsumer*> consumers;
	uint32_t consumers_start = 0;


	inline bool releaseConsumer( ICacheConsumer * cur_consumer ){

		ICacheConsumer * consumer;

		{
			std::scoped_lock lk (syncMutex );
			if( consumers.size() <= consumers_start ) return false;

			if( isFIFO ){
				consumer = consumers.back();
				if( consumer == cur_consumer ) return false;
				consumers.pop_back();
			}
			else {
				consumer = consumers[consumers_start];
				if( consumer == cur_consumer ) return false;
				consumers_start++;
			}
		}


		if( consumer != nullptr )
			consumer->FreeCache();

		return true;
	}

};


//struct MemoryLocker{
//	MemoryLocker( MemoryManager memory_manager, uint32_t lock_size = 0 )
//		:memory_manager(memory_manager)
//	{
//		if( lock_size ) {

//		}
//	}

//	lock
//private:
//	MemoryManager &memory_manager;
//	uint64_t locked_size;

//};

#endif // SRC_CPP_CACHE_MANAGER_HPP_
