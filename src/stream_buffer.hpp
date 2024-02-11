// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//    http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef SRC_CPP_STREAM_BUFFER_HPP_
#define SRC_CPP_STREAM_BUFFER_HPP_

#include "disk.hpp"

struct StreamBuffer{

	StreamBuffer( const StreamBuffer & other ) = delete;

	explicit StreamBuffer( uint32_t size = BUF_SIZE ) : size_(size){	}

	StreamBuffer( uint8_t* buf, uint32_t size, uint32_t used = 0 )
		: buffer( buf ), size_(size), used_size( used )	{	}

	inline uint8_t* get(){
		return buffer == nullptr ?
			( buffer = Util::NewSafeBuffer( size_ ) ) : buffer;
	}

	inline uint8_t* getEnd(){ return get() + used_size; }

	inline StreamBuffer & add( const uint8_t* data, uint32_t data_size) {
		assert( used_size + data_size <= size_ );
		memcpy( get() + used_size, data, data_size );
		used_size += data_size;
		return *this;
	}
	inline StreamBuffer & add( const uint8_t byte ){
		assert( used_size + 1 <= size_ );
		get()[used_size] = byte;
		used_size++;
		return *this;
	}

	inline StreamBuffer & setUsed( uint32_t new_used ) {
		assert( new_used <= size_ );
		used_size = new_used;
		return *this;
	}
	inline StreamBuffer & addUsed( uint32_t add ) {
		assert( used_size + add <= size_ );
		used_size += add;
		return *this;
	}

	inline uint32_t used() const { return used_size; }
	inline uint32_t size() const { return size_; }
	inline bool isFull() const { return used_size >= size_; }

	inline operator bool() const { return buffer != nullptr; }


	inline StreamBuffer & ensureSize( uint32_t minSize, bool withCleanUsed = true ){
		if( minSize > size_ ){
			if( buffer != nullptr ){
				if( used_size > 0 && !withCleanUsed ){
					auto new_buf = Util::NewSafeBuffer( minSize );
					memcpy( new_buf, buffer, used_size );
					delete [] buffer;
					buffer = new_buf;
				} else {
					delete [] buffer;
					buffer = nullptr;
				}

			}
			size_ = minSize;
		}
		if( withCleanUsed ) used_size = 0;
		return *this;
	}

	inline StreamBuffer & swap( StreamBuffer & other ){
		auto s = other.size_;
		other.size_ = size_;
		size_ = s;

		s = other.used_size;
		other.used_size = used_size;
		used_size = s;

		auto buf = other.buffer;
		other.buffer = buffer;
		buffer = buf;

		return *this;
	}

	inline void reset( uint8_t* new_buf, uint32_t new_size, uint32_t used = 0 ){
		reset();
		if( new_buf != nullptr ){
			assert( new_size > 0 );
			assert( new_size >= used );
			buffer = new_buf;
			size_ = new_size;
			used_size = used;
		}
	}

	inline void reset(){
		if( buffer != nullptr ){
			delete[]buffer;
			buffer = nullptr;
			used_size = 0;
		}
	}

	// Warning replacement should be same size as inner buffer!!!
	inline uint8_t* release( uint8_t * replacement = nullptr){
		auto res = buffer;
		buffer = replacement;
		used_size = 0;
		return res;
	}

	~StreamBuffer(){
		if( buffer != nullptr ) delete[]buffer;
	}

private:
	uint8_t * buffer = nullptr;
	uint32_t size_;
	uint32_t used_size = 0;
};


#endif //SRC_CPP_STREAM_BUFFER_HPP_
