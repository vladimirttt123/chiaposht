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

#ifndef SRC_CPP_UTIL_HPP_
#define SRC_CPP_UTIL_HPP_

#include <atomic>
#include <cassert>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#ifndef _WIN32
#include <sys/resource.h>
#endif

template <typename Int, typename Int2>
constexpr inline Int cdiv(Int a, Int2 b) { return (a + b - 1) / b; }

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <processthreadsapi.h>
#include "uint128_t.h"
#else
// __uint__128_t is only available in 64 bit architectures and on certain
// compilers.
typedef __uint128_t uint128_t;

// Allows printing of uint128_t
std::ostream &operator<<(std::ostream &strm, uint128_t const &v)
{
    strm << "uint128(" << (uint64_t)(v >> 64) << "," << (uint64_t)(v & (((uint128_t)1 << 64) - 1))
         << ")";
    return strm;
}

#endif

// compiler-specific byte swap macros.
#if defined(_MSC_VER)

#include <cstdlib>

// https://docs.microsoft.com/en-us/cpp/c-runtime-library/reference/byteswap-uint64-byteswap-ulong-byteswap-ushort?view=msvc-160
inline uint16_t bswap_16(uint16_t x) { return _byteswap_ushort(x); }
inline uint32_t bswap_32(uint32_t x) { return _byteswap_ulong(x); }
inline uint64_t bswap_64(uint64_t x) { return _byteswap_uint64(x); }

#elif defined(__clang__) || defined(__GNUC__)

inline uint16_t bswap_16(uint16_t x) { return __builtin_bswap16(x); }
inline uint32_t bswap_32(uint32_t x) { return __builtin_bswap32(x); }
inline uint64_t bswap_64(uint64_t x) { return __builtin_bswap64(x); }

#else
#error "unknown compiler, don't know how to swap bytes"
#endif

/* Platform-specific cpuid include. */
#if defined(_WIN32)
#include <intrin.h>
#elif defined(__x86_64__)
#include <cpuid.h>
#endif

uint64_t GetTotalBytesWritten();

class Timer {
public:
    Timer()
    {
				write_byte_start = GetTotalBytesWritten();
        wall_clock_time_start_ = std::chrono::steady_clock::now();
#if _WIN32
        ::GetProcessTimes(::GetCurrentProcess(), &ft_[3], &ft_[2], &ft_[1], &ft_[0]);
#else
        cpu_time_start_ = clock();
#endif
    }

    static char *GetNow()
    {
        auto now = std::chrono::system_clock::now();
        auto tt = std::chrono::system_clock::to_time_t(now);
        return ctime(&tt);  // ctime includes newline
    }

    void PrintElapsed(const std::string &name) const
    {
        auto end = std::chrono::steady_clock::now();
        auto wall_clock_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                                 end - this->wall_clock_time_start_)
                                 .count();

#if _WIN32
        FILETIME nowft_[6];
        nowft_[0] = ft_[0];
        nowft_[1] = ft_[1];

        ::GetProcessTimes(::GetCurrentProcess(), &nowft_[5], &nowft_[4], &nowft_[3], &nowft_[2]);
        ULARGE_INTEGER u[4];
        for (size_t i = 0; i < 4; ++i) {
            u[i].LowPart = nowft_[i].dwLowDateTime;
            u[i].HighPart = nowft_[i].dwHighDateTime;
        }
        double user = (u[2].QuadPart - u[0].QuadPart) / 10000.0;
        double kernel = (u[3].QuadPart - u[1].QuadPart) / 10000.0;
        double cpu_time_ms = user + kernel;
#else
        double cpu_time_ms =
            1000.0 * (static_cast<double>(clock()) - this->cpu_time_start_) / CLOCKS_PER_SEC;
#endif

        double cpu_ratio = static_cast<int>(10000 * (cpu_time_ms / wall_clock_ms)) / 100.0;

        std::cout << name << " " << (wall_clock_ms / 1000.0) << " seconds. CPU (" << cpu_ratio
									<< "%), written " << ( GetTotalBytesWritten() - write_byte_start)/1024.0/1024/1024 << " GiB. "
									<< Timer::GetNow();
    }

private:
    std::chrono::time_point<std::chrono::steady_clock> wall_clock_time_start_;
#if _WIN32
    FILETIME ft_[4];
#else
    clock_t cpu_time_start_;
#endif
		uint64_t write_byte_start;
};





#define MEM_SAFE_BUF_SIZE  7

#ifndef __GNUC__
#define NO_HUGE_PAGES
#endif

#ifndef NO_HUGE_PAGES
#include <linux/mman.h> // MAP_HUGE_1GB
#include <sys/mman.h> // mmap, munmap
const uint64_t HUGE_MEM_PAGE_BITS = 21;
const uint64_t HUGE_MEM_PAGE_SIZE = 1UL << HUGE_MEM_PAGE_BITS;
const uint64_t HUGE_1GB_PAGE_BITS = 30;
const uint64_t HUGE_1GB_PAGE_SIZE = 1UL << HUGE_1GB_PAGE_BITS;
#else
#define HUGE_MEM_PAGE_SIZE BUF_SIZE
#endif // NO_HUGE_PAGES

namespace Util {
	enum class AllocationType : uint8_t { HUGE_1G, HUGE_2M_INSTEADOF_1G, HUGE_2M, THP_INSTEADOF_1G, THP, NORMAL, UNDEFINED };

	struct {
		std::atomic_uint64_t current_usage = 0, top_usage = 0,
				current_usage_2M = 0, top_usage_2M = 0,
				possible_usage_2M = 0, top_possible_usage_2M = 0,
				current_usage_1G = 0, top_usage_1G = 0,
				possible_usage_1G = 0, top_possible_usage_1G = 0;
		void change( AllocationType type, int64_t size ){
			switch (type) {
				case AllocationType::NORMAL:
					setMax( current_usage.fetch_add( size, std::memory_order::relaxed ) + size, top_usage );
					break;
				case AllocationType::HUGE_1G:
					setMax( current_usage_1G.fetch_add( size, std::memory_order_relaxed ) + size, top_usage_1G );
					break;
				case AllocationType::HUGE_2M:
					setMax( current_usage_2M.fetch_add( size, std::memory_order_relaxed ) + size, top_usage_2M );
					break;
				case AllocationType::HUGE_2M_INSTEADOF_1G:
					setMax( current_usage_2M.fetch_add( size, std::memory_order_relaxed ) + size, top_usage_2M );
					setMax( possible_usage_1G.fetch_add( size, std::memory_order_relaxed ) + size + current_usage_1G.load(std::memory_order::relaxed), top_possible_usage_1G );
					break;
				case AllocationType::THP_INSTEADOF_1G:
					setMax( current_usage.fetch_add( size, std::memory_order_relaxed ) + size, top_usage );
					setMax( possible_usage_1G.fetch_add( size, std::memory_order_relaxed ) + size + current_usage_1G.load(std::memory_order::relaxed), top_possible_usage_1G );
					break;
				case AllocationType::THP:
					setMax( current_usage.fetch_add( size, std::memory_order_relaxed ) + size, top_usage );
					setMax( possible_usage_2M.fetch_add( size, std::memory_order_relaxed ) + size + current_usage_2M.load(std::memory_order::relaxed), top_possible_usage_2M );
					break;
				case AllocationType::UNDEFINED: break;
			}
		}

		inline void setMax( uint64_t cur, std::atomic_uint64_t&to ){
			assert( ((int64_t)cur) >= 0 );

			uint64_t exp = to.load(std::memory_order::relaxed);
			assert( ((int64_t)exp) >= 0 );
			while( exp < cur && !to.compare_exchange_weak( exp, cur, std::memory_order::relaxed, std::memory_order::relaxed ) )
				assert( ((int64_t)exp) >= 0 );
		}

		void print( bool withClear = true ){
			std::cout << "Memory alloctions stats: top: " << top_usage
#ifndef NO_HUGE_PAGES
								<< "; top 1G pages: " << (top_usage_1G>>HUGE_1GB_PAGE_BITS);
			if( top_possible_usage_1G > 0 )
					std::cout << "; possible 1G pages: " << (top_possible_usage_1G>>HUGE_1GB_PAGE_BITS);
			std::cout << "; top 2M pages: " << (top_usage_2M>>HUGE_MEM_PAGE_BITS);
			if( top_possible_usage_2M > 0 )
					std::cout << "; possible 2M pages: " << (top_possible_usage_2M>>HUGE_MEM_PAGE_BITS);
			std::cout
#endif
								<< std::endl;
			if( withClear ) clear();
		}

		void clear(){
			//current_usage = current_usage_2M = current_usage_1G =
			top_usage = top_usage_2M =
				possible_usage_2M = top_possible_usage_2M =
				top_usage_1G = possible_usage_1G = top_possible_usage_1G = 0;
		}
	} MemAllocationStats;


	template <typename T>
	struct Deleter{
#ifdef NO_HUGE_PAGES
		Deleter( AllocationType t = AllocationType::UNDEFINED ){};
		void operator()(T *p){ delete [] p; }
		void swap(Deleter<T>& other){}
#else // NO_HUGE_PAGES

		Deleter( AllocationType type = AllocationType::UNDEFINED, uint64_t size = 0 )
				: size(size), type(type) {
			MemAllocationStats.change( type, size );
		}

		void operator()(T *p){
			if( size == 0 || p == NULL ) return;
			MemAllocationStats.change( type, -size );
			switch( type ) {
				case AllocationType::HUGE_1G:
				case AllocationType::HUGE_2M_INSTEADOF_1G:
				case AllocationType::HUGE_2M:
					//std::cout << "unmap" << std::endl;
					munmap( (void*)p, size );
					break;
				case AllocationType::THP_INSTEADOF_1G:
				case AllocationType::THP:
					//std::cout << "free" << std::endl;
					std::free( (void*)p );
					break;
				case AllocationType::NORMAL:
					//std::cout << "delete" << std::endl;
					delete [] p;
					break;
				case AllocationType::UNDEFINED: break;
			}
		}
		void swap(Deleter<T>& other)
		{
			auto t = other.type; other.type = type; type = t;
			auto s = other.size; other.size = size; size = s;
		}
	private:
		uint64_t size;
		AllocationType type;
#endif // else NO_HUGE_PAGES
	};
	template <typename T>
	inline std::unique_ptr<T, Deleter<T>> allocate( uint64_t count, double huge_page_fraction = 0.8, double huge_1gb_page_min_use = 0.9 ){
		if( count == 0 )
			return std::unique_ptr<T, Deleter<T>>( NULL, Deleter<T>() );

#ifndef NO_HUGE_PAGES
		auto size =  sizeof(T)*count;
		if( size >= HUGE_MEM_PAGE_SIZE*huge_page_fraction ){
			// Try to allocate huge pages
			uint64_t h1gb_size = (( size + HUGE_1GB_PAGE_SIZE - 1 )>>HUGE_1GB_PAGE_BITS)<<HUGE_1GB_PAGE_BITS;
			bool suits_for_1G =  size >= HUGE_1GB_PAGE_SIZE*huge_1gb_page_min_use; // (size / (double)h1gb_size ) > huge_1gb_page_min_use;
			if( suits_for_1G ) {
				auto ptr = mmap( NULL, h1gb_size , PROT_READ | PROT_WRITE,
												MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB | MAP_HUGE_1GB,
												-1, 0);
				if( ptr != MAP_FAILED )
					return std::unique_ptr<T, Deleter<T>>( (T*)ptr, Deleter<T>( AllocationType::HUGE_1G, h1gb_size ) );

				std::cout << "Cannot allocate 1GiB " << (h1gb_size >> HUGE_1GB_PAGE_BITS) << " huge pages ( currently allocated " << MemAllocationStats.current_usage_1G << " pages): " <<  errno << " " <<::strerror(errno) << std::endl;
			}

			uint64_t hsize = ( size <= HUGE_MEM_PAGE_SIZE ? 1 : ((size + HUGE_MEM_PAGE_SIZE-1)>>HUGE_MEM_PAGE_BITS) ) << HUGE_MEM_PAGE_BITS;
			auto ptr = mmap( NULL, hsize , PROT_READ | PROT_WRITE,
											MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB | MAP_HUGE_2MB,
											-1, 0);
			if( ptr != MAP_FAILED )
				return 	std::unique_ptr<T, Deleter<T>>( (T*)ptr, Deleter<T>( suits_for_1G ? AllocationType::HUGE_2M_INSTEADOF_1G : AllocationType::HUGE_2M, hsize ) );

			std::cout << "Cannot allocate " << (hsize >> HUGE_MEM_PAGE_BITS) << " huge pages " << (suits_for_1G?" instead of 1GiB ":"")
								<< " ( currently allocated " << MemAllocationStats.current_usage_2M << " pages ): " <<  errno << " " <<::strerror(errno) << std::endl;

			// Try to allocate THP
			if( posix_memalign( &ptr, HUGE_MEM_PAGE_SIZE, hsize ) == 0 && ptr != nullptr ){
				madvise( ptr, hsize, MADV_HUGEPAGE );
				return std::unique_ptr<T, Deleter<T>>( (T*)ptr, Deleter<T>( suits_for_1G ? AllocationType::THP_INSTEADOF_1G : AllocationType::THP, hsize ) );
			}
		}
#endif // NO_HUGE_PAGES

		// std:: cout << "allocate regular page" << std::endl;
		return std::unique_ptr<T, Deleter<T>>( new T[count], Deleter<T>( AllocationType::NORMAL, sizeof(T)*count ));
	}

	// Converting and extracting bits of last element of buffer
	// usually needs additional 7 bytes at the end.
	// this funciton is adds and inits this bytes.
	inline uint8_t * NewSafeBuffer( const uint64_t &size ){
		auto res = new uint8_t[size+8];
		// init additional ram - this remove valgrind warning of accessing not inited and saving pointer there allow check buffers
		((uint64_t*)(res + size))[0] = reinterpret_cast<std::uint64_t>(res);
		return res;
	}

    template <typename X>
		inline X Mod(X i, X n) { return (i % n + n) % n; }

    inline uint32_t ByteAlign(uint32_t num_bits) { return (num_bits + (8 - ((num_bits) % 8)) % 8); }

    inline std::string HexStr(const uint8_t *data, size_t len)
    {
        std::stringstream s;
        s << std::hex;
        for (size_t i = 0; i < len; ++i)
            s << std::setw(2) << std::setfill('0') << static_cast<int>(data[i]);
        s << std::dec;
        return s.str();
    }

    inline void IntToTwoBytes(uint8_t *result, const uint16_t input)
    {
        uint16_t r = bswap_16(input);
        memcpy(result, &r, sizeof(r));
    }

    // Used to encode deltas object size
    inline void IntToTwoBytesLE(uint8_t *result, const uint16_t input)
    {
        result[0] = input & 0xff;
        result[1] = input >> 8;
    }

    inline uint16_t TwoBytesToInt(const uint8_t *bytes)
    {
        uint16_t i;
        memcpy(&i, bytes, sizeof(i));
        return bswap_16(i);
    }

		/* Converts a 64 bit int to bytes. */
    inline void IntToEightBytes(uint8_t *result, const uint64_t input)
    {
			*((uint64_t*)result) = bswap_64( input );
    }

		/* Converts a byte array to a 64 bit int. */
    inline uint64_t EightBytesToInt(const uint8_t *bytes)
    {
			return bswap_64( *((const uint64_t*)bytes) );
    }

		inline void IntTo16Bytes(uint8_t *result, const uint128_t input)
    {
			// WARNING this implementation for BIGENDIANS only!!!
			((uint64_t*)result)[0] = bswap_64( ((const uint64_t*)&input)[1] );
			((uint64_t*)result)[1] = bswap_64( ((const uint64_t*)&input)[0] );
		}

    /*
     * Retrieves the size of an integer, in Bits.
     */
    inline uint8_t GetSizeBits(uint128_t value)
    {
        uint8_t count = 0;
        while (value) {
            count++;
            value >>= 1;
        }
        return count;
    }



    // 'bytes' points to a big-endian 64 bit value (possibly truncated, if
    // (start_bit % 8 + num_bits > 64)). Returns the integer that starts at
    // 'start_bit' that is 'num_bits' long (as a native-endian integer).
    //
    // Note: requires that 8 bytes after the first sliced byte are addressable
    // (regardless of 'num_bits'). In practice it can be ensured by allocating
    // extra 7 bytes to all memory buffers passed to this function.
    inline uint64_t SliceInt64FromBytes(
        const uint8_t *bytes,
        uint32_t start_bit,
        const uint32_t num_bits)
    {
        uint64_t tmp;

        if (start_bit + num_bits > 64) {
            bytes += start_bit / 8;
            start_bit %= 8;
        }

        tmp = Util::EightBytesToInt(bytes);
        tmp <<= start_bit;
        tmp >>= 64 - num_bits;
        return tmp;
    }

		inline uint64_t SliceInt64FromBytes(
				const uint8_t *bytes,
				const uint32_t &num_bits)
		{
			return Util::EightBytesToInt( bytes ) >> (64-num_bits);
		}

    inline uint64_t SliceInt64FromBytesFull(
        const uint8_t *bytes,
        uint32_t start_bit,
        uint32_t num_bits)
    {
        uint32_t last_bit = start_bit + num_bits;
        uint64_t r = SliceInt64FromBytes(bytes, start_bit, num_bits);
        if (start_bit % 8 + num_bits > 64)
            r |= bytes[last_bit / 8] >> (8 - last_bit % 8);
        return r;
    }

    inline uint128_t SliceInt128FromBytes(
        const uint8_t *bytes,
        const uint32_t start_bit,
        const uint32_t num_bits)
    {
        if (num_bits <= 64)
            return SliceInt64FromBytesFull(bytes, start_bit, num_bits);

        uint32_t num_bits_high = num_bits - 64;
        uint64_t high = SliceInt64FromBytesFull(bytes, start_bit, num_bits_high);
        uint64_t low = SliceInt64FromBytesFull(bytes, start_bit + num_bits_high, 64);
        return ((uint128_t)high << 64) | low;
    }

    inline void GetRandomBytes(uint8_t *buf, uint32_t num_bytes)
    {
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_int_distribution<int> dist(0, 255);
        for (uint32_t i = 0; i < num_bytes; i++) {
            buf[i] = dist(mt);
        }
    }

    inline uint64_t ExtractNum(
        const uint8_t *bytes,
        uint32_t len_bytes,
        uint32_t begin_bits,
        uint32_t take_bits)
    {
        if ((begin_bits + take_bits) / 8 > len_bytes - 1) {
            take_bits = len_bytes * 8 - begin_bits;
        }
        return Util::SliceInt64FromBytes(bytes, begin_bits, take_bits);
    }

		//Warning! this function can access beyond the buffer because it needs at least 8 bytes after begin_bits/8 also when extracting only 16bit
		inline uint64_t ExtractNum64( const uint8_t *bytes, const uint32_t begin_bits, const uint32_t take_bits ){
			assert( (begin_bits&7) + take_bits <= 64 );
			auto moved = bswap_64( ((uint64_t*)(bytes + (begin_bits>>3) ))[0] ) >> (64-take_bits - (begin_bits&7));
			auto mask = (((uint64_t)1)<<take_bits)-1;
			return moved&mask;
		}

		// extracts exact 32 bit - access at most 5 bytes if begin_bits > 0 and only 4 when begin_bits == 0
		// begin bits should be less than 8
		inline uint32_t ExtractNum32( const uint8_t *bytes, const uint32_t begin_bits ){
			assert( begin_bits < 8 );

			return begin_bits == 0 ? bswap_32( ((const uint32_t*)bytes)[0] )
					: ( (bswap_32( ((const uint32_t*)bytes)[0] ) << begin_bits ) | (bytes[4]>>(8-begin_bits)) );
		}

    // The number of memory entries required to do the custom SortInMemory algorithm, given the
    // total number of entries to be sorted.
    inline uint64_t RoundSize(uint64_t size)
    {
        size *= 2;
        uint64_t result = 1;
        while (result < size) result *= 2;
        return result + 50;
    }

    /*
     * Like memcmp, but only compares starting at a certain bit.
     */
    inline int MemCmpBits(
				const uint8_t *left_arr,
				const uint8_t *right_arr,
        uint32_t len,
        uint32_t bits_begin)
    {
	if( bits_begin&7 ){
    	    uint32_t start_byte = bits_begin >> 3;
	    uint8_t l = (left_arr[start_byte] << (bits_begin&7)), r = ( right_arr[start_byte] << (bits_begin&7) );
            if( l != r ) return l - r;
        }
        //int rn = l - r;
        //if( rn ) return rn;
//        uint8_t mask = ((1 << (8 - (bits_begin % 8))) - 1);
//        if ((left_arr[start_byte] & mask) != (right_arr[start_byte] & mask)) {
//            return (left_arr[start_byte] & mask) - (right_arr[start_byte] & mask);
//        }
        
        for (uint32_t i = (bits_begin>>3) + ((bits_begin&7)?1:0); i < len; i++) {
//        for (uint32_t i = start_byte + 1; i < len; i++) {
//    	    int re = left_arr[i] - right_arr[i];
//	    if( re ) return re;
            if (left_arr[i] != right_arr[i])
                return left_arr[i] - right_arr[i];
        }
        return 0;
    }

    inline double RoundPow2(double a)
    {
        // https://stackoverflow.com/questions/54611562/truncate-float-to-nearest-power-of-2-in-c-performance
        int exp;
        double frac = frexp(a, &exp);
        if (frac > 0.0)
            frac = 0.5;
        else if (frac < 0.0)
            frac = -0.5;
        double b = ldexp(frac, exp);
        return b;
    }

#if defined(_WIN32) || defined(__x86_64__)
    void CpuID(uint32_t leaf, uint32_t *regs)
    {
#if defined(_WIN32)
        __cpuid((int *)regs, (int)leaf);
#else
        __get_cpuid(leaf, &regs[0], &regs[1], &regs[2], &regs[3]);
#endif /* defined(_WIN32) */
    }

    bool HavePopcnt(void)
    {
        // EAX, EBX, ECX, EDX
        uint32_t regs[4] = {0};

        CpuID(1, regs);
        // Bit 23 of ECX indicates POPCNT instruction support
        return (regs[2] >> 23) & 1;
    }
#endif /* defined(_WIN32) || defined(__x86_64__) */

    inline uint64_t PopCount(uint64_t n)
    {
#if defined(_WIN32)
        return __popcnt64(n);
#elif defined(__x86_64__)
        uint64_t r;
        __asm__("popcnt %1, %0" : "=r"(r) : "r"(n));
        return r;
#else
        return __builtin_popcountl(n);
#endif /* defined(_WIN32) ... defined(__x86_64__) */
    }


		void setOpenFilesLimit( uint32_t forNumberOfBuckets ){
#ifndef _WIN32
			// Increases the open files limit, in case it is too low.
			struct rlimit the_limit;// = { need_limit , need_limit };
			if( getrlimit(RLIMIT_NOFILE, &the_limit ) < 0 )
					std::cout << "Warning: cannot read files limit... skipping" << std::endl;
			else{
				rlim_t need_limit = 100 + forNumberOfBuckets*3;
				if( the_limit.rlim_cur < need_limit ){
					if( the_limit.rlim_max < need_limit ){
						std::cout << "Warning: max open files limit " << the_limit.rlim_max << " is less than sugested " << need_limit << std::endl;
						need_limit = the_limit.rlim_max;
					}
					the_limit.rlim_cur = need_limit;
					if( -1 == setrlimit(RLIMIT_NOFILE, &the_limit) ) {
						std::cout << "Warning: set opened files limit failed" << std::endl;
					}
				}
			}
#endif
		}

		// used for debug
		uint32_t CheckSort( uint8_t *memory, uint32_t entry_len, uint32_t const bits_begin, uint32_t entries_number ){
			for( uint32_t i = entry_len; i < entries_number*entry_len; i+=entry_len ){
				if( Util::MemCmpBits( memory + i - entry_len, memory + i, entry_len, bits_begin) > 0)
					return i/entry_len;

				//if( Util::MemCmpBits( memory + i - entry_len, memory + i, entry_len, bits_begin) == 0 ) std::cout << "equals!!" << std::endl;
			}

			return 0;
		}
} // end of namespac Util

#endif  // SRC_CPP_UTIL_HPP_
