#ifndef SRC_CPP_CMP_TOOLS_HPP_
#define SRC_CPP_CMP_TOOLS_HPP_

#include <atomic>
#include <cstdint>
#include <cstring>
#ifndef _WIN32
#include <unistd.h>
#endif
#include <string>
#include <iostream>
#include <fstream>

#include "util.hpp"
#include "pos_constants.hpp"

namespace TCompress {

// This is replacement made for kFormatDescription in compress plots
const std::string tFormatDescription = "t0.1";
const char* plotMagicFrase = "Proof of Space Plot";

const char * program_header = "*** Chia plot compressing software made by Vladimir T\n"
			"*** If this compression is helpfull for you, please consider donate\n"
			"***   xch1ch6s3q0enuj9wtemn473gkkvj0u8vlggypr375mk547e7aa48hmsql74e8\n";


struct LinePointCacheEntry{
public:
	uint64_t partial_id = 0, position = 0;
	uint128_t line_point = 0;
	uint8_t importance = 0, k = 0, table_no = 0;
};

struct LinePointsCache{
public:
	LinePointsCache( uint32_t size = 4096 ) : size(size), cache( new LinePointCacheEntry[size] ){	}

	uint128_t AddLinePoint( uint8_t k, const uint8_t * id, uint8_t table_no, uint64_t position, uint128_t line_point, bool important = false ){
		if( line_point == 0 ) return 0;
		std::lock_guard<std::mutex> lk(mut);
		while( cache[next_pos].importance != 0 ){
			cache[next_pos].importance--;
			next_pos = (next_pos+1)%size;
		}

		cache[next_pos].k = k;
		cache[next_pos].table_no = table_no;
		cache[next_pos].position = position;
		cache[next_pos].line_point = line_point;
		cache[next_pos].importance = table_no > 1 ? 0 : ( important ? 3 : 1 );
		cache[next_pos++].partial_id = ((uint64_t*)id)[0];
		if( count < next_pos ) count = next_pos;
		next_pos %= size;

		return line_point;
	}

	uint128_t GetLinePoint( uint8_t k, const uint8_t * id, uint8_t table_no, uint64_t position ){
		std::lock_guard<std::mutex> lk(mut);
		for( uint32_t i = 1; i <= count; i++ ){
			uint32_t idx = (next_pos - i)%size;
			if( cache[idx].k == k && cache[idx].table_no == table_no
					&& cache[idx].position == position && cache[idx].partial_id == ((uint64_t*)id)[0] )
				return cache[idx].line_point;
		}
		return 0; // NOT FOUND
	}

	void IncreaseSize() {
		std::lock_guard<std::mutex> lk(mut);
		if( count < size ) return; // not a time to increase ( may be increase from other thread )

		uint32_t new_size = size * 2;
		LinePointCacheEntry *new_cache = new LinePointCacheEntry[new_size];
		memcpy( new_cache, cache, sizeof(LinePointCacheEntry)*count );
		next_pos = size;
		size = new_size;
		delete [] cache;
		cache = new_cache;
	}

	uint32_t GetCount() const { return count; }
	uint32_t GetSize() const { return size; }
	uint32_t GetImportantCount() const {
		uint32_t res = 0;
		for( uint32_t i = 0; i < count; i++ )
			if( cache[i].importance > 0 ) res++;
		return res;
	}

	~LinePointsCache(){ delete[] cache; }
private:
	uint32_t size;
	LinePointCacheEntry *cache;
	uint32_t next_pos = 0, count = 0;
	std::mutex mut;
};


struct DeltasStorage{
public:
	const uint32_t parks_count;

	uint64_t total_size = 0;
	uint16_t *all_sizes = nullptr;

	DeltasStorage( uint64_t parks_count ) : parks_count(parks_count){
		if( parks_count > 0xffffffffUL )
			throw std::runtime_error("too many parks " + std::to_string(parks_count) );

		all_sizes = new uint16_t[ parks_count ];
	}

	void Add( uint64_t idx, uint16_t size ){
		assert( idx < parks_count );
		total_size += (all_sizes[idx] = size);
	}

	void TotalEndToBuf( uint64_t idx, uint8_t *buf ){
		buf[0] = total_size >> ((idx&1)?24:16);
		buf[1] = total_size>>8;
		buf[2] = total_size;
	}

	// check all partially saved position could be restored
	bool IsDeltasPositionRestorable(){
		uint64_t initial_total_size = total_size;
		uint64_t park_avg_size = total_size/parks_count;
		uint8_t buf[6];
		uint8_t *buf_prev = buf, *buf_cur = buf + 3;
		uint64_t delta_pos;
		uint16_t delta_size;

		total_size = 0; // to work with function TotalEndToBuf

		for( uint64_t i = 0; i < parks_count; i++ ){
			total_size += all_sizes[i];
			TotalEndToBuf( i, buf_cur );
			RestoreParkPositionAndSize( park_avg_size, i, buf_prev, buf_cur, delta_pos, delta_size );
			if( delta_pos != (total_size-all_sizes[i]) || delta_size != all_sizes[i] ){
				std::cout << "Couldn't restore parkd index " << i << " with position " << (total_size-all_sizes[i])
									<< " and size " << all_sizes[i] << ". Predicted position " << delta_pos
									<< " predicted size " << delta_size << std::endl;

				total_size = initial_total_size;
				return false;
			}

			buf_prev = buf_cur;
			buf_cur = buf + 3*(i&1);
		}

		return true;
	}

	static void RestoreParkPositionAndSize( uint64_t park_avg_size, uint64_t park_idx, uint8_t *prev_buf, uint8_t *cur_buf,
																				 uint64_t &delta_position, uint16_t &encoded_delta_size ){
		uint64_t cur_pos = (((uint16_t)cur_buf[1])<<8) | cur_buf[2];

		if( park_idx == 0 ){
			delta_position = 0;
			encoded_delta_size = cur_pos;
			return;
		}


		uint64_t prev_pos = (((uint16_t)prev_buf[1])<<8) | prev_buf[2];
		encoded_delta_size = cur_pos - prev_pos;

		if( park_idx&1 ){
			prev_pos |= ((uint64_t)prev_buf[0])<<16;
			bool is_jump = (prev_pos+encoded_delta_size)>>24; // is adding current delta will lead to jump to next value in bits24+
			prev_pos |= (((uint64_t)(cur_buf[0]-(is_jump?1:0)) )&0xffUL)<<24;
		}else{
			cur_pos = (cur_pos | 0x1000000 | (((uint64_t)cur_buf[0])<<16)) - encoded_delta_size;
			prev_pos |= (cur_pos&0xff0000) | (((uint64_t)prev_buf[0])<<24);
		}

		uint64_t expected_by_magic = park_avg_size*park_idx;
		if( (expected_by_magic&0xffffffff) < 0x50000000	&& ( prev_pos&0xffffffff ) > 0xa0000000 )
			expected_by_magic -= 0x100000000UL;
		else if( (expected_by_magic&0xffffffff) > 0xa0000000	&& ( prev_pos&0xffffffff ) < 0x50000000 )
			expected_by_magic += 0x100000000UL;
		prev_pos |= expected_by_magic&~0xffffffffUL;

		delta_position = prev_pos;
	}

	~DeltasStorage(){
		if( all_sizes != nullptr ) delete[]all_sizes;
	}

private:
	static uint64_t compound( uint8_t top, uint8_t high, uint8_t med, uint8_t low )
	{
		return (((uint64_t)top)<<24) | (((uint64_t)high)<<16) | (((uint64_t)med)<<8)  | (((uint64_t)low));
	}
};


struct ReadFileWrapper{
public:
	explicit ReadFileWrapper( const std::string& filename )
			: disk_file( new std::ifstream( filename, std::ios::in | std::ios::binary ) ), isExternal(false)
	{
		if( !disk_file->is_open() )
			throw std::invalid_argument( "Cannot open file " + filename);
	}

	explicit ReadFileWrapper( std::ifstream * file )
			: disk_file( file ), isExternal(true)
	{
		if( !disk_file->is_open() )
			throw std::invalid_argument( "File is not opened" );
	}

	void Read( uint8_t * buffer, uint64_t size ){
		int64_t pos = disk_file->tellg();
		disk_file->read(reinterpret_cast<char*>(buffer), size);

		if (disk_file->fail()) {
			std::cout << "Failed to read input disk at position " << pos << " size " << size << std::endl;
			throw std::runtime_error("disk read file " +
															 std::to_string(size) + " at position " + std::to_string(pos));
		}
	}

	void Read( uint64_t position, uint8_t * buffer, uint64_t size ){
		Seek( position );
		Read( buffer, size );
	}

	void Seek( uint64_t position ){
		disk_file->seekg( position );

		if( disk_file->fail() ) {
			std::cout << "Disk seek FAILED to position " << position << std::endl;
			throw std::runtime_error("disk seek failed " + std::to_string(position));
		}
	}

	int64_t Size(){
		int64_t pos = disk_file->tellg();
		disk_file->seekg( 0, std::ios_base::end );
		if (disk_file->fail()) {
			std::cout << "Disk seek FAILED to end of file" << std::endl;
			throw std::runtime_error("disk seek failed to end of file" );
		}
		int64_t size = disk_file->tellg();
		Seek( pos );
		return size;
	}

	~ReadFileWrapper(){
		if( !isExternal ) delete disk_file;
	}
private:
	std::ifstream *disk_file;
	bool isExternal;
};


struct ParkReader{
public:
	const bool is_compressed;
	const uint64_t line_point_size, line_point_size_bits, stubs_size, stub_size_bits;
	const uint128_t first_line_point;
	uint16_t deltas_size;
	uint8_t *buf;
	std::vector<uint8_t> deltas;

	ParkReader( uint8_t *buf, uint64_t line_point_size_bits )
			: is_compressed(false), line_point_size(0), line_point_size_bits(line_point_size_bits)
			, stubs_size(0),  stub_size_bits(0)
			, first_line_point( Util::SliceInt128FromBytes( buf, 0, line_point_size_bits ) )
	{	}

	ParkReader( uint8_t *buf, uint64_t line_point_size_bits, uint64_t stub_size_bits, uint32_t max_deltas_size, uint8_t table_no )
			: is_compressed(false), line_point_size( (line_point_size_bits+7)/8 ), line_point_size_bits(line_point_size_bits)
			, stubs_size( Util::ByteAlign(stub_size_bits*(kEntriesPerPark-1))/8), stub_size_bits(stub_size_bits)
			, first_line_point( Util::SliceInt128FromBytes( buf, 0, line_point_size_bits ) )
			, deltas_size( ((uint32_t)buf[line_point_size+stubs_size]) | (((uint32_t)buf[line_point_size+stubs_size+ 1])<<8) )
			, buf(buf), cur_stub_buf( buf + line_point_size )

	{
		if( deltas_size&0x8000 )
			throw std::runtime_error( "uncompressed deltas is not supported yet" );
		if( deltas_size > max_deltas_size )
			throw std::runtime_error( "incorrect deltas size " + std::to_string( deltas_size ) );

		uint8_t *deltas_buf = buf + line_point_size + stubs_size + 2;
		deltas = Encoding::ANSDecodeDeltas( deltas_buf , deltas_size, kEntriesPerPark - 1, kRValues[table_no-1] );
	}

	ParkReader( uint8_t *buf, uint8_t *deltas_buf, uint16_t deltas_size, uint64_t line_point_size_bits, uint64_t stub_size_bits, uint8_t table_no )
			: is_compressed(true), line_point_size( (line_point_size_bits+7)/8 ), line_point_size_bits(line_point_size_bits)
			, stubs_size( Util::ByteAlign(stub_size_bits*(kEntriesPerPark-1))/8)
			, stub_size_bits(stub_size_bits), first_line_point( Util::SliceInt128FromBytes( buf + 3, 0, line_point_size_bits ) )
			, deltas_size(deltas_size), buf(buf), cur_stub_buf( buf + line_point_size + 3 /*for previous deltas pos*/ )

	{
		if( deltas_buf != nullptr && deltas_size > 0 )
			deltas = Encoding::ANSDecodeDeltas( deltas_buf , deltas_size, kEntriesPerPark - 1, kRValues[table_no] );
	}

	uint8_t * deltas_buf() const {
		assert( !is_compressed ); // this works for uncompressed only
		return buf + line_point_size + stubs_size + 2;
	}

	inline uint32_t DeltasSize() const { return deltas.size(); }
	inline uint32_t GetNextIdx() const { return next_idx; }
	inline uint32_t HasNext() const { return next_idx <= deltas.size(); }
	inline uint32_t GetSameDeltasCount() const { return same_deltas_count; }

	uint128_t NextLinePoint( uint32_t skip = 0 ){
		if( next_idx == 0 ){
			next_idx = 1;
			if( skip == 0 )
				return first_line_point;
			skip--; // it should be positive here
		}

		for( uint32_t i = 0; i <= skip; i++ )
			if( !MoveNext() ) return 0;

		return first_line_point + (((uint128_t)deltas_sum)<<stub_size_bits) + stubs_sum;
	}

	~ParkReader(){ if( is_compressed ) delete[]buf;	}

private:
	uint64_t stubs_sum = 0, deltas_sum = 0, last_delta = 0;
	uint8_t *cur_stub_buf;
	uint8_t stubs_start_bit = 0;
	uint32_t next_idx = 0, same_deltas_count = 0;

	inline bool MoveNext() {
		if( !HasNext() ) return false;

		if( stub_size_bits > 0 ){
			uint64_t stub = Util::EightBytesToInt( cur_stub_buf )<<stubs_start_bit;
			stub >>= 64 - stub_size_bits;
			stubs_start_bit += stub_size_bits;
			cur_stub_buf += stubs_start_bit/8;
			stubs_start_bit %= 8;

			stubs_sum += stub;
		}
		deltas_sum += deltas[next_idx-1];

		uint64_t new_delta = (deltas_sum<<stub_size_bits) + stubs_sum;
		same_deltas_count = new_delta == last_delta ? (same_deltas_count + 1) : 0;

		next_idx++;
		return true;
	}
};



struct LinePointInfo{
public:
	uint128_t orig_line_point = 0, full_line_point = 0;
	uint64_t position = (uint64_t)-1; // default unknown
	uint8_t skip_points = 0;
	LinePointInfo(){};
	LinePointInfo( uint64_t position, uint128_t orig_line_point, uint128_t full_line_point = 0 )
			: orig_line_point(orig_line_point), full_line_point(full_line_point), position(position){};

	void Set( uint8_t position, uint128_t orig_line_point, uint8_t skip_points ){
		this->position = position;
		this->orig_line_point = orig_line_point;
		this->skip_points = skip_points;
	}
};


struct LinePointTable2{
	LinePointInfo left;
	LinePointInfo right[kEntriesPerPark];
	uint16_t right_count = 0;
};

// 1st get input and find all appropriate points.
// 2nd collect all point to ranges includes ys[first] and min:max x[second]
//		every point could create number of ranges.
//	sort by x[second].min
//	evaluate buckets and find matches.
} // namespace tcompress


#endif // SRC_CPP_CMP_TOOLS_HPP_
