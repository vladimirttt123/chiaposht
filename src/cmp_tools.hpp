#ifndef SRC_CPP_CMP_TOOLS_HPP_
#define SRC_CPP_CMP_TOOLS_HPP_

#include "calculate_bucket.hpp"
#include "encoding.hpp"
#include <atomic>
#include <cstdint>
#include <cstring>
#include <memory>
#include <thread>
#ifndef _WIN32
#include <unistd.h>
#endif
#include <string>
#include <iostream>
#include <fstream>

#include "util.hpp"
#include "pos_constants.hpp"



namespace TCompress {

#ifdef THREADS_PER_LINE_POINT
const uint32_t THREADS_PER_LP = THREADS_PER_LINE_POINT;
#else
const uint32_t THREADS_PER_LP = 4;
#endif

// This is replacement made for kFormatDescription in compress plots
const std::string tFormatDescription = "t0.1";
const char* plotMagicFrase = "Proof of Space Plot";

const char * program_header = "*** Chia plot compressing software made by Vladimir T\n"
			"*** If this compression is helpfull for you, please consider donate\n"
			"***   xch1ch6s3q0enuj9wtemn473gkkvj0u8vlggypr375mk547e7aa48hmsql74e8\n";


const uint32_t overdraftPointerSize = 3;

struct ThreadDeleter{
	void operator()(std::thread * p) const
	{
		p->join();
		delete p;
	};
};

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

	void TotalEndToBuf( uint64_t idx, uint8_t *buf, bool new_align = true ){
		assert( overdraftPointerSize == 3 );
		if( new_align ){
			buf[0] = total_size >> ((idx&1)?32:16);
			buf[1] = total_size >> ((idx&1)?24:8);
			buf[2] = total_size;
		}else {
			buf[0] = total_size >> ((idx&1)?24:16);
			buf[1] = total_size>>8;
			buf[2] = total_size;
		}
	}

	// check all partially saved position could be restored
	bool IsDeltasPositionRestorable( bool is_new_align = true ){
		uint64_t initial_total_size = total_size;
		uint64_t park_avg_size = total_size/parks_count;
		uint8_t buf[overdraftPointerSize*2];
		uint8_t *buf_prev = buf, *buf_cur = buf + overdraftPointerSize;
		uint64_t delta_pos;
		uint16_t delta_size;


		total_size = 0; // to work with function TotalEndToBuf

		for( uint64_t i = 0; i < parks_count; i++ ){
			total_size += all_sizes[i];
			TotalEndToBuf( i, buf_cur );
			RestoreParkPositionAndSize( is_new_align, park_avg_size, i, buf_prev, buf_cur, delta_pos, delta_size );
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

	void showStats(){
		uint32_t min_park_size = (uint32_t)-1, max_park_size = 0, number_of_zeros = 0;
		for( uint64_t i = 0; i < parks_count-1; i++ ){
			if( min_park_size > all_sizes[i] ) min_park_size = all_sizes[i];
			if( max_park_size < all_sizes[i] ) max_park_size = all_sizes[i];
			if( all_sizes[i] == 0 ) number_of_zeros++;
		}
		std::cout<< " min_park_size: " << min_park_size << "; max_park_size: " << max_park_size
							<< ", number of zeros: " << number_of_zeros << "(" << ((number_of_zeros/(double) parks_count)*100) << "%)" << std::endl;
	}

	static void RestoreParkPositionAndSize( bool is_new_align, uint64_t park_avg_size, uint64_t park_idx, uint8_t *prev_buf, uint8_t *cur_buf,
																				 uint64_t &delta_position, uint16_t &encoded_delta_size ){
		if( is_new_align ){
			uint64_t cur_pos = cur_buf[2];

			if( park_idx == 0 ){
				delta_position = 0;
				encoded_delta_size = cur_pos;
				return;
			}

			uint64_t prev_pos = prev_buf[2];
			encoded_delta_size = (uint8_t)(cur_pos - prev_pos);

			// now we restore 5 bottom bytes of position it is nought for more than 1Tb pointers
			if( park_idx&1 ){ // top bytes after known and low before we need full position before
				prev_pos |= (((uint64_t)prev_buf[0])<<16) | (((uint64_t)prev_buf[1])<<8);
				delta_position = ((((uint64_t)cur_buf[0])<<32) | (((uint64_t)cur_buf[1])<<24) | ((prev_pos+encoded_delta_size)&0xffffff)) - encoded_delta_size;
			}else{
				cur_pos |= (((uint64_t)cur_buf[0])<<16) | (((uint64_t)cur_buf[1])<<8);
				delta_position = (((uint64_t)prev_buf[0])<<32) | (((uint64_t)prev_buf[1])<<24) | ((cur_pos - encoded_delta_size)&0xffffff);
			}
		}else { // old align
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

	std::ifstream * getStream() const { return disk_file; }
private:
	std::ifstream *disk_file;
	bool isExternal;
};


struct ParkReader{
public:
	const bool is_compressed;
	const uint64_t line_point_size, line_point_size_bits, stubs_size, stub_size_bits, overdraft_pos = 0;
	const uint128_t first_line_point;
	uint16_t overdraft_size = 0, deltas_size;
	uint8_t *buf;
	std::vector<uint8_t> deltas;

	// used to create park reader with single point
	ParkReader( uint8_t *buf, uint64_t line_point_size_bits )
			: is_compressed(false), line_point_size(0), line_point_size_bits(line_point_size_bits)
			, stubs_size(0),  stub_size_bits(0)
			, first_line_point( Util::SliceInt128FromBytes( buf, 0, line_point_size_bits ) )
	{	}

	// used for uncompressed plots
	ParkReader( uint8_t *buf, uint64_t line_point_size_bits, uint64_t stub_size_bits, uint32_t max_deltas_size, uint8_t table_no )
			: is_compressed(false), line_point_size( (line_point_size_bits+7)/8 ), line_point_size_bits(line_point_size_bits)
			, stubs_size( Util::ByteAlign(stub_size_bits*(kEntriesPerPark-1))/8), stub_size_bits(stub_size_bits)
			, first_line_point( Util::SliceInt128FromBytes( buf, 0, line_point_size_bits ) )
			, deltas_size( ((uint32_t)buf[line_point_size+stubs_size]) | (((uint32_t)buf[line_point_size+stubs_size+ 1])<<8) )
			, buf(buf), cur_stub_buf( buf + line_point_size ), src_stubs_buf(cur_stub_buf), src_check_point_buf(buf)
	{
		uint8_t *deltas_buf = buf + line_point_size + stubs_size + 2/* uint16 of deltas sizes */;

		if( deltas_size&0x8000 ){ // this is uncompressed deltas
			deltas_size &= 0x7fff; // this is real size of uncompressed deltas
			if( deltas_size > 0 ) {
				if( deltas_size > 2 )
					throw std::runtime_error( "Unsupported case of uncompressed deltas" );

				// this could be only on last park
				// the ANS does not encode such lengthes than we add artifitially one or 2 entries just to compress them
				// this additional entries would never be mentioned from higer table and than shouldn't do any harm.
				deltas.push_back( deltas_buf[0] );
				deltas.push_back( deltas_buf[1] );
				deltas.push_back( 0 );
				uint32_t new_size = Encoding::ANSEncodeDeltas( deltas, kRValues[table_no-1], deltas_buf );
				if( new_size == 0 || new_size > max_deltas_size )
					throw std::runtime_error( "Unsupported case of uncompressed deltas with short length" );
				deltas_size = new_size;
			}
		} else {
			if( deltas_size > max_deltas_size )
				throw std::runtime_error( "incorrect deltas size " + std::to_string( deltas_size ) );

			if( deltas_size > 0 )
				deltas = Encoding::ANSDecodeDeltas( deltas_buf , deltas_size, kEntriesPerPark - 1, kRValues[table_no-1] );
		}
	}

	// used for compressed plots
	// WARNING variable buf is stored and deleted in destructor
	ParkReader( uint8_t *buf, uint8_t*check_point_buf, uint8_t*stubs_buf, uint8_t *deltas_buf, uint16_t deltas_size,
						 uint64_t line_point_size_bits, uint64_t stub_size_bits, uint8_t table_no,
						 uint64_t overdraft_pos = 0, uint16_t overdraf_size = 0 )
			: is_compressed(true), line_point_size( (line_point_size_bits+7)/8 ), line_point_size_bits(line_point_size_bits)
			, stubs_size( Util::ByteAlign(stub_size_bits*(kEntriesPerPark-1))/8), stub_size_bits(stub_size_bits)
			, overdraft_pos(overdraft_pos), first_line_point( Util::SliceInt128FromBytes( check_point_buf, 0, line_point_size_bits ) )
			, overdraft_size(overdraf_size), deltas_size(deltas_size), buf(buf), cur_stub_buf( stubs_buf )
			, src_deltas_buf( deltas_buf ), src_stubs_buf(cur_stub_buf), src_check_point_buf( check_point_buf )

	{
		if( deltas_buf != nullptr && deltas_size > 0 )
			deltas = Encoding::ANSDecodeDeltas( deltas_buf , deltas_size, kEntriesPerPark - 1, kRValues[table_no] );
	}

	inline uint8_t * deltas_buf() const { return is_compressed ? src_deltas_buf : (buf + line_point_size + stubs_size + 2);}
	inline const uint8_t* stubs_buf() const { return src_stubs_buf; }
	inline uint8_t * stubs_overdraft_buf() {return src_stubs_buf + stubs_size - overdraft_size; }

	const uint8_t* check_point_buf() const { return src_check_point_buf; }

	inline uint16_t StubsLinePointsCount() const { return (stubs_size-overdraft_size)*8/stub_size_bits; }
	inline uint32_t DeltasSize() const { return deltas.size(); }
	inline uint32_t GetNextIdx() const { return next_idx; }
	inline bool HasNext() const { return next_idx <= deltas.size(); }
	inline bool HasNextInStub() const { return HasNext() && next_idx <= StubsLinePointsCount(); }
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

	~ParkReader(){ if( is_compressed ) delete[]buf; }

private:
	uint64_t stubs_sum = 0, deltas_sum = 0, last_delta = 0;
	uint8_t *cur_stub_buf, *src_deltas_buf = nullptr, *src_stubs_buf, *src_check_point_buf;
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




uint128_t FindNextLinePoint( uint128_t line_point, uint8_t removed_bits_no, uint8_t k_size, F1Calculator &f1, FxCalculator &f ) {
	assert( removed_bits_no >= kBatchSizes );

	auto x1x2 = Encoding::LinePointToSquare( line_point );
	uint64_t firstY = f1.CalculateF( Bits( x1x2.first, k_size ) ).GetValue();
	uint64_t firstYkBC = firstY / kBC;

	const uint64_t BATCH_SIZE = 1U << kBatchSizes; //256
	uint64_t batch[BATCH_SIZE];

	std::vector<PlotEntry> bucket_L(1);
	std::vector<PlotEntry> bucket_R(1);

	uint64_t b = x1x2.second, left_to_do = (1ULL << removed_bits_no) - (line_point&((1ULL << removed_bits_no)-1));
	while( left_to_do > 0 ){
		f1.CalculateBuckets( b, BATCH_SIZE, batch );
		uint64_t cur_batch_size = std::min( left_to_do, std::min( BATCH_SIZE, x1x2.first - b ) ); // possible to minimize with left_to_do
		for( uint32_t i = 0; i < cur_batch_size; i++ ){
			uint64_t cdiff = firstYkBC - batch[i] / kBC;
			if( cdiff == 1 ){
				bucket_L[0].y = batch[i];
				bucket_R[0].y = firstY;
			}else if( cdiff == (uint64_t)-1 ){
				bucket_L[0].y = firstY;
				bucket_R[0].y = batch[i];
			}else continue;

			if( f.FindMatches( bucket_L, bucket_R, nullptr, nullptr ) == 1)
				return Encoding::SquareToLinePoint( x1x2.first, b+i );
		}
		left_to_do -= cur_batch_size;
		b += cur_batch_size;
		if( b == x1x2.first ){
			firstY = f1.CalculateF( Bits( ++(x1x2.first), k_size ) ).GetValue();
			firstYkBC = firstY / kBC;
			b = 0;
		}
	}

	return 0; // NOT FOUND
}


inline uint128_t FindNextLinePoint( uint128_t line_point, uint8_t removed_bits_no, uint8_t k_size, const uint8_t * plot_id ) {
	assert( removed_bits_no >= kBatchSizes );

	F1Calculator f1( k_size, plot_id );
	FxCalculator f( k_size, 2 );
	return FindNextLinePoint( line_point, removed_bits_no, k_size, f1, f );
}

// Checks matches at table 2 level
struct LinePointMatcher{
	const uint8_t k_size, *plot_id;
	F1Calculator f1;
	FxCalculator f_d2, f_d3;

	LinePointMatcher( uint8_t k_size, const uint8_t * plot_id, uint128_t lp1 = 0 )
			: k_size(k_size), plot_id(plot_id), f1(k_size,plot_id), f_d2( k_size, 2 ), f_d3( k_size, 3 ), ys1(CalculateYs( lp1) )
	{	}

	uint64_t ResetLP1( uint128_t lp1 ) { return ys1 = CalculateYs(lp1); }
	uint64_t ResetLP1( uint64_t x1, uint64_t x2 ) { return ys1 = CalculateYs(x1, x2); }


	uint64_t CalculateYs( uint128_t lp ){
		auto x1x2 = Encoding::LinePointToSquare( lp );
		return CalculateYs( x1x2.first, x1x2.second );
	}

	uint64_t CalculateYs( uint64_t x1, uint64_t x2 ){
		Bits xs1( x1, k_size ), xs2( x2, k_size );
		auto ys1 = f1.CalculateF(xs1), ys2 = f1.CalculateF(xs2);

		return (ys1.GetValue() < ys2.GetValue() ? f_d2.CalculateBucket( ys1, xs1, xs2 )
																						: f_d2.CalculateBucket( ys2, xs2, xs1 )).first.GetValue();

	}

	inline bool CheckMatch( uint128_t lp2 ){
		auto x1x2 = Encoding::LinePointToSquare( lp2 );
		return CheckMatch( x1x2.first, x1x2.second );
	}

	inline bool CheckMatch( uint64_t x1, uint64_t x2 ){
		return isYsMatch( ys1, CalculateYs( x1, x2 ) );
	}
	inline bool isYsMatch( uint64_t ys2 ){ return isYsMatch( ys1, ys2 );	}

	bool isYsMatch( uint64_t ys1, uint64_t ys2 ){
		uint64_t cdiff = ys1/kBC - ys2/kBC;

		if( cdiff != 1 && cdiff != (uint64_t)-1 ) return false;

		std::vector<PlotEntry> bucket_L(1), bucket_R(1);
		bucket_L[0].y = std::min( ys1, ys2 );
		bucket_R[0].y = std::max( ys1, ys2 );
		if( f_d3.FindMatches( bucket_L, bucket_R, nullptr, nullptr ) != 1)
			return false;

		return true;
	}

private:
	uint64_t ys1;
};

inline uint128_t FindNextLinePoint( uint128_t line_point, uint8_t removed_bits_no, LinePointMatcher &validator ) {
	assert( removed_bits_no >= kBatchSizes );
	return FindNextLinePoint( line_point, removed_bits_no, validator.k_size, validator.f1, validator.f_d2 );
}
inline uint128_t RestoreLinePoint( uint128_t cutted_line_point, uint8_t removed_bits_no, LinePointMatcher &validator ) {
	return FindNextLinePoint( cutted_line_point<<removed_bits_no, removed_bits_no, validator );
}

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


struct LinePointToMatch{
	uint16_t orig_idx = 0;
	uint64_t x1 = 0, x2 = 0, ys = 0;
	inline uint128_t LinePoint() const { return  Encoding::SquareToLinePoint( x1, x2 );}
};

struct MatchVector{
	const uint16_t max_size;
	std::unique_ptr<LinePointToMatch[]> points;
	std::atomic_uint_fast16_t size = 0;
	MatchVector( uint16_t max_size ) :max_size(max_size), points( new LinePointToMatch[max_size] ){	}

	void Add( uint16_t idx, uint64_t x1, uint64_t x2, uint64_t ys ){
		auto cur_idx = size.fetch_add(1, std::memory_order_relaxed);
		if( cur_idx < max_size )
			points[cur_idx] = { idx, x1, x2, ys };
		else
			std::cout << "WARNING!!! MatchVector overflow. skipping insert of value" << std::endl;
			// throw std::overflow_error( "The vector cannot hold more than " + std::to_string(max_size) + " values " );
	}

	void Add( uint16_t idx, uint128_t lp, uint64_t ys ){
		auto x1x2 = Encoding::LinePointToSquare( lp );
		Add( idx, x1x2.first, x1x2.second, ys );
	}
};

struct Table2MatchData{
	MatchVector left, right;
	uint128_t matched_left = 0, matched_right = 0;
	uint16_t matched_right_idx = 0;

	Table2MatchData() :left(5), right(kEntriesPerPark+100) {}

	// add line point to matching data
	// and returns if status of matched changed after this adding
	bool Add( int16_t idx, uint64_t x1, uint64_t x2, LinePointMatcher &validator ){
		if( matched_left != 0 ) return false; // may not so good for threads...
		if( idx < 0 ) AddLeft( x1, x2, validator );
		else AddRight( idx, x1, x2, validator );
		return matched_left != 0;
	}

	bool AddLeft( uint128_t lp, LinePointMatcher &validator ){
		if( lp == 0 ) return false;
		auto x1x2 = Encoding::LinePointToSquare(lp);
		AddLeft( x1x2.first, x1x2.second, validator );
		return true;
	}

	void AddLeft( uint64_t x1, uint64_t x2, LinePointMatcher &validator ){
		if( matched_left != 0 ) return; // match done than do not add to matched
		left.Add( -1, x1, x2, validator.ResetLP1( x1, x2 ) );
		// validate left to all exists right
		// some points could be validated twice but this is the price for not perform sync
		for( uint32_t i = 0; i < right.size && matched_left == 0; i++ )
			if( validator.isYsMatch( right.points[i].ys ) ){
				matched_left = Encoding::SquareToLinePoint( x1, x2 );
				matched_right = right.points[i].LinePoint();
				matched_right_idx = right.points[i].orig_idx;
			}
	}

	bool AddRight( uint16_t idx, uint128_t lp, LinePointMatcher &validator ){
		if( lp == 0 ) return false;
		auto x1x2 = Encoding::LinePointToSquare(lp);
		AddRight( idx, x1x2.first, x1x2.second, validator );
		return true;
	}
	void AddRight( uint16_t idx, uint64_t x1, uint64_t x2, LinePointMatcher &validator ){
		if( matched_left != 0 ) return;// do not add to mathed
		right.Add( idx, x1, x2, validator.ResetLP1( x1, x2 ) );
		for( uint32_t i = 0; i < left.size && matched_left == 0; i++ )
			if( validator.isYsMatch( left.points[i].ys ) ){
				matched_right = Encoding::SquareToLinePoint( x1, x2 );
				matched_left = left.points[i].LinePoint();
				matched_right_idx = idx;
			}
	}


	bool AddBulk( LinePointInfo &left_LP, uint16_t first_idx, std::vector<uint128_t> src_line_points,
								uint8_t removed_bits_no, uint64_t table2_pos, uint64_t table1_start_pos,
								LinePointMatcher &pvalidator, uint32_t threads_no = THREADS_PER_LP ){
		std::atomic_uint_fast16_t next_point_idx = 0;

		auto thread_func = [this, &next_point_idx, &table1_start_pos, &table2_pos, &pvalidator, &src_line_points, &removed_bits_no, &first_idx](){
			LinePointMatcher validator( pvalidator.k_size, pvalidator.plot_id );

			for( uint16_t i = next_point_idx.fetch_add(1, std::memory_order_relaxed);
					 matched_left == 0/*not found yet*/ && i < src_line_points.size();
					 i = next_point_idx.fetch_add(1, std::memory_order_relaxed) ){
				if( (i+1) < src_line_points.size() && src_line_points[i] == src_line_points[i+1] ) continue; // equally cut points evaluated in the same thread
				auto restored_lp = RestoreLinePoint( src_line_points[i], removed_bits_no, validator );
				//if( CheckRestored( restored_lp, x1x2.second + i, position, false ) ){
				if( restored_lp == 0 ){
					std::cout << "WARNING!!! Cannot restore line point at pos " << (table1_start_pos+i) << " from table2 pos " << table2_pos << " - SKIPPING" << std::endl;
				}else {
					AddRight( i + first_idx, restored_lp, validator );
					if( i > 0 && src_line_points[i-1] == src_line_points[i] )
						while(  matched_left == 0/*not found yet*/ && (restored_lp = FindNextLinePoint( restored_lp+1, removed_bits_no, validator) ) )
							AddRight( --i, restored_lp, validator );
				}
			};
		};

		// check if we need to restore left line point
		std::unique_ptr<std::thread,ThreadDeleter> lp1_thread;
		if( left.size == 0 && left_LP.orig_line_point != 0 ) {
			lp1_thread.reset( new std::thread( [ this, &table2_pos, &left_LP, &pvalidator, &removed_bits_no ](){
				uint128_t lp = RestoreLinePoint( left_LP.orig_line_point, removed_bits_no, pvalidator );
				for( uint i = 0; i < left_LP.skip_points; i++ )
					lp = FindNextLinePoint( lp + 1, removed_bits_no, pvalidator );
				if( lp == 0 ) {
					std::cout << "WARNING!!! Cannot restore left line point at pos " << left_LP.position << " from table2 pos "
										<< table2_pos << " - SKIPPING" << std::endl;
				}
				AddLeft( lp, pvalidator );
			} ) );
		}
		{
			std::unique_ptr<std::thread, ThreadDeleter> threads[THREADS_PER_LP];

			for( uint32_t i = 0; i < THREADS_PER_LP; i++ )
				threads[i].reset( new std::thread( thread_func ) );
		}

		return matched_left != 0;
	}

	bool RunSecondRound( uint8_t removed_bits_no, uint64_t table2_pos, uint64_t table1_start_pos,
											LinePointMatcher &validator, uint32_t threads_no = THREADS_PER_LP ){

		// check more left points
		for( uint128_t lp = FindNextLinePoint( left.points[0].LinePoint() + 1, removed_bits_no, validator );
				 lp != 0; lp = FindNextLinePoint( lp + 1, removed_bits_no, validator ) ){
			AddLeft( lp, validator );
			if( matched_left != 0 ) return true;// return if match found
		}
		// TODO make by threads
		// check more right points - it could be long process and seems need threads!!!
		for( uint32_t size = right.size, i = 0; i < size; i++ )
			for( uint128_t lp = FindNextLinePoint( left.points[i].LinePoint() + 1, removed_bits_no, validator );
					 lp != 0; lp = FindNextLinePoint( lp + 1, removed_bits_no, validator ) ){
				AddRight( right.points[i].orig_idx, lp, validator );
				if( matched_left != 0 ) return true;// return if match found
			}

		return matched_left != 0;
	}
};

} // namespace tcompress


#endif // SRC_CPP_CMP_TOOLS_HPP_
