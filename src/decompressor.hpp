#ifndef SRC_CPP_DECOMPRESSOR_HPP_
#define SRC_CPP_DECOMPRESSOR_HPP_

#include "compressor.hpp"
#include <atomic>
#include <cstdint>
#include <cstring>
#ifndef _WIN32
#include <unistd.h>
#endif
#include <string>
#include <iostream>
#include <fstream>
#include <future>

#include "util.hpp"
#include "pos_constants.hpp"
#include "entry_sizes.hpp"
#include "encoding.hpp"
#include "compressor.hpp"

namespace TCompress {

struct LinePointMatcher{
public:
	const uint8_t k_size;
	LinePointMatcher( uint8_t k_size, uint8_t * plot_id, uint128_t lp1 )
			: k_size(k_size), f1(k_size,plot_id), f_d2( k_size, 2 ), f_d3( k_size, 3 ), ys1(CalculateYs( lp1) )
	{	}

	void ResetLP1( uint128_t lp1 ) { ys1 = CalculateYs(lp1); }

	uint64_t CalculateYs( uint128_t lp ){
		auto x1x2 = Encoding::LinePointToSquare( lp );
		Bits xs1( x1x2.first, k_size ), xs2( x1x2.second, k_size );
		auto ys1 = f1.CalculateF(xs1), ys2 = f1.CalculateF(xs2);

		return (ys1.GetValue() < ys2.GetValue() ? f_d2.CalculateBucket( ys1, xs1, xs2 )
																						: f_d2.CalculateBucket( ys2, xs2, xs1 )).first.GetValue();
	}

	bool CheckMatch( uint128_t lp2 ){
		auto ys2 = CalculateYs( lp2 );

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
	F1Calculator f1;
	FxCalculator f_d2, f_d3;
	uint64_t ys1;
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

// this is a cache for restored line points
LinePointsCache LPCache;

class Decompressor{
public:
	Decompressor() {}
	Decompressor(Decompressor& other) noexcept{
		k_size = other.k_size;
		memcpy( plot_id, other.plot_id, 32 );
		stubs_size = other.stubs_size;
		bits_cut_no = other.bits_cut_no;
		memcpy( table_pointers, other.table_pointers, 11*8 );
		memcpy( avg_delta_sizes, other.avg_delta_sizes, 7*8 );
		memcpy( parks_counts, other.parks_counts, 7*4 );
	}

	static Decompressor* CheckForCompressed( const std::string& filename, uint16_t memo_size, uint8_t k, const uint8_t *plot_id ){
		ReadFileWrapper disk_file( filename );
		uint8_t buf[4];
		disk_file.Read( 54, buf, 4 );
		if( memcmp( buf, tFormatDescription.c_str(), 4 ) == 0 ){
			Decompressor *res = new Decompressor();
			res->init( disk_file, memo_size, k, plot_id );
			return res;
		}
		return nullptr;
	}

	void init( std::ifstream& file, uint16_t memo_size, uint8_t k, const uint8_t *plot_id ){
		ReadFileWrapper disk_file( &file );
		init( disk_file, memo_size, k, plot_id );
	}

	void init( ReadFileWrapper& disk_file, uint16_t memo_size, uint8_t k, const uint8_t *plot_id ){
		k_size = k;
		memcpy( this->plot_id, plot_id, 32 );
		uint32_t line_point_size = EntrySizes::CalculateLinePointSize(k);
		stubs_size = EntrySizes::CalculateStubsSize(k);

		disk_file.Read( 60 + memo_size, (uint8_t*)table_pointers, 80 );
		disk_file.Read( &bits_cut_no, 1 );
		table2_cut = (0x80 & bits_cut_no) != 0;
		bits_cut_no = 0x7f & bits_cut_no;

		disk_file.Read( (uint8_t*)parks_counts, 28 );
		for( uint i = 0; i < 10; i++ )
			table_pointers[i] = bswap_64(table_pointers[i]);
		table_pointers[10] = disk_file.Size();

		for( uint i = 0; i < 7; i++ ){
			parks_counts[i] = bswap_32( parks_counts[i] );
			uint64_t table_size = table_pointers[ i == 6 ? 10 : (i+1)] - table_pointers[ i == 6 ? 9 : i];
			uint64_t static_size = GetStubsSize( i ); // i == 6 ? 0 : ( ( i == 0 ? first_table_stubs_size : stubs_size ) + line_point_size);
			avg_delta_sizes[i] = (table_size - (static_size+3)*parks_counts[i])/parks_counts[i];
		}
	}

	uint8_t GetNumberOfRemovedBits() { return bits_cut_no; }

	uint16_t ReadC3Park( std::ifstream& file, uint64_t idx, uint8_t *buf, uint16_t max_to_read ){
		if( idx >= parks_counts[6] )
			throw std::runtime_error( "too big C3 park index " + std::to_string(idx)
															 + " when max is " + std::to_string( parks_counts[6] ) );
		ReadFileWrapper disk_file( &file );

		uint8_t ibuf[6];
		disk_file.Read( table_pointers[9] + (idx-1)*3, ibuf, 6 );
		uint64_t pos;
		uint16_t delta_size;
		DeltasStorage::RestoreParkPositionAndSize( avg_delta_sizes[6], idx, ibuf, ibuf+3, pos, delta_size );

		if( delta_size > EntrySizes::CalculateC3Size(k_size) )
			throw std::runtime_error( "too big delta " + std::to_string( delta_size ) + ". possibly corrupted file ");
		if( delta_size <= max_to_read )
			disk_file.Read( table_pointers[9] + parks_counts[6]*3 + pos, buf, delta_size );

		return delta_size;
	}

	uint128_t ReadLinePoint( std::ifstream& file, // this need for parallel reading - will support it later
													uint8_t table_no /*0 is a first table*/,
													uint64_t position ){
		uint128_t line_point = LPCache.GetLinePoint( k_size, plot_id, table_no, position );
		if( line_point != 0 ) return line_point;

		switch( table_no ){
		case 0:
			line_point = ReadRealLinePoint( file, table_no, position );
			if( bits_cut_no > 0 ){
				std::cout << "tcompress: Missing line point cache. plot: ( k: " << (int)k_size
									<< ", id: " << Util::HexStr( plot_id, 32 ) << ", table: " << (int)table_no << ", position: "
									<< position << "), cache: (size: " << LPCache.GetSize() << ", count: "
									<< LPCache.GetCount() << ", important count: " << LPCache.GetImportantCount() << ")"
									<< std::endl;
				LPCache.IncreaseSize();
				line_point = RestoreLinePoint( line_point );
			}
			return line_point; // do not add here to cache sine it could be incorrect

		case 1:
			line_point = ReadRealLinePoint( file, table_no, position );
			if( table2_cut ){
				auto x1x2 = Encoding::LinePointToSquare( line_point << 11 );

				if( x1x2.second + kEntriesPerPark < x1x2.second )
					throw std::runtime_error( "unsupported case of moving first" );

				uint128_t valid_lp1 = LPCache.GetLinePoint( k_size, plot_id, 0, x1x2.first );
				std::vector<std::thread> threads;
				if( valid_lp1 == 0 )
					threads.emplace_back( [this, &valid_lp1](uint128_t lp ){ valid_lp1 = RestoreLinePoint(lp);}, ReadRealLinePoint( file, 0, x1x2.first )  );

				uint128_t valid_lp2 = 0;
				std::vector<uint128_t> restored_points, restored_points_first;
				{ // time to read second parks
					{
						uint32_t pos_in_park = x1x2.second % kEntriesPerPark;
						auto pReader = GetParkReader( file, 0, x1x2.second );
						restored_points.push_back( pReader.NextLinePoint( pos_in_park ) );
						while( pReader.HasNext() ) restored_points.push_back( pReader.NextLinePoint() );
					}
					if( restored_points.size() < kEntriesPerPark ){
						auto pReader = GetParkReader( file, 0, x1x2.second + kEntriesPerPark );
						while( restored_points.size() < kEntriesPerPark && pReader.HasNext() )
							restored_points.push_back( pReader.NextLinePoint() );
					}
				}

				for (auto& t : threads)  t.join(); // wait for valid_lp1
				threads.clear();

				//run some thread to read points ahead
				const uint32_t threads_number = 3; // equal to threads cound
				std::atomic_uint32_t evaluated_points = 0;
				uint8_t points_flags[restored_points.size()];
				memset( points_flags, 0, restored_points.size() );
				auto thread_func = [this, &restored_points, &evaluated_points, &valid_lp2, &points_flags](){
						for( uint32_t i = evaluated_points.fetch_add(1, std::memory_order_relaxed); valid_lp2 == 0 && i < restored_points.size();
								 i = evaluated_points.fetch_add(1, std::memory_order_relaxed) ){
							restored_points[i] = RestoreLinePoint( restored_points[i] );
							points_flags[i] = 1;
						}
					};
				for( uint32_t i = 0; i < threads_number; i++ )
					threads.emplace_back( thread_func );

				LinePointMatcher validator( k_size, plot_id, valid_lp1 );

				for( uint32_t i = 0; valid_lp2 == 0 && i < restored_points.size(); i++ ){
					while( points_flags[i] == 0 ) usleep( 1000 );
					// restored_points[i] = RestoreLinePoint( restored_points[i] );
					if( validator.CheckMatch( restored_points[i] ) ){
						valid_lp2 = restored_points[i];
						line_point = Encoding::SquareToLinePoint( x1x2.first, x1x2.second += i );
						for (auto& t : threads)  t.join(); // wait for threads to not fail...
						//std::cout << "evaluated overhead " << (evaluated_points-i) << std::endl;
					}
				}

				if( valid_lp2 == 0 ){ // match didn't succeed
					restored_points_first.push_back( valid_lp1 );
					while( valid_lp2 == 0 && (valid_lp1 = FindNextLinePoint( valid_lp1+1, bits_cut_no ) ) != 0 ){
						validator.ResetLP1( valid_lp1 );
						for( uint32_t i = 0; valid_lp2 == 0 && i < restored_points.size(); i++ ){
							if( validator.CheckMatch( restored_points[i] ) ) {
								valid_lp2 = restored_points[i];
								line_point = Encoding::SquareToLinePoint( x1x2.first, x1x2.second += i );
							}
						}
						if( valid_lp1 != 0 && valid_lp2 == 0 )
							restored_points_first.push_back( valid_lp1 );
					}
				}
				if( valid_lp2 == 0 ){ // match didn't succeed... now we need pass all restored points and find additional for each of them
					for( uint32_t i = 0; i < restored_points.size(); i++ ){
						uint128_t itmp = restored_points[i];
						while( valid_lp2 == 0 && itmp != 0 ){
							itmp = FindNextLinePoint( itmp + 1, bits_cut_no );
							if( itmp != 0 ) {
								for( uint32_t i = 0; valid_lp2 != 0 && i < restored_points_first.size(); i++ ){
									validator.ResetLP1( restored_points_first[i] );
									if( validator.CheckMatch( itmp ) ){
										valid_lp2 = itmp;
										line_point = Encoding::SquareToLinePoint( x1x2.first, x1x2.second + i );
									}
								}
							}
						}
					}
				}
				if( valid_lp2 == 0 )
					throw std::runtime_error( "cannot restore table 2 line point at position " + std::to_string(position) ); // TODO write all in exception include k and id

				LPCache.AddLinePoint( k_size, plot_id, 0, x1x2.first, valid_lp1, restored_points_first.size() > 0 );
				LPCache.AddLinePoint( k_size, plot_id, 0, x1x2.second, valid_lp2, restored_points_first.size() > 0 );

			}	else if( bits_cut_no > 0 )
			{ // Now we need to fine line_point from table 1 to cache them
				auto x1x2 = Encoding::LinePointToSquare( line_point );
				uint128_t valid_lp1 = LPCache.GetLinePoint( k_size, plot_id, 0, x1x2.first ),
						valid_lp2 = LPCache.GetLinePoint( k_size, plot_id, 0, x1x2.second );

				if( valid_lp1 == 0 || valid_lp2 == 0 ) { // one can survive in cache if it is important
					std::vector<uint128_t> lps1, lps2;
					std::vector<std::thread> threads;
					if( valid_lp1 == 0 )
						threads.emplace_back( [this, &valid_lp1](uint128_t lp ){ valid_lp1 = RestoreLinePoint(lp);}, ReadRealLinePoint( file, 0, x1x2.first )  );
					if( valid_lp2 == 0 )
						threads.emplace_back( [this, &valid_lp2](uint128_t lp ){ valid_lp2 = RestoreLinePoint(lp);}, ReadRealLinePoint( file, 0, x1x2.second )  );
					for (auto& t : threads)  t.join();

					bool match = CheckMatch( valid_lp1, valid_lp2 );
					while( !match && valid_lp1 ){
						lps1.push_back( valid_lp1 );
						valid_lp1 = FindNextLinePoint( valid_lp1 + 1, bits_cut_no );
						match = valid_lp1 != 0 && CheckMatch( valid_lp1, valid_lp2 );
					}
					while( !match && valid_lp2 != 0 ){
						lps2.push_back( valid_lp2 );
						valid_lp2 = FindNextLinePoint( valid_lp2 + 1, bits_cut_no );
						if( valid_lp2 != 0)
							for( uint i = 0; !match && i < lps1.size(); i++ )
								match = CheckMatch( valid_lp1 = lps1[i], valid_lp2 );
					}

					if( match ){
						LPCache.AddLinePoint( k_size, plot_id, 0, x1x2.first, valid_lp1, lps1.size() > 0 );
						LPCache.AddLinePoint( k_size, plot_id, 0, x1x2.second, valid_lp2, lps2.size() > 0 );
					} else {
						std::cout << "tcompress: No match at table2 position " << position << std::endl;
					}
				}
			}
			return LPCache.AddLinePoint( k_size, plot_id, table_no, position, line_point );

		default:
			return LPCache.AddLinePoint( k_size, plot_id, table_no, position, ReadRealLinePoint( file, table_no, position ) );
		}
	}

private:
	uint8_t k_size, plot_id[32];
	uint32_t stubs_size;
	uint8_t bits_cut_no;
	bool table2_cut = false;
	uint64_t table_pointers[11], avg_delta_sizes[7];
	uint32_t parks_counts[7];

	// match 2 lp from table 1 and supposed they are correct.
	bool CheckMatch( uint128_t lp1, uint128_t lp2 ){
		F1Calculator f1( k_size, plot_id );
		FxCalculator f_d2( k_size, 2 );
		FxCalculator f_d3( k_size, 3 );

		return CheckMatch( lp1, lp2, f1, f_d2, f_d3 );
	}
	// match 2 lp from table 1 and supposed they are correct.
	bool CheckMatch( uint128_t lp1, uint128_t lp2, F1Calculator &f1, FxCalculator &f_d2, FxCalculator &f_d3 ){

		auto xs1 = Encoding::LinePointToSquare( lp1 ), xs2 = Encoding::LinePointToSquare( lp2 );
		std::vector<Bits> xs;
		xs.emplace_back( xs1.first, k_size );
		xs.emplace_back( xs1.second, k_size );
		xs.emplace_back( xs2.first, k_size );
		xs.emplace_back( xs2.second, k_size );
		std::vector<Bits> ys;
		for( uint32_t i = 0; i < 4; i ++ )
			ys.push_back( f1.CalculateF(xs[i]) );

		// change order to proof order ?
		if( ys[0].GetValue() > ys[1].GetValue() ) {
			vectorswap( xs, 0, 1 );
			vectorswap( ys, 0, 1 );
		}
		if( ys[2].GetValue() > ys[3].GetValue() ) {
			vectorswap( xs, 2, 3 );
			vectorswap( ys, 2, 3 );
		}


		auto new_ys1 = f_d2.CalculateBucket( ys[0], xs[0], xs[1] ).first.GetValue();
		auto new_ys2 = f_d2.CalculateBucket( ys[2], xs[2], xs[3] ).first.GetValue();

		uint64_t cdiff = new_ys1/kBC - new_ys2/kBC;

		if( cdiff != 1 && cdiff != (uint64_t)-1 ) return false;

		std::vector<PlotEntry> bucket_L(1), bucket_R(1);
		bucket_L[0].y = std::min( new_ys1, new_ys2 );
		bucket_R[0].y = std::max( new_ys1, new_ys2 );
		if( f_d3.FindMatches( bucket_L, bucket_R, nullptr, nullptr ) != 1)
			return false;

		return true;
	}

	void vectorswap( std::vector<Bits> &v, uint32_t i, uint32_t j ){
		auto val = v[i]; 		v[i] = v[j]; 		v[j] = val;
	}

	inline uint32_t GetStubsSize( uint8_t table_no ){
		if( table_no == 6 ) return 0;
		if( table_no > (table2_cut?1:0) ) return stubs_size;
		if( table_no == 1 ) return ((k_size - kStubMinusBits - 11)*(kEntriesPerPark-1)+7)/8;
		return ((k_size - kStubMinusBits - bits_cut_no)*(kEntriesPerPark-1)+7)/8;
	}

	ParkReader GetParkReader( std::ifstream& file, // this need for parallel reading - will support it later
													 uint8_t table_no /*0 is a first table*/,
													 uint64_t position, bool for_single_point = false ){
		if( table_no >=6 )
			throw std::invalid_argument( "table couldn't be bigger than 5 for reading line point: " + std::to_string(table_no) );

		const uint64_t park_idx = position/kEntriesPerPark;
		assert( park_idx < parks_counts[table_no] );

		const uint8_t cur_bits_cut = (table_no?0:bits_cut_no) + ((table_no==1&&table2_cut)?11:0);
		const uint32_t cur_line_point_size_bits = k_size*2 - cur_bits_cut;
		const uint32_t cur_stub_size_bits = k_size - kStubMinusBits - cur_bits_cut;
		const uint32_t cur_line_point_size = (cur_line_point_size_bits+7)/8;
		const uint32_t stubs_size = (cur_stub_size_bits*(kEntriesPerPark-1)+7)/8;
		const uint64_t plps_size = 3 + cur_line_point_size + stubs_size;


		ReadFileWrapper disk_file( &file );

		if( for_single_point ){
			// Simplest case read first line point at the begining of park
			uint8_t line_point_buf[cur_line_point_size + 7];
			disk_file.Read( table_pointers[table_no] + plps_size*park_idx, line_point_buf, cur_line_point_size );
			return ParkReader( line_point_buf, cur_line_point_size_bits );
		}

		uint8_t * stubs_buf = new uint8_t[plps_size + 3 /*because we read 2 deltas pointers*/];
		disk_file.Read( table_pointers[table_no] + plps_size*park_idx - 3 /* to read prev delta pointer */, stubs_buf, 3+plps_size );
		uint64_t deltas_pos;
		uint16_t deltas_size;
		DeltasStorage::RestoreParkPositionAndSize( avg_delta_sizes[table_no], park_idx, stubs_buf, stubs_buf+plps_size, deltas_pos, deltas_size );

		if( deltas_size > EntrySizes::CalculateMaxDeltasSize( k_size, table_no + 1 ) )
			throw std::runtime_error( "incorrect deltas size " + std::to_string( deltas_size ) );

		uint8_t deltas_buf[deltas_size];
		disk_file.Read( table_pointers[table_no] + plps_size*parks_counts[table_no] + deltas_pos, deltas_buf, deltas_size );

		return ParkReader( stubs_buf, deltas_buf, deltas_size, cur_line_point_size_bits, cur_stub_size_bits, table_no );
	}

	uint128_t ReadRealLinePoint( std::ifstream& file, // this need for parallel reading - will support it later
												 uint8_t table_no /*0 is a first table*/,
												 uint64_t position ){
		const uint64_t pos_in_park = position%kEntriesPerPark;
		return GetParkReader( file, table_no, position, pos_in_park == 0 ).NextLinePoint( pos_in_park );

		// if( table_no >=6 )
		// 	throw std::invalid_argument( "table couldn't be bigger than 5 for reading line point: " + std::to_string(table_no) );

		// const uint32_t cur_line_point_size_bits = k_size*2 - (table_no?0:bits_cut_no);
		// const uint32_t cur_stub_size_bits = k_size - kStubMinusBits - (table_no?0:bits_cut_no);
		// const uint32_t cur_line_point_size = (cur_line_point_size_bits+7)/8;
		// const uint64_t plps_size = 3 + cur_line_point_size + (table_no ? stubs_size : first_table_stubs_size );
		// const uint64_t park_idx = position/kEntriesPerPark;

		// assert( park_idx < parks_counts[table_no] );

		// ReadFileWrapper disk_file( &file );

		// if( pos_in_park == 0 ){
		// 	// Simplest case read first line point at the begining of park
		// 	uint8_t line_point_buf[cur_line_point_size + 7];
		// 	disk_file.Read( table_pointers[table_no] + plps_size*park_idx, line_point_buf, cur_line_point_size );
		// 	return table_no == 0 ? Util::SliceInt128FromBytes( line_point_buf, 0, cur_line_point_size_bits )
		// 											: Util::SliceInt128FromBytes( line_point_buf, 0, cur_line_point_size_bits );
		// }

		// uint8_t stubs_buf[plps_size + 3 /*because we read 2 deltas pointers*/];
		// disk_file.Read( table_pointers[table_no] + plps_size*park_idx - 3 /* to read prev delta pointer */, stubs_buf, 3+plps_size );
		// uint64_t deltas_pos;
		// uint16_t deltas_size;
		// DeltasStorage::RestoreParkPositionAndSize( avg_delta_sizes[table_no], park_idx, stubs_buf, stubs_buf+plps_size, deltas_pos, deltas_size );

		// auto line_point = Util::SliceInt128FromBytes( stubs_buf+3, 0, cur_line_point_size_bits ); // base line point

		// std::vector<uint8_t> deltas;
		// {
		// 	if( deltas_size > EntrySizes::CalculateMaxDeltasSize( k_size, table_no + 1 ) )
		// 		throw std::runtime_error( "incorrect deltas size " + std::to_string( deltas_size ) );

		// 	uint8_t deltas_buf[deltas_size];
		// 	disk_file.Read( table_pointers[table_no] + plps_size*parks_counts[table_no] + deltas_pos, deltas_buf, deltas_size );
		// 	const double R = kRValues[table_no];
		// 	deltas = Encoding::ANSDecodeDeltas(deltas_buf, deltas_size, kEntriesPerPark - 1, R);
		// }

		// uint64_t sum_deltas = 0, sum_stubs = 0;

		// for( uint32_t i = 0, start_bit = (3+cur_line_point_size)*8;
		// 		 i < std::min((uint32_t)(position % kEntriesPerPark), (uint32_t)deltas.size());
		// 		 i++) {
		// 	uint64_t stub = cur_stub_size_bits == 0 ? 0 : Util::EightBytesToInt(stubs_buf + start_bit / 8); // seems max k is 56
		// 	stub <<= start_bit % 8;
		// 	stub >>= 64 - cur_stub_size_bits;

		// 	sum_stubs += stub;
		// 	start_bit += cur_stub_size_bits;

		// 	sum_deltas += deltas[i];
		// }

		// uint128_t big_delta = ((uint128_t)sum_deltas << cur_stub_size_bits) + sum_stubs;

		// auto pR = GetParkReader( file, table_no, position, pos_in_park == 0 );
		// auto byParkReader = pR.NextLinePoint( pos_in_park );

		// assert( byParkReader == line_point + big_delta );

		// return line_point + big_delta;
	}


	uint128_t RestoreLinePoint( uint128_t cutted_line_point ){

		auto lp = FindNextLinePoint( cutted_line_point<<bits_cut_no, bits_cut_no );
		if( lp != 0 ) return lp;

		throw std::runtime_error( "Cannot restore linepoint " + std::to_string((uint64_t)(cutted_line_point>>64) )
														 + ":" + std::to_string( (uint64_t)cutted_line_point ));
	}

	uint128_t FindNextLinePoint( uint128_t line_point, uint8_t removed_bits_no ) const {
		assert( removed_bits_no >= kBatchSizes );

		auto x1x2 = Encoding::LinePointToSquare( line_point );

		F1Calculator f1( k_size, plot_id );
		FxCalculator f( k_size, 2 );

		uint64_t firstY = f1.CalculateF( Bits( x1x2.first, k_size ) ).GetValue();
		const uint64_t BATCH_SIZE = 1U << kBatchSizes; //256
		uint64_t batch[BATCH_SIZE];

		std::vector<PlotEntry> bucket_L(1);
		std::vector<PlotEntry> bucket_R(1);

		uint64_t b = x1x2.second, left_to_do = (1U << removed_bits_no) - (line_point&((1ULL << removed_bits_no)-1));
		while( left_to_do > 0 ){
			f1.CalculateBuckets( b, BATCH_SIZE, batch );
			uint64_t cur_batch_size = std::min( left_to_do, std::min( BATCH_SIZE, x1x2.first - b ) ); // possible to minimize with left_to_do
			for( uint32_t i = 0; i < cur_batch_size; i++ ){
				uint64_t cdiff = firstY / kBC - batch[i] / kBC;
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
				b = 0;
			}
		}

		return 0; // NOT FOUND
	}

};

} // endof namespace

#endif  // SRC_CPP_COMPRESSOR_HPP_
