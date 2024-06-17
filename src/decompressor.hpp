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

struct AtomicAdder{
	const uint32_t inited_at;
	std::atomic_uint32_t &src;
	AtomicAdder( std::atomic_uint32_t &src )
			:inited_at( src.fetch_add( 1, std::memory_order_relaxed ) ), src(src)
	{	}

	~AtomicAdder(){
		src.fetch_sub( 1, std::memory_order_relaxed );
	}
};

struct LinePointCacheEntry{
public:
	uint64_t partial_id = 0, position = 0;
	uint128_t line_point = 0;
};

struct LinePointsCache{
public:
	const uint32_t size;
	LinePointsCache( uint32_t size = 1024) : size(size), cache(new LinePointCacheEntry[size] ){

	}
	void AddLinePoint( const uint8_t * id, uint64_t position, uint128_t line_point ){
		if( line_point == 0 ) return;
		std::lock_guard<std::mutex> lk(mut);
		cache[next_pos].line_point = line_point;
		cache[next_pos].position = position;
		cache[next_pos++].partial_id = ((uint64_t*)id)[0];
		if( count < next_pos ) count = next_pos;
		next_pos %= size;
	}

	uint128_t GetLinePoint( const uint8_t * id, uint64_t position ){
		std::lock_guard<std::mutex> lk(mut);
		for( uint32_t i = 1; i <= count; i++ ){
			uint32_t idx = (next_pos - i)%size;
			if( cache[idx].position == position && cache->partial_id == ((uint64_t*)id)[0] )
				return cache[idx].line_point;
		}
		return 0; // NOT FOUND
	}

	~LinePointsCache(){ delete[] cache; }
private:
	LinePointCacheEntry *cache;
	uint32_t next_pos = 0, count = 0;
	std::mutex mut;
};

LinePointsCache LPCache;

class Decompressor{
public:
	Decompressor() {}
	Decompressor(Decompressor& other) noexcept{
		k_size = other.k_size;
		memcpy( plot_id, other.plot_id, 32 );
		stubs_size = other.stubs_size;
		first_table_stubs_size = other.first_table_stubs_size;
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
		disk_file.Read( (uint8_t*)parks_counts, 28 );
		for( uint i = 0; i < 10; i++ )
			table_pointers[i] = bswap_64(table_pointers[i]);
		table_pointers[10] = disk_file.Size();

		first_table_stubs_size = Util::ByteAlign((kEntriesPerPark - 1) * (k - kStubMinusBits - bits_cut_no)) / 8;

		for( uint i = 0; i < 7; i++ ){
			parks_counts[i] = bswap_32( parks_counts[i] );
			uint64_t table_size = table_pointers[ i == 6 ? 10 : (i+1)] - table_pointers[ i == 6 ? 9 : i];
			uint64_t static_size = i == 6 ? 0 : ( ( i == 0 ? first_table_stubs_size : stubs_size ) + line_point_size);
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
		uint128_t line_point;
		switch( table_no ){
		case 0:
			line_point = LPCache.GetLinePoint( plot_id, position );
			if( line_point == 0 ){
				std::cout << "Missing line point cache - qualit/proof could be incerroct at position " << position << std::endl;
				return ReadLinePointReal( file, table_no, position );
			}
			return line_point;

		case 1:
			line_point = ReadLinePointReal( file, table_no, position );
			// Now we need to fine line_point from table 1 to cache them
			{
				auto x1x2 = Encoding::LinePointToSquare( line_point );
				if( LPCache.GetLinePoint( plot_id, x1x2.first ) == 0 || LPCache.GetLinePoint( plot_id, x1x2.second ) == 0 ) {
					std::vector<uint128_t> lps1, lps2;
					uint128_t valid_lp1 = ReadLinePointReal( file, 0, x1x2.first ),
							valid_lp2 = ReadLinePointReal( file, 0, x1x2.second );
					bool match = CheckMatch( valid_lp1, valid_lp2 );
					while( !match && valid_lp1 ){
						lps1.push_back( valid_lp1 );
						valid_lp1 = RestoreNextLinePoint( valid_lp1 + 1, bits_cut_no );
						match = valid_lp1 != 0 && CheckMatch( valid_lp1, valid_lp2 );
					}
					while( !match && valid_lp2 != 0 ){
						lps2.push_back( valid_lp2 );
						valid_lp2 = RestoreNextLinePoint( valid_lp2 + 1, bits_cut_no );
						if( valid_lp2 != 0)
							for( uint i = 0; !match && i < lps1.size(); i++ )
								match = CheckMatch( valid_lp1 = lps1[i], valid_lp2 );
					}

					if( match ){
						LPCache.AddLinePoint( plot_id, x1x2.first, valid_lp1 );
						LPCache.AddLinePoint( plot_id, x1x2.second, valid_lp2 );
					} else {
						std::cout << "no match at table2 position " << position << std::endl;
					}
				}
			}
			return line_point;

		default:
			return ReadLinePointReal( file, table_no, position );
		}
	}

private:
	uint8_t k_size, plot_id[32];
	uint32_t stubs_size, first_table_stubs_size;
	uint8_t bits_cut_no;
	uint64_t table_pointers[11], avg_delta_sizes[7];
	uint32_t parks_counts[7];
	std::atomic_uint32_t threads_count = 0;

	// match 2 lp from table 1 and supposed they are correct.
	bool CheckMatch( uint128_t lp1, uint128_t lp2 ){
		F1Calculator f1( k_size, plot_id );
		FxCalculator f_d2( k_size, 2 );
		FxCalculator f_d3( k_size, 3 );


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
		// TODO check match...
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

	uint128_t ReadLinePointReal( std::ifstream& file, // this need for parallel reading - will support it later
												 uint8_t table_no /*0 is a first table*/,
												 uint64_t position ){

		if( table_no >=6 )
			throw std::invalid_argument( "table couldn't be bigger than 5 for reading line point: " + std::to_string(table_no) );

		const uint32_t cur_line_point_size_bits = k_size*2 - (table_no?0:bits_cut_no);
		const uint32_t cur_stub_size_bits = k_size - kStubMinusBits - (table_no?0:bits_cut_no);
		const uint32_t cur_line_point_size = (cur_line_point_size_bits+7)/8;
		const uint64_t plps_size = 3 + cur_line_point_size + (table_no ? stubs_size : first_table_stubs_size );
		const uint64_t park_idx = position/kEntriesPerPark;

		assert( park_idx < parks_counts[table_no] );

		ReadFileWrapper disk_file( &file );

		if( (position%kEntriesPerPark) == 0 ){
			// Simplest case read first line point at the begining of park
			uint8_t line_point_buf[cur_line_point_size + 7];
			disk_file.Read( table_pointers[table_no] + plps_size*park_idx, line_point_buf, cur_line_point_size );
			return (table_no == 0 && bits_cut_no > 0 ) ?
					RestoreLinePoint( Util::SliceInt128FromBytes( line_point_buf, 0, cur_line_point_size_bits ) )
																								: Util::SliceInt128FromBytes( line_point_buf, 0, cur_line_point_size_bits );
		}

		uint8_t stubs_buf[plps_size + 3 /*because we read 2 deltas pointers*/];
		disk_file.Read( table_pointers[table_no] + plps_size*park_idx - 3 /* to read prev delta pointer */, stubs_buf, 3+plps_size );
		uint64_t deltas_pos;
		uint16_t deltas_size;
		DeltasStorage::RestoreParkPositionAndSize( avg_delta_sizes[table_no], park_idx, stubs_buf, stubs_buf+plps_size, deltas_pos, deltas_size );

		auto line_point = Util::SliceInt128FromBytes( stubs_buf+3, 0, cur_line_point_size_bits ); // base line point

		std::vector<uint8_t> deltas;
		{
			if( deltas_size > EntrySizes::CalculateMaxDeltasSize( k_size, table_no + 1 ) )
				throw std::runtime_error( "incorrect deltas size " + std::to_string( deltas_size ) );

			uint8_t deltas_buf[deltas_size];
			disk_file.Read( table_pointers[table_no] + plps_size*parks_counts[table_no] + deltas_pos, deltas_buf, deltas_size );
			const double R = kRValues[table_no];
			deltas = Encoding::ANSDecodeDeltas(deltas_buf, deltas_size, kEntriesPerPark - 1, R);
		}

		uint64_t sum_deltas = 0, sum_stubs = 0;

		for( uint32_t i = 0, start_bit = (3+cur_line_point_size)*8;
				 i < std::min((uint32_t)(position % kEntriesPerPark), (uint32_t)deltas.size());
				 i++) {
			uint64_t stub = cur_stub_size_bits == 0 ? 0 : Util::EightBytesToInt(stubs_buf + start_bit / 8); // seems max k is 56
			stub <<= start_bit % 8;
			stub >>= 64 - cur_stub_size_bits;

			sum_stubs += stub;
			start_bit += cur_stub_size_bits;

			sum_deltas += deltas[i];
		}

		uint128_t big_delta = ((uint128_t)sum_deltas << cur_stub_size_bits) + sum_stubs;

		if( table_no == 0 && bits_cut_no > 0 )
			return RestoreLinePoint( line_point + big_delta );

		return line_point + big_delta;
	}


	uint128_t RestoreLinePoint( uint128_t line_point, bool allow_threads = false ){

		AtomicAdder threads_no( threads_count );
		if( allow_threads && threads_no.inited_at < 8 && bits_cut_no > 14 ){
			std::atomic_bool not_found = true;
			auto zero = std::async( std::launch::async, &Decompressor::RestoreLinePointPart, this, line_point<<1, bits_cut_no - 1, &not_found );
			auto one = std::async( std::launch::async, &Decompressor::RestoreLinePointPart, this, (line_point<<1)+1, bits_cut_no - 1, &not_found );
			uint128_t val = zero.get();
			if( val != 0 ) return val;
			val = one.get();
			if( val != 0 ) return val;
		} else {
			auto lp = RestoreLinePointPart( line_point, bits_cut_no );
			if( lp != 0 ) return lp;
		}

		throw std::runtime_error( "Cannot restore linepoint " + std::to_string((uint64_t)(line_point>>64) )
														 + ":" + std::to_string( (uint64_t)line_point ));
	}

	uint128_t RestoreLinePointPart( uint128_t line_point, uint8_t removed_bits_no, std::atomic_bool *not_found = nullptr ){
		assert( removed_bits_no >= kBatchSizes );
		return RestoreNextLinePoint( line_point << removed_bits_no, removed_bits_no, not_found );
	}

	uint128_t RestoreNextLinePoint( uint128_t line_point, uint8_t removed_bits_no, std::atomic_bool *not_found = nullptr ){
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
		while( left_to_do > 0 && ( not_found == nullptr || not_found->load(std::memory_order_relaxed) ) ){
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

				if( f.FindMatches( bucket_L, bucket_R, nullptr, nullptr ) == 1){
					if( not_found != nullptr )
						not_found->store( false, std::memory_order_relaxed );

					return Encoding::SquareToLinePoint( x1x2.first, b+i );
				}
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
