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

class Decompressor{
public:
	Decompressor() {}

	void init( std::ifstream& file, uint16_t memo_size, uint8_t k, const uint8_t *plot_id ){
		k_size = k;
		memcpy( this->plot_id, plot_id, 32 );
		line_point_size = EntrySizes::CalculateLinePointSize(k);
		stubs_size = EntrySizes::CalculateStubsSize(k);
		ReadFileWrapper disk_file( &file );

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
			uint64_t stub = Util::EightBytesToInt(stubs_buf + start_bit / 8); // seems max k is 56
			stub <<= start_bit % 8;
			stub >>= 64 - cur_stub_size_bits;

			sum_stubs += stub;
			start_bit += cur_stub_size_bits;

			sum_deltas += deltas[i];
		}

		uint128_t big_delta = ((uint128_t)sum_deltas << cur_stub_size_bits) + sum_stubs;

		if( table_no == 0 && bits_cut_no > 0 )
			return RestoreLinePoint( line_point + big_delta );
		// if( table_no == -1 ) {
		// 	uint128_t lp = line_point + big_delta;
		// 	checkLinePoint( table_no, position, lp );
		// 	const uint128_t mask = 0xffff;
		// 	const uint64_t shift = k_size-kStubMinusBits-16;
		// 	lp -= lp&(mask<<shift);
		// 	uint32_t match_count = 0;
		// 	for( uint64_t i = 0; i < mask+1; i++ ){
		// 		if( checkLinePoint( table_no, position, lp + (i<<shift) ) )
		// 			match_count ++;
		// 	}
		// 	if( match_count != 1 )
		// 		std::cout << "!!!!!!!!!number of matches: " << match_count << std::endl;
		// }

		return line_point + big_delta;
	}

	bool checkLinePoint( uint8_t table_no, uint64_t position, uint128_t line_point ){

		auto x1x2 = Encoding::LinePointToSquare( line_point );
		F1Calculator f1(k_size, plot_id);
		std::vector<Bits> ys;

		std::pair<Bits, Bits> results = f1.CalculateBucket( Bits(x1x2.second,k_size) );
		ys.push_back(std::get<0>(results));
		results = f1.CalculateBucket( Bits(x1x2.first,k_size) );
		ys.push_back(std::get<0>(results));


		FxCalculator f(k_size, 2);
		std::vector<Bits> new_ys;
		std::vector<Bits> new_metadata;

		PlotEntry l_plot_entry{};
		PlotEntry r_plot_entry{};
		bool is_swap = ys[0].GetValue() > ys[1].GetValue();
		l_plot_entry.y = ys[is_swap?1:0].GetValue();
		r_plot_entry.y = ys[is_swap?0:1].GetValue(); // for this need the alt position?
		std::vector<PlotEntry> bucket_L = {l_plot_entry};
		std::vector<PlotEntry> bucket_R = {r_plot_entry};

		// If there is no match, fails.
		uint64_t cdiff = r_plot_entry.y / kBC - l_plot_entry.y / kBC;
		if (cdiff != 1) {
			return false;
		} else {
			if(f.FindMatches(bucket_L, bucket_R, nullptr, nullptr) != 1) {
				return false;
			}
		}

		return true;
	}

private:
	uint8_t k_size, plot_id[32];
	uint32_t line_point_size, stubs_size, first_table_stubs_size;
	uint8_t bits_cut_no;
	uint64_t table_pointers[11], avg_delta_sizes[7];
	uint32_t parks_counts[7];
#define MAX_RECENT_LINE_POINTS 100U
	std::mutex mut_recent;
	uint128_t recent_line_points[MAX_RECENT_LINE_POINTS];
	uint32_t cur_recent_line_point = 0, recent_line_points_size = 0;
	std::atomic_uint32_t threads_count = 0;

	uint128_t RestoreLinePoint( uint128_t line_point, uint8_t removed_bits_no = 0 ){

		if( removed_bits_no == 0 ) removed_bits_no = bits_cut_no;

		// clear low bits of line_point
		uint128_t base_line_point = line_point << removed_bits_no;

		{
			std::lock_guard<std::mutex> lk(mut_recent);
			for( uint32_t i = 0; i < recent_line_points_size; i++ ){
				if( (recent_line_points[i]>>removed_bits_no) == line_point )
					return recent_line_points[i];
			}
		}

		AtomicAdder threads_no( threads_count );
		if( threads_no.inited_at < 8 && removed_bits_no > 14 ){
			auto zero = std::async( std::launch::async, &Decompressor::RestoreLinePoint, this, line_point<<1, removed_bits_no - 1 );
			auto one = std::async( std::launch::async, &Decompressor::RestoreLinePoint, this, (line_point<<1)+1, removed_bits_no - 1 );
			uint128_t val = zero.get();
			if( val != 0 ) return val;
			val = one.get();
			if( val != 0 ) return val;
			if( removed_bits_no == bits_cut_no )
				throw std::runtime_error( "Cannot restore linepoint " + std::to_string((uint64_t)(base_line_point>>64) )
															 + ":" + std::to_string( (uint64_t)base_line_point ));
			return 0;
		}

		F1Calculator f1( k_size, plot_id );
		std::vector<Bits> ys;

		FxCalculator f(k_size, 2);
		Bits first_bits, second_bits;

		uint64_t last_first = 0;
		// uint128_t res_line_point = 0, res_line_point_diff;

		for( uint64_t i = 0; i < (1UL<<removed_bits_no); i++ ){
			line_point = (base_line_point | (uint128_t)i);

			auto x1x2 = Encoding::LinePointToSquare( line_point );

			std::pair<Bits, Bits> results = f1.CalculateBucket( Bits(x1x2.second,k_size) );
			second_bits = std::get<0>(results);
			if( last_first != x1x2.first ){ // the first value usually not changed and could be not reavaluated
				results = f1.CalculateBucket( Bits(x1x2.first,k_size) );
				first_bits = std::get<0>(results);
				last_first = x1x2.first;
			}

			PlotEntry l_plot_entry{};
			PlotEntry r_plot_entry{};
			bool is_swap = second_bits.GetValue() > first_bits.GetValue();
			l_plot_entry.y = (is_swap?first_bits:second_bits).GetValue();
			r_plot_entry.y = (is_swap?second_bits:first_bits).GetValue();
			std::vector<PlotEntry> bucket_L = {l_plot_entry};
			std::vector<PlotEntry> bucket_R = {r_plot_entry};

			// If there is no match, fails.
			uint64_t cdiff = r_plot_entry.y / kBC - l_plot_entry.y / kBC;
			if( cdiff == 1 && f.FindMatches(bucket_L, bucket_R, nullptr, nullptr) == 1){
				// add to recent line points
				{
					std::lock_guard<std::mutex> lk(mut_recent);
					recent_line_points[(cur_recent_line_point++)%MAX_RECENT_LINE_POINTS] = line_point;
					recent_line_points_size = std::min(cur_recent_line_point, MAX_RECENT_LINE_POINTS);
				}
				return line_point;

				// if(res_line_point == 0){
				// 	res_line_point = line_point;
				// 	res_line_point_diff = r_plot_entry.y - l_plot_entry.y;
				// }
				// else {
				// 	// TODO compare line points?
				// 	if( res_line_point_diff > (r_plot_entry.y - l_plot_entry.y) ){
				// 		res_line_point = line_point;
				// 		res_line_point_diff = r_plot_entry.y - l_plot_entry.y;
				// 	}
				// }
			}
		}

		if( removed_bits_no == bits_cut_no )
			throw std::runtime_error( "cannot restore linepoint " + std::to_string((uint64_t)(base_line_point>>64) )
														 + ":" + std::to_string( (uint64_t)base_line_point ));

		return 0;
	}
};

} // endof namespace

#endif  // SRC_CPP_COMPRESSOR_HPP_
