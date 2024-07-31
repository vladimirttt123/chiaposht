#ifndef SRC_CPP_DECOMPRESSOR_HPP_
#define SRC_CPP_DECOMPRESSOR_HPP_

#include <atomic>
#include <cstdint>
#include <cstring>
#include <thread>
#ifndef _WIN32
#include <unistd.h>
#endif
#include <string>
#include <iostream>
#include <fstream>

#include "util.hpp"
#include "pos_constants.hpp"
#include "entry_sizes.hpp"
#include "encoding.hpp"
#include "cmp_tools.hpp"
#include "cmp_reconstrutor.hpp"

namespace TCompress {

// this is a cache for restored line points
LinePointsCache LPCache;

class Decompressor{
public:
	Decompressor( const std::string & filename) : filename(filename) {}
	Decompressor(Decompressor& other) noexcept
	{
		filename = other.filename;
		k_size = other.k_size;
		memcpy( plot_id, other.plot_id, 32 );
		bits_cut_no = other.bits_cut_no;
		table2_cut = other.table2_cut;
		improved_file_allign = other.improved_file_allign;
		memcpy( table_pointers, other.table_pointers, 11*sizeof(uint64_t) );
		memcpy( avg_delta_sizes, other.avg_delta_sizes, 7*sizeof(uint64_t) );
		memcpy( parks_counts, other.parks_counts, 7*sizeof(uint32_t) );
		memcpy( min_deltas_sizes, other.min_deltas_sizes, 7*sizeof(uint16_t) );
	}

	// Check if plot represented by this filename is compressed. If it is the case than returnes decompressor for it.
	static Decompressor* CheckForCompressed( const std::string& filename, uint16_t memo_size, uint8_t k, const uint8_t *plot_id ){
		ReadFileWrapper disk_file( filename );
		uint8_t buf[4];
		disk_file.Read( 54, buf, 4 );
		if( memcmp( buf, tFormatDescription.c_str(), 4 ) == 0 ){
			Decompressor *res = new Decompressor( filename );
			res->init( disk_file, memo_size, k, plot_id );
			return res;
		}
		return nullptr;
	}

	// Show info about the current plot
	void ShowInfo( bool is_short = true ){
		if( is_short ) std::cout << program_header << std::endl;

		std::cout << "Compression: " << (table2_cut?"table2 + " : "" ) << (int)bits_cut_no << " bits, "
							<< (improved_file_allign? "improved" : "old") << " file alingnment" << std::endl ;
		if( !is_short ){
			std::cout << "Tables pointers: ";
			for( uint i = 0; i < 10; i++ ) std::cout << table_pointers[i] << (i<9?", ":"");
			if( improved_file_allign ){
				std::cout << std::endl << "Tables io_parameters: ";
				for( uint i = 0; i < 7; i++ ) std::cout << min_deltas_sizes[i] << (i<6?", ":"");
			}
			std::cout <<std::endl;
		}
		std::cout <<std::endl;
	}

	void init( std::ifstream& file, uint16_t memo_size, uint8_t k, const uint8_t *plot_id ){
		ReadFileWrapper disk_file( &file );
		init( disk_file, memo_size, k, plot_id );
	}

	// read the static data from header of compressed plot
	void init( ReadFileWrapper& disk_file, uint16_t memo_size, uint8_t k, const uint8_t *plot_id ){
		k_size = k;
		memcpy( this->plot_id, plot_id, 32 );

		disk_file.Read( 60 + memo_size, (uint8_t*)table_pointers, 80 );
		disk_file.Read( &bits_cut_no, 1 );
		table2_cut = (0x80 & bits_cut_no) != 0;
		improved_file_allign = (0x40 & bits_cut_no ) != 0;
		bits_cut_no = 0x3f & bits_cut_no;

		disk_file.Read( (uint8_t*)parks_counts, 7*4 );
		for( uint i = 0; i < 10; i++ )
			table_pointers[i] = bswap_64(table_pointers[i]);
		table_pointers[10] = disk_file.Size();

		if(improved_file_allign)
			disk_file.Read( (uint8_t*)min_deltas_sizes, 7*2 );

		// now we now how to fill avg_delta_sizes
		for( uint i = 0; i < 7; i++ ){
			parks_counts[i] = bswap_32( parks_counts[i] );
			uint64_t table_size = table_pointers[ i == 6 ? 10 : (i+1)] - table_pointers[ i == 6 ? 9 : i];
			uint64_t static_size = GetMainParkSize( i );
			avg_delta_sizes[i] = (table_size - (static_size+(improved_file_allign?overdraftPointerSize:3))*parks_counts[i])/parks_counts[i];
		}

	}

	uint8_t GetNumberOfRemovedBits() const { return bits_cut_no; }
	bool isImprovedFileFormat() const { return improved_file_allign; }
	bool isTable2Cutted() const { return table2_cut; }
	uint64_t getParksCount( uint8_t table_no ) const { return parks_counts[table_no<6?table_no:6]; }


	uint16_t ReadC3Park( std::ifstream& file, uint64_t idx, uint8_t *buf, uint16_t max_to_read ){
		ReadFileWrapper disk_file( &file );
		return ReadC3Park( disk_file, idx, buf, max_to_read );
	}

	uint16_t ReadC3Park( ReadFileWrapper& disk_file, uint64_t idx, uint8_t *buf, uint16_t max_to_read ){
		if( idx >= parks_counts[6] )
			throw std::runtime_error( "too big C3 park index " + std::to_string(idx)
															 + " when max is " + std::to_string( parks_counts[6] ) );
		uint16_t delta_size = 0;
		uint64_t pos;

		if( improved_file_allign ){
			const uint16_t min_ds = min_deltas_sizes[6];
			const uint32_t main_park_size = min_ds+overdraftPointerSize;
			uint8_t rbuf[main_park_size + overdraftPointerSize];

			disk_file.Read( table_pointers[9] + idx*main_park_size - overdraftPointerSize, rbuf, main_park_size+overdraftPointerSize );
			DeltasStorage::RestoreParkPositionAndSize( true, avg_delta_sizes[6], idx, rbuf, rbuf+main_park_size, pos, delta_size );
			if( delta_size == 0 ){ // No overdraft
				// define real deltas size
				memcpy( buf, rbuf+overdraftPointerSize, delta_size = getNonZerosSize( rbuf+overdraftPointerSize, min_ds ) );
			} else { // Read overdrafted part
				memcpy( buf, rbuf+overdraftPointerSize, min_ds );
				disk_file.Read( table_pointers[9] + main_park_size*(uint64_t)parks_counts[6] + pos, buf+min_ds, delta_size );
				delta_size += min_ds;
			}

		}else {
			uint8_t ibuf[6];
			disk_file.Read( table_pointers[9] + (idx-1)*3, ibuf, 6 );
			DeltasStorage::RestoreParkPositionAndSize( false, avg_delta_sizes[6], idx, ibuf, ibuf+3, pos, delta_size );

			if( delta_size > EntrySizes::CalculateC3Size(k_size) )
				throw std::runtime_error( "too big delta " + std::to_string( delta_size ) + ". possibly corrupted file ");
			if( delta_size <= max_to_read )
				disk_file.Read( table_pointers[9] + parks_counts[6]*3ULL + pos, buf, delta_size );
		}

		return delta_size;
	}

	uint128_t ReadLinePoint( std::ifstream& file, // this need for parallel reading - will support it later
													uint8_t table_no /*0 is a first table*/,
													uint64_t position ){
		uint128_t line_point = LPCache.GetLinePoint( k_size, plot_id, table_no, position );
		if( line_point != 0 ) return line_point;

		switch( table_no ){
		case 0: // The points from first table in case of compression should be taken from cache.
						// beenig in this case means exception that could lead to incorrect proof
			line_point = ReadRealLinePoint( file, table_no, position );
			if( bits_cut_no > 0 ){
				std::cout << "tcompress: Missing line point cache. plot: ( k: " << (int)k_size
									<< ", id: " << Util::HexStr( plot_id, 32 ) << ", table: " << (int)table_no << ", position: "
									<< position << "), cache: (size: " << LPCache.GetSize() << ", count: "
									<< LPCache.GetCount() << ", important count: " << LPCache.GetImportantCount() << ")"
									<< std::endl;
				LPCache.IncreaseSize();
				CheckRestored( line_point = RestoreLinePoint( line_point ), position );
			}
			return line_point; // do not add here to cache sine it could be incorrect

		case 1:
			line_point = ReadRealLinePoint( file, table_no, position );
			if( table2_cut ){
				auto x1x2 = Encoding::LinePointToSquare( line_point << 11 );

				if( x1x2.second + kEntriesPerPark >= x1x2.first )
					throw std::runtime_error( "unsupported case of moving first" ); // this case could be supported too... hope it is rare enough

				uint128_t valid_lp1 = LPCache.GetLinePoint( k_size, plot_id, 0, x1x2.first );
				// if it is possible that left point in cache may be we need to look for right point in cache too?

#define NEW_METHOD_
#ifdef NEW_METHOD_

				// Table2MatchData mdata;
				// LinePointMatcher validator( k_size, plot_id );
				// LinePointInfo left_lp;
				// if( valid_lp1 != 0 )
				// 	mdata.left.Add( -1, valid_lp1, validator.CalculateYs(valid_lp1 ) );
				// else left_lp = ReadLinePointFull( file, 0, x1x2.first );
				// uint16_t added_count = 0;
				// std::vector<uint16_t> rejects;

				// auto show_rejects = [&rejects, &position](){
				// 	if( rejects.size() > 0 ){
				// 		std::cout << "WARNING!!! couldn't restore points on table2 posistion " << position << ": ";
				// 		for( auto v : rejects ) std::cout << v << " ";
				// 		std::cout << std::endl;
				// 	}
				// };
				// auto get_found = [this, &x1x2, &line_point, &position, &mdata, &show_rejects](){
				// 	show_rejects();
				// 	x1x2.second += mdata.matched_right_idx;
				// 	line_point =  Encoding::SquareToLinePoint( x1x2.first, x1x2.second );
				// 	LPCache.AddLinePoint( k_size, plot_id, 0, x1x2.first, mdata.matched_left );
				// 	LPCache.AddLinePoint( k_size, plot_id, 0, x1x2.second, mdata.matched_right );
				// 	return LPCache.AddLinePoint( k_size, plot_id, 1, position, line_point );
				// };

				// while( added_count < kEntriesPerPark ){
				// 	uint16_t pos_in_park = (x1x2.second + added_count)% kEntriesPerPark;
				// 	auto pReader = GetParkReader( file, 0, x1x2.second + added_count, std::max( 2, pos_in_park + 1 ) );
				// 	std::vector<uint128_t> line_points;

				// 	for( line_points.push_back( pReader.NextLinePoint( pos_in_park ) );
				// 			 (line_points.size()+added_count) < kEntriesPerPark && pReader.HasNextInStub(); )
				// 		line_points.push_back( pReader.NextLinePoint() );

				// 	if( mdata.AddBulk( left_lp, added_count, line_points, bits_cut_no, rejects, validator ) )
				// 		return get_found();
				// 	left_lp.orig_line_point = 0; // clear to not add any more
				// 	added_count += line_points.size();

				// 	if( pReader.overdraft_size > 0 && added_count < kEntriesPerPark ){
				// 		// read overdraft and check its points
				// 		ReadFileWrapper disk_file( &file );
				// 		disk_file.Read( pReader.overdraft_pos, pReader.stubs_overdraft_buf(), pReader.overdraft_size );

				// 		line_points.clear();
				// 		while( (added_count+line_points.size() ) < kEntriesPerPark && pReader.HasNext() )
				// 			line_points.push_back( pReader.NextLinePoint() );

				// 		if( mdata.AddBulk( left_lp, added_count, line_points, bits_cut_no, rejects, validator ) )
				// 			return get_found();
				// 		added_count += line_points.size();
				// 	}
				// }
				// if( mdata.RunSecondRound( bits_cut_no, validator ) )
				// 	return get_found();

				// show_rejects();
				// ============================

				Reconstructor rec( k_size, plot_id, bits_cut_no, position, x1x2.second, valid_lp1 );
				if( valid_lp1 == 0 )
					rec.setLeftLinePoint( ReadLinePointFull( file, 0, x1x2.first ) );

				auto get_found = [this, &x1x2, &line_point, &position, &rec](){
					x1x2.second += rec.match_data.matched_right_idx;
					line_point =  Encoding::SquareToLinePoint( x1x2.first, x1x2.second );
					LPCache.AddLinePoint( k_size, plot_id, 0, x1x2.first, rec.match_data.matched_left );
					LPCache.AddLinePoint( k_size, plot_id, 0, x1x2.second, rec.match_data.matched_right );
					return LPCache.AddLinePoint( k_size, plot_id, 1, position, line_point );
				};


				while( rec.isNeedMoreLinePoints() ){
					uint16_t pos_in_park = (x1x2.second + rec.processed_lps_count)% kEntriesPerPark;
					auto pReader = GetParkReader( file, 0, x1x2.second + rec.processed_lps_count, std::max( 2, pos_in_park + 1 ) );

					for( rec.addRightLinePoint( pReader.NextLinePoint( pos_in_park ) );
							 rec.isNeedMoreLinePoints() && pReader.HasNextInStub(); )
						rec.addRightLinePoint( pReader.NextLinePoint() );

					if( rec.Run() ) return get_found();

					if( pReader.overdraft_size > 0 && rec.isNeedMoreLinePoints() ){
						// read overdraft and check its points
						ReadFileWrapper disk_file( &file );
						disk_file.Read( pReader.overdraft_pos, pReader.stubs_overdraft_buf(), pReader.overdraft_size );

						while( rec.isNeedMoreLinePoints() && pReader.HasNext() )
							rec.addRightLinePoint( pReader.NextLinePoint() );

						if( rec.Run() )  return get_found();
					}
				}

				if( rec.RunSecondRound() ) return get_found();
#else


				Table2MatchData match_data;
				LinePointMatcher validator( k_size, plot_id );

				// Read left line point
				std::unique_ptr<std::thread,ThreadDeleter> lp1_thread;
				if( valid_lp1 == 0 )
					lp1_thread.reset( new std::thread( [this, &match_data, &validator, &x1x2]( LinePointInfo lpi ){
									uint128_t lp = RestoreLinePoint( lpi.orig_line_point, validator );
									for( uint i = 0; i < lpi.skip_points; i++ )
										lp = FindNextLinePoint( lp + 1, bits_cut_no, validator );
									CheckRestored( lp, x1x2.first );
									match_data.AddLeft( lp, validator );
						}, ReadLinePointFull( file, 0, x1x2.first )  ) );
				else match_data.AddLeft( valid_lp1, validator );


				// Read right line points
				auto return_found = [&match_data, this, &x1x2, &line_point, &position](){
					x1x2.second += match_data.matched_right_idx;
					line_point =  Encoding::SquareToLinePoint( x1x2.first, x1x2.second );
					LPCache.AddLinePoint( k_size, plot_id, 0, x1x2.first, match_data.matched_left );
					LPCache.AddLinePoint( k_size, plot_id, 0, x1x2.second, match_data.matched_right );
					return LPCache.AddLinePoint( k_size, plot_id, 1, position, line_point );
				};
				auto run_threads = [this, &match_data, &position, &x1x2]( uint128_t line_points[], uint16_t from_idx, uint16_t line_points_count ){
					if( line_points_count == from_idx ) return; // nothing to do
					std::atomic_uint_fast16_t next_point_idx = from_idx;
					auto thread_func = [this, &next_point_idx, &match_data, &line_points_count, &line_points, &position, &x1x2](){
						LinePointMatcher validator( k_size, plot_id );
						for( uint16_t i = next_point_idx.fetch_add(1, std::memory_order_relaxed);
								 match_data.matched_left == 0/*not found yet*/ && i < line_points_count;
								 i = next_point_idx.fetch_add(1, std::memory_order_relaxed) ){
							if( (i+1) < line_points_count && line_points[i] == line_points[i+1] ) continue; // equally cut points evaluated in the same thread
							auto restored_lp = RestoreLinePoint( line_points[i], validator );
							if( CheckRestored( restored_lp, x1x2.second + i, position, false ) ){
								match_data.AddRight( i, restored_lp, validator );
								if( i > 0 && line_points[i-1] == line_points[i] )
									while(  match_data.matched_left == 0/*not found yet*/
												 && (restored_lp = FindNextLinePoint(restored_lp+1, bits_cut_no, validator) ) )
										match_data.AddRight( --i, restored_lp, validator );
							}
						};
					};

					std::unique_ptr<std::thread, ThreadDeleter> threads[THREADS_PER_LP];

					for( uint32_t i = 0; i < THREADS_PER_LP; i++ )
						threads[i].reset( new std::thread( thread_func ) );
				};

				uint128_t line_points[kEntriesPerPark];
				for( uint16_t count_lps = 0; count_lps < kEntriesPerPark; ){
					uint16_t first_lp = count_lps;
					uint16_t pos_in_park = (x1x2.second + count_lps)% kEntriesPerPark;
					auto pReader = GetParkReader( file, 0, x1x2.second + count_lps, std::max( 2, pos_in_park + 1 ) );

					for( line_points[count_lps++] = pReader.NextLinePoint( pos_in_park );
							 count_lps < kEntriesPerPark && pReader.HasNextInStub(); )
						line_points[count_lps++] = pReader.NextLinePoint();

					run_threads( line_points, first_lp, count_lps ); // find in first part
					if( match_data.matched_left != 0 ) return return_found();// return if match found

					first_lp = count_lps;

					if( pReader.overdraft_size > 0 && count_lps < kEntriesPerPark ){
						// read overdraft and check its points
						ReadFileWrapper disk_file( &file );
						disk_file.Read( pReader.overdraft_pos, pReader.stubs_overdraft_buf(), pReader.overdraft_size );

						line_points[count_lps++] = pReader.NextLinePoint( );
						while( count_lps < kEntriesPerPark && pReader.HasNext() )
							line_points[count_lps++] = pReader.NextLinePoint();

						run_threads( line_points, first_lp, count_lps );
						if( match_data.matched_left != 0 ) return return_found();// return if match found
					}
				}

				// check more left points
				for( uint128_t lp = FindNextLinePoint( match_data.left.points[0].LinePoint() + 1, bits_cut_no, validator );
						 lp != 0; lp = FindNextLinePoint( lp + 1, bits_cut_no, validator ) ){
					match_data.AddLeft( lp, validator );
					if( match_data.matched_left != 0 ) return return_found();// return if match found
				}
				// check more right points - it could be long process and seems need threads!!!
				for( uint32_t size = match_data.right.size, i = 0; i < size; i++ )
					for( uint128_t lp = FindNextLinePoint( match_data.left.points[i].LinePoint() + 1, bits_cut_no, validator );
							 lp != 0; lp = FindNextLinePoint( lp + 1, bits_cut_no, validator ) ){
						match_data.AddRight( match_data.right.points[i].orig_idx, lp, validator );
						if( match_data.matched_left != 0 ) return return_found();// return if match found
					}

#endif // NEW_METHOD_
				// IF we are here than we can't restore this line point
				throw std::runtime_error( "Cannot restore line point of table 2 at position " + std::to_string(position) );

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
					CheckRestored( valid_lp1, x1x2.first, position );
					CheckRestored( valid_lp2, x1x2.second, position );

					bool match = CheckMatch( valid_lp1, valid_lp2 );
					while( !match && valid_lp1 ){
						lps1.push_back( valid_lp1 );
						valid_lp1 = FindNextLinePoint( valid_lp1 + 1, bits_cut_no, k_size, plot_id );
						match = valid_lp1 != 0 && CheckMatch( valid_lp1, valid_lp2 );
					}
					while( !match && valid_lp2 != 0 ){
						lps2.push_back( valid_lp2 );
						valid_lp2 = FindNextLinePoint( valid_lp2 + 1, bits_cut_no, k_size, plot_id );
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


	ParkReader GetParkReader( std::ifstream& file, // this need for parallel reading - will support it later
													 uint8_t table_no /*0 is a first table*/,
													 uint64_t position, int16_t need_num_entries = -1 /* -1 is automatic by position */ ){
		ReadFileWrapper disk_file( &file );
		return GetParkReader( disk_file, table_no, position, need_num_entries );
	}

	ParkReader GetParkReader( ReadFileWrapper& disk_file, // this need for parallel reading - will support it later
													 uint8_t table_no /*0 is a first table*/,
													 uint64_t position, int16_t need_num_entries = -1 /* -1 is automatic by position */ ){
		if( table_no >=6 )
			throw std::invalid_argument( "table couldn't be bigger than 5 for reading line point: " + std::to_string(table_no) );

		if( need_num_entries < 0 ) // need to fix number of entries to read
			need_num_entries = (position%kEntriesPerPark)+1;

		const uint64_t park_idx = position/kEntriesPerPark;
		assert( park_idx < parks_counts[table_no] );

		const uint8_t cur_bits_cut = (table_no?0:bits_cut_no) + ((table_no==1&&table2_cut)?11:0);
		const uint32_t cur_line_point_size_bits = k_size*2 - cur_bits_cut;
		const uint32_t cur_stub_size_bits = k_size - kStubMinusBits - cur_bits_cut;
		const uint32_t cur_line_point_size = (cur_line_point_size_bits+7)/8;
		const uint32_t cur_stubs_size = (cur_stub_size_bits*(kEntriesPerPark-1)+7)/8;
		const uint32_t end_pointer_size = improved_file_allign ? overdraftPointerSize : 3;
		const uint64_t main_park_size = cur_line_point_size + cur_stubs_size + end_pointer_size
																		+ ( improved_file_allign ? min_deltas_sizes[table_no] : 0 );
		const uint32_t max_deltas_size = EntrySizes::CalculateMaxDeltasSize( k_size, table_no + 1 );


		if( need_num_entries == 1 ){ // Simplest case read first line point at the begining of park
			uint8_t line_point_buf[cur_line_point_size + 7];
			disk_file.Read( table_pointers[table_no] + main_park_size*park_idx, line_point_buf, cur_line_point_size );
			return ParkReader( line_point_buf, cur_line_point_size_bits );
		}

		uint8_t * full_buf = new uint8_t[main_park_size + end_pointer_size /*because we read 2 deltas pointers*/
																		 + std::max(max_deltas_size, 255U) /*for overdraft or deltas*/];
		uint8_t * line_point_buf = full_buf + end_pointer_size;
		uint8_t * stubs_buf = line_point_buf + cur_line_point_size;
		uint8_t * deltas_buf = improved_file_allign?(line_point_buf+cur_line_point_size):(full_buf + main_park_size);

		disk_file.Read( table_pointers[table_no] + main_park_size*park_idx - end_pointer_size /* to read prev delta pointer */,
									 full_buf, main_park_size + end_pointer_size );

		uint64_t overdraft_pos;
		uint16_t overdraft_size, deltas_size;
		DeltasStorage::RestoreParkPositionAndSize( improved_file_allign, avg_delta_sizes[table_no], park_idx,
																							full_buf, full_buf+main_park_size, overdraft_pos, overdraft_size );

		if( !improved_file_allign ) {
			if( overdraft_size > max_deltas_size )
				throw std::runtime_error( "incorrect deltas size " + std::to_string( overdraft_size ) );

			// all deltas in overdraft in original file format than read them
			disk_file.Read( table_pointers[table_no] + main_park_size*parks_counts[table_no] + overdraft_pos, deltas_buf, overdraft_size );
			deltas_size = overdraft_size;
			overdraft_size = 0; // clear to say parker all read.
		} else {
			stubs_buf += min_deltas_sizes[table_no] + overdraft_size; // now fix start of stubs buf
			// here for improved file format
			if( overdraft_size == 0 ) // no additional read need just evaluate real delstas size
				deltas_size = getNonZerosSize( deltas_buf, min_deltas_sizes[table_no] );
			else {
				deltas_size = min_deltas_sizes[table_no] + overdraft_size;
				overdraft_pos += table_pointers[table_no] + main_park_size*parks_counts[table_no];// fix overdraft position to real in file.

				// time to define do we need to read overdraft or not
				uint32_t real_stubs_size = cur_stubs_size - overdraft_size;
				uint32_t number_of_entries_in_stubs = real_stubs_size*8/cur_stub_size_bits;
				if( (number_of_entries_in_stubs+1/*include first line point*/) < (uint32_t)need_num_entries ) {
					disk_file.Read( overdraft_pos, stubs_buf + cur_stubs_size - overdraft_size, overdraft_size );
					overdraft_size = 0; // clear to say parker all read
				}
			}
		}

		return ParkReader( full_buf, line_point_buf, stubs_buf, deltas_buf, deltas_size,
											cur_line_point_size_bits, cur_stub_size_bits, table_no /* used to unpack deltas */,
											overdraft_pos, overdraft_size );
	}

private:
	uint8_t k_size, plot_id[32];
	uint8_t bits_cut_no;
	bool table2_cut = false, improved_file_allign = false;
	uint64_t table_pointers[11], avg_delta_sizes[7];
	uint32_t parks_counts[7];
	uint16_t min_deltas_sizes[7];
	std::string filename;


	inline bool CheckRestored( uint128_t lp, uint64_t table1_pos, uint64_t table2_pos = 0, bool with_exception = true ){
		if( lp == 0 ){
			auto msg = "Cannot restore line point on table 1 position " + std::to_string(table1_pos)
								 + (table2_pos == 0 ? "" : (", table 2 position " + std::to_string(table2_pos) ) )
								 + ", plotid " + Util::HexStr(plot_id,32);
			if( with_exception ) throw std::runtime_error( msg );

			std::cout << msg << " - SKIPPING " << std::endl;
			return false;
		}
		return true;
	}

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

	// this return size of park without end pointer
	inline uint32_t GetMainParkSize( uint8_t table_no ){
		if( table_no == 6 ) return improved_file_allign ? min_deltas_sizes[6] : 0;

		uint16_t removed = table_no == 0 ? bits_cut_no : ((table_no==1&&table2_cut)?11:0);
		uint16_t line_point_size = ( k_size*2 + 7 - removed)/8;
		uint16_t single_stub_size_bits = k_size - kStubMinusBits - removed;
		uint16_t stubs_size = (single_stub_size_bits*(kEntriesPerPark-1)+7)/8;

		return line_point_size + stubs_size + (improved_file_allign?min_deltas_sizes[table_no]:0);
	}


	static uint16_t getNonZerosSize( uint8_t* buf, uint16_t max_size ){
		while( max_size > 0 && buf[max_size-1] == 0 ) max_size--;
		return max_size;
	}

	uint128_t ReadRealLinePoint( std::ifstream& file, // this need for parallel reading - will support it later
												 uint8_t table_no /*0 is a first table*/,
												 uint64_t position ){
		return ReadLinePointFull( file, table_no, position).orig_line_point;
	}

	LinePointInfo ReadLinePointFull(std::ifstream& file, // this need for parallel reading - will support it later
																	uint8_t table_no /*0 is a first table*/,
																	uint64_t position ){
		const uint64_t pos_in_park = position%kEntriesPerPark;
		auto park = GetParkReader( file, table_no, position );

		LinePointInfo res(position, park.NextLinePoint( pos_in_park ) );
		res.skip_points = park.GetSameDeltasCount();
		return res;
	}

	inline uint128_t RestoreLinePoint( uint128_t cutted_line_point ){
		return FindNextLinePoint( cutted_line_point<<bits_cut_no, bits_cut_no, k_size, plot_id );
	}
	inline uint128_t RestoreLinePoint( uint128_t cutted_line_point, LinePointMatcher &validator ){
		return FindNextLinePoint( cutted_line_point<<bits_cut_no, bits_cut_no, validator );
	}

};

} // endof namespace

#endif  // SRC_CPP_COMPRESSOR_HPP_
