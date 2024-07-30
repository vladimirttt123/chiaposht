#ifndef SRC_CPP_CMP_RECONSTRUCTOR_HPP_
#define SRC_CPP_CMP_RECONSTRUCTOR_HPP_

#include <vector>
#include <thread>
#include <chrono>
using namespace std::chrono_literals; // for operator""min;



#include "util.hpp"
#include "pos_constants.hpp"
#include "entry_sizes.hpp"
#include "encoding.hpp"
#include "cmp_tools.hpp"
#include "cmp_network.hpp"

namespace TCompress{



// class ReconstructorsManager{
// 	std::atomic_int32_t free_locals_count;
// 	std::mutex mut;
// 	std::vector<int> clients;

// public:

// 	struct SourceDataLocker {
// 		const int32_t socket = -1;
// 		bool bad_connection = false;

// 		explicit SourceDataLocker( std::atomic_int32_t *src ) : release_to(src){}
// 		SourceDataLocker( int32_t clienSocket, ReconstructorsManager *mngr )
// 				: socket(clienSocket), mngr(mngr){}

// 		inline bool isLocal() const { return socket < 0; }
// 		~SourceDataLocker() {
// 			if( isLocal() ){
// 				release_to->fetch_add(1, std::memory_order_relaxed );
// 				release_to->notify_all();
// 			} else {
// 				if( bad_connection )
// 					try{ close( socket ); } catch(...){}
// 				else
// 					mngr->AddRemoteSource( socket );
// 			}
// 		}

// 	private:
// 		std::atomic_int32_t *release_to = nullptr;
// 		ReconstructorsManager * mngr = nullptr;
// 	};


// 	ReconstructorsManager( uint32_t locals_count ):free_locals_count(locals_count){}

// 	uint32_t getClientsCount() const { return clients.size(); }
// 	void AddRemoteSource( int clientSocket ){
// 		std::lock_guard<std::mutex> lk(mut);
// 		clients.push_back( clientSocket );
// 	}

// 	SourceDataLocker getSource(){
// 		while( true ){
// 			{
// 				std::lock_guard<std::mutex> lk(mut);
// 				if( clients.size() > 0 ){
// 					int32_t sock = clients[clients.size()-1];
// 					clients.pop_back();
// 					return SourceDataLocker( sock, this );
// 				}
// 			}
// 			if( free_locals_count.fetch_sub( 1, std::memory_order_relaxed ) > 0 )
// 				return SourceDataLocker( &free_locals_count );
// 			auto old = free_locals_count.fetch_add( 1, std::memory_order_relaxed );
// 			if( old < 0 )
// 				free_locals_count.wait( old + 1, std::memory_order::relaxed );
// 			//std::this_thread::sleep_for( 10ns ); // TODO implement by atomic wait with exists resources amount but than need to upgrade C++ version
// 		}
// 	}
// };

// ReconstructorsManager RManager(2);

// // --------------------------------------------------------
// // Class to run point reconstructions locally and remotely
// // --------------------------------------------------------
// class Reconstructor{
// 	const uint8_t k_size, *plot_id, removed_bits_no;
// 	Table2MatchData match_data;
// 	LinePointInfo src_left_line_point;
// 	std::vector<uint128_t> src_line_points;
// 	uint16_t processed_line_points = 0;
// 	LinePointMatcher validator;
// 	F1Calculator f1;
// 	FxCalculator f;

// public:
// 	// this constructor could block up to achieving reconstraction source
// 	Reconstructor( uint8_t k_size, const uint8_t * plot_id, uint8_t removed_bits_no, uint128_t left_line_point )
// 			: k_size(k_size), plot_id(plot_id), removed_bits_no(removed_bits_no)
// 			, validator( k_size, plot_id ), f1( k_size, plot_id ), f( k_size, 2 )
// 	{
// 		if( left_line_point ) match_data.AddLeft( left_line_point, validator );
// 	}

// 	inline void setLeftLinePoint( LinePointInfo left_line_point ) {
// 		assert( src_left_line_point.orig_line_point == 0 );
// 		src_left_line_point = left_line_point;
// 	}

// 	inline void addRightLinePoint( uint128_t line_point ){
// 		src_line_points.push_back(line_point);
// 	}

// 	inline uint128_t left_LP() const { return match_data.matched_left; }
// 	inline uint128_t right_LP() const {return match_data.matched_right; }
// 	inline uint16_t right_LP_idx() const { return match_data.matched_right_idx; }
// 	inline uint16_t get_processed_count() const { return processed_line_points; }

// 	// Running line point restore by provided data
// 	bool Run( uint64_t table2_pos, uint64_t table1_start_pos ){
// 		if( processed_line_points == src_line_points.size() ) return match_data.matched_left != 0; // nothing to do - no points were added since last time
// 		std::atomic_uint_fast16_t next_point_idx = processed_line_points;

// 		auto rSrc  = RManager.getSource();
// 		if( rSrc.isLocal() ){
// 			auto thread_func = [this, &next_point_idx, &table1_start_pos, &table2_pos](){
// 				LinePointMatcher validator( k_size, plot_id );
// 				F1Calculator f1( k_size, plot_id );
// 				FxCalculator f( k_size, 2 );

// 				for( uint16_t i = next_point_idx.fetch_add(1, std::memory_order_relaxed);
// 						 match_data.matched_left == 0/*not found yet*/ && i < src_line_points.size();
// 						 i = next_point_idx.fetch_add(1, std::memory_order_relaxed) ){
// 					if( (i+1) < src_line_points.size() && src_line_points[i] == src_line_points[i+1] ) continue; // equally cut points evaluated in the same thread
// 					auto restored_lp = RestoreLinePoint( src_line_points[i], f1, f );
// 					//if( CheckRestored( restored_lp, x1x2.second + i, position, false ) ){
// 					if( restored_lp == 0 ){
// 						std::cout << "WARNING!!! Cannot restore line point at pos " << (table1_start_pos+i) << " from table2 pos " << table2_pos << " - SKIPPING" << std::endl;
// 					}else {
// 						match_data.AddRight( i, restored_lp, validator );
// 						if( i > 0 && src_line_points[i-1] == src_line_points[i] )
// 							while(  match_data.matched_left == 0/*not found yet*/ && (restored_lp = FindNextLinePoint( restored_lp+1, f1, f ) ) )
// 								match_data.AddRight( --i, restored_lp, validator );
// 					}
// 				};
// 			};

// 			// check if we need to restore left line point
// 			std::unique_ptr<std::thread,ThreadDeleter> lp1_thread;
// 			if( match_data.left.size == 0 && src_left_line_point.orig_line_point != 0 ) {
// 				lp1_thread.reset( new std::thread( [this, &table2_pos](){
// 					uint128_t lp = RestoreLinePoint( src_left_line_point.orig_line_point, f1, f );
// 					for( uint i = 0; i < src_left_line_point.skip_points; i++ )
// 						lp = FindNextLinePoint( lp + 1, f1, f );
// 					if( lp == 0 ) {
// 						std::cout << "WARNING!!! Cannot restore left line point at pos " << src_left_line_point.position << " from table2 pos " << table2_pos << " - SKIPPING" << std::endl;
// 					}
// 					match_data.AddLeft( lp, validator );
// 				} ) );
// 			}
// 			std::unique_ptr<std::thread, ThreadDeleter> threads[THREADS_PER_LP];

// 			for( uint32_t i = 0; i < THREADS_PER_LP; i++ )
// 				threads[i].reset( new std::thread( thread_func ) );

// 		} else {
// 			// Prepare buffer to send
// 			uint8_t buf[0xffff];
// 			buf[0] = 3; // request for restore
// 			buf[3] = k_size;
// 			buf[4] = removed_bits_no;
// 			memcpy( buf+5, plot_id, 32 );
// 			ParkBits bits;
// 			if( match_data.left.size ){
// 				bits.AppendValue( 255, 8 ); // uncompressed left line point
// 				bits.AppendValue( match_data.left.points[0].LinePoint(), k_size*2 );
// 			} else {
// 				bits.AppendValue( src_left_line_point.skip_points, 8 );
// 				bits.AppendValue( src_left_line_point.orig_line_point, k_size*2 - removed_bits_no );
// 			}

// 			bits.AppendValue( next_point_idx, 16 );
// 			while( next_point_idx < src_line_points.size() )
// 				bits.AppendValue( src_left_line_point[next_point_idx++], k_size*2 - removed_bits_no );

// 			uint16_t buf_size = 34 + (bits.GetSize()+7)/8;
// 			bits.ToBytes( buf + 37 );
// 			buf[1] = buf_size >> 8;
// 			buf[2] = buf_size;

// 			try{ sendData( rSrc.socket, buf, buf_size  + 3 ) } catch(...){
// 				rSrc.bad_connection = true;
// 				return Run( table2_pos, table1_start_pos );
// 			}

// 			// recieve results
// 			throw std::runtime_error("not implemented" );
// 		}

// 		processed_line_points = src_line_points.size();

// 		return match_data.matched_left != 0;
// 	}

// 	bool RunSecondRound( uint64_t table2_pos, uint64_t table1_start_pos ){
// 		// TODO allow by resource manager.

// 		// check more left points
// 		for( uint128_t lp = FindNextLinePoint( match_data.left.points[0].LinePoint() + 1, f1, f );
// 				 lp != 0; lp = FindNextLinePoint( lp + 1, f1, f ) ){
// 			match_data.AddLeft( lp, validator );
// 			if( match_data.matched_left != 0 ) return true;// return if match found
// 		}
// 		// TODO make by threads
// 		// check more right points - it could be long process and seems need threads!!!
// 		for( uint32_t size = match_data.right.size, i = 0; i < size; i++ )
// 			for( uint128_t lp = FindNextLinePoint( match_data.left.points[i].LinePoint() + 1, f1, f);
// 					 lp != 0; lp = FindNextLinePoint( lp + 1, f1, f ) ){
// 				match_data.AddRight( match_data.right.points[i].orig_idx, lp, validator );
// 				if( match_data.matched_left != 0 ) return true;// return if match found
// 			}

// 		return match_data.matched_left != 0;
// 	}
// private:
// 	inline uint128_t RestoreLinePoint( uint128_t cutted_line_point, F1Calculator &f1, FxCalculator &f ) const {
// 		return FindNextLinePoint( cutted_line_point<<removed_bits_no, f1, f );
// 	}

// 	uint128_t FindNextLinePoint( uint128_t line_point, F1Calculator &f1, FxCalculator &f ) const {
// 		assert( removed_bits_no >= kBatchSizes );

// 		auto x1x2 = Encoding::LinePointToSquare( line_point );
// 		uint64_t firstY = f1.CalculateF( Bits( x1x2.first, k_size ) ).GetValue();
// 		uint64_t firstYkBC = firstY / kBC;

// 		const uint64_t BATCH_SIZE = 1U << kBatchSizes; //256
// 		uint64_t batch[BATCH_SIZE];

// 		std::vector<PlotEntry> bucket_L(1);
// 		std::vector<PlotEntry> bucket_R(1);

// 		uint64_t b = x1x2.second, left_to_do = (1ULL << removed_bits_no) - (line_point&((1ULL << removed_bits_no)-1));
// 		while( left_to_do > 0 ){
// 			f1.CalculateBuckets( b, BATCH_SIZE, batch );
// 			uint64_t cur_batch_size = std::min( left_to_do, std::min( BATCH_SIZE, x1x2.first - b ) ); // possible to minimize with left_to_do
// 			for( uint32_t i = 0; i < cur_batch_size; i++ ){
// 				uint64_t cdiff = firstYkBC - batch[i] / kBC;
// 				if( cdiff == 1 ){
// 					bucket_L[0].y = batch[i];
// 					bucket_R[0].y = firstY;
// 				}else if( cdiff == (uint64_t)-1 ){
// 					bucket_L[0].y = firstY;
// 					bucket_R[0].y = batch[i];
// 				}else continue;

// 				if( f.FindMatches( bucket_L, bucket_R, nullptr, nullptr ) == 1)
// 					return Encoding::SquareToLinePoint( x1x2.first, b+i );
// 			}
// 			left_to_do -= cur_batch_size;
// 			b += cur_batch_size;
// 			if( b == x1x2.first ){
// 				firstY = f1.CalculateF( Bits( ++(x1x2.first), k_size ) ).GetValue();
// 				firstYkBC = firstY / kBC;
// 				b = 0;
// 			}
// 		}

// 		return 0; // NOT FOUND
// 	}
// };




}
#endif // SRC_CPP_CMP_RECONSTRUCTOR_HPP_
