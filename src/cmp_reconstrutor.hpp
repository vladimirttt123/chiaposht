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



class ReconstructorsManager{
	std::atomic_int32_t free_locals_count;
	std::mutex mut;
	std::vector<int> clients;

public:

	struct SourceDataLocker {
		const int32_t socket = -1;

		explicit SourceDataLocker( std::atomic_int32_t *src ) : release_to(src){}
		SourceDataLocker( int32_t clienSocket, ReconstructorsManager *mngr )
				: socket(clienSocket), mngr(mngr){}

		inline bool isLocal() const { return socket < 0; }
		inline void setBadConnection() {
			bad_connection = true;
			try{ close( socket ); } catch(...){}
		}

		~SourceDataLocker() {
			if( isLocal() ){
				release_to->fetch_add(1, std::memory_order_relaxed );
				release_to->notify_all();
			} else {
				if( !bad_connection ) mngr->AddRemoteSource( socket );
			}
		}

	private:
		bool bad_connection = false;
		std::atomic_int32_t *release_to = nullptr;
		ReconstructorsManager * mngr = nullptr;
	};


	ReconstructorsManager( uint32_t locals_count ):free_locals_count(locals_count){}

	uint32_t getClientsCount() const { return clients.size(); }

	void AddRemoteSource( int clientSocket ){
		std::lock_guard<std::mutex> lk(mut);
		clients.push_back( clientSocket );
	}

	SourceDataLocker getSource(){
		while( true ){
			{
				std::lock_guard<std::mutex> lk(mut);
				if( clients.size() > 0 ){
					int32_t sock = clients[clients.size()-1];
					clients.pop_back();
					return SourceDataLocker( sock, this );
				}
			}
			if( free_locals_count.fetch_sub( 1, std::memory_order_relaxed ) > 0 )
				return SourceDataLocker( &free_locals_count );
			auto old = free_locals_count.fetch_add( 1, std::memory_order_relaxed );
			if( old < 0 )
				free_locals_count.wait( old + 1, std::memory_order::relaxed );
			//std::this_thread::sleep_for( 10ns ); // TODO implement by atomic wait with exists resources amount but than need to upgrade C++ version
		}
	}
};

ReconstructorsManager RManager(2);

// --------------------------------------------------------
// Class to run point reconstructions locally and remotely
// --------------------------------------------------------
struct Reconstructor{
	const uint8_t k_size, *plot_id, removed_bits_no;
	const uint64_t table2_pos, table1_init_pos;
	LinePointInfo left_lp;
	Table2MatchData match_data;
	uint16_t processed_lps_count = 0;
	std::vector<uint128_t> to_process_lp;
	LinePointMatcher validator;
	std::vector<uint16_t> rejects;

	// this constructor could block up to achieving reconstraction source
	Reconstructor( uint8_t k_size, const uint8_t * plot_id, uint8_t removed_bits_no, uint64_t table2_pos,
								uint64_t table1_pos, uint128_t left_line_point )
			: k_size(k_size), plot_id(plot_id), removed_bits_no(removed_bits_no)
			, table2_pos(table2_pos), table1_init_pos(table1_pos), validator( k_size, plot_id )
	{
		if( left_line_point ) match_data.AddLeft( left_line_point, validator );
	}

	inline void setLeftLinePoint( LinePointInfo left_line_point ) {
		assert( left_lp.orig_line_point == 0 );
		left_lp = left_line_point;
	}

	inline bool isNeedMoreLinePoints() const { return (to_process_lp.size()+processed_lps_count) < kEntriesPerPark; }
	inline void addRightLinePoint( uint128_t line_point ){ to_process_lp.push_back(line_point);	}

	// Running line point restore by provided data
	bool Run(){
		if( to_process_lp.size() == 0 ) return match_data.matched_left != 0; // nothing to do - no points were added since last time

		auto rSrc  = RManager.getSource();
		if( rSrc.isLocal() ){
			if( match_data.AddBulk( left_lp, processed_lps_count, to_process_lp, removed_bits_no, rejects, validator ) )
				return true;
		} else {
			// // Prepare buffer to send
			uint8_t buf[0xffff];
			buf[0] = NET_REQUEST_RESTORE; // request for restore
			buf[3] = k_size;
			buf[4] = removed_bits_no;
			memcpy( buf+5, plot_id, 32 );
			ParkBits bits;
			if( match_data.left.size ){
				bits.AppendValue( 255, 8 ); // uncompressed left line point
				bits.AppendValue( match_data.left.points[0].LinePoint(), k_size*2 );
			} else {
				bits.AppendValue( left_lp.skip_points, 8 );
				bits.AppendValue( left_lp.orig_line_point, k_size*2 - removed_bits_no );
			}

			bits.AppendValue( to_process_lp.size(), 16 );
			for( uint32_t i = 0; i < to_process_lp.size(); i++ )
				bits.AppendValue( to_process_lp[i], k_size*2 - removed_bits_no );

			uint16_t buf_size = 34 + (bits.GetSize()+7)/8;
			bits.ToBytes( buf + 37 );
			buf[1] = buf_size >> 8;
			buf[2] = buf_size;

			uint16_t rsize;
			try{
				sendData( rSrc.socket, buf, buf_size  + 3 );
				getData( rSrc.socket, buf, 3, 1000 );
				rsize = (((uint16_t)buf[1])<<8) + buf[2];
				getData( rSrc.socket, buf + 3, rsize );
			} catch(...){
				rSrc.setBadConnection();
				return Run(); // recursion to try another source
			}

			BufValuesReader r( buf+3, rsize );
			if( buf[0] == NET_RESTORED ){
				match_data.matched_left = r.Next( k_size*2 );
				match_data.matched_right_idx = (uint16_t)r.Next(16);
				match_data.matched_right = r.Next( k_size*2 );
				// TODO check is it really match
			} else if( buf[0] == NET_NOT_RESTORED ){
				uint16_t num_points = (uint16_t)r.Next(16);
				uint16_t idx = (uint16_t)r.Next(16);
				if( idx != (uint16_t)-1 ){
					std::cerr << "Incorrect left point index " << idx << std::endl;
					rSrc.setBadConnection();
					return Run();
				}
				uint128_t lp = r.Next(k_size*2);
				match_data.left.Add( idx, lp, validator.CalculateYs( lp ) );
				for( uint16_t i = 0; i < num_points; i++ ){
					idx = (uint16_t)r.Next(16);
					lp = r.Next( k_size*2 );
					match_data.right.Add( idx, lp, validator.CalculateYs(lp) );
				}

			} else {
				std::cerr << "Incorrect client response type " << (int)buf[0] << std::endl;
				rSrc.setBadConnection();
				return Run(); // recursion to try another source
			}
			// TODO read rejects
		}

		processed_lps_count += to_process_lp.size();
		to_process_lp.clear();

		return match_data.matched_left != 0;
	}

	bool RunSecondRound(){
		// TODO add remote posibility
		return match_data.RunSecondRound( removed_bits_no, validator );
	}

	~Reconstructor(){
		// TODO show rejected
	}

};

}
#endif // SRC_CPP_CMP_RECONSTRUCTOR_HPP_
