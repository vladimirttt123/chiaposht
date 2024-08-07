#ifndef SRC_CPP_CMP_RECONSTRUCTOR_HPP_
#define SRC_CPP_CMP_RECONSTRUCTOR_HPP_

#include <vector>

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

	struct ClientError : public std::exception {
		std::string s;
		ClientError(std::string ss) : s(ss) {}
		~ClientError() throw() {}  // Updated
		const char* what() const throw() { return s.c_str(); }
	};

	struct SourceDataLocker {
		const int32_t socket = -1;

		explicit SourceDataLocker( std::atomic_int32_t *src ) : release_to(src){}
		SourceDataLocker( int32_t clienSocket, ReconstructorsManager *mngr )
				: socket(clienSocket), mngr(mngr){}

		inline bool isLocal() const { return socket < 0; }
		inline void setBadConnection() {
			std::cout << "Connection broken" << std::endl;
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

	inline void StartServer( uint16_t port, bool restart = false ){
#ifdef TCOMPERESS_WITH_NETWORK
		std::lock_guard<std::mutex> lk(mut);
		if( !defaultServer || restart )
			defaultServer.reset( new Server( port, *this ) );
#endif // TCOMPERESS_WITH_NETWORK
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

#ifdef TCOMPERESS_WITH_NETWORK
private:
	std::unique_ptr<Server<ReconstructorsManager>> defaultServer;
#endif 	// TCOMPERESS_WITH_NETWORK
};

#ifdef TCOMPRESS_LOCAL_CLIENTS
ReconstructorsManager RManager(TCOMPRESS_LOCAL_CLIENTS);
#else // TCOMPRESS_LOCAL_CLIENTS
ReconstructorsManager RManager(2);
#endif

// --------------------------------------------------------
// Class to run point reconstructions locally and remotely
// --------------------------------------------------------
struct Reconstructor{
	const uint8_t k_size, *plot_id, removed_bits_no;
	const uint64_t table2_pos, table1_init_pos;
	LinePointInfo left_lp;
	Table2MatchData match_data;
	uint16_t first_right_lp_index = 0;
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

	inline uint16_t addRightLinePoint( uint128_t line_point ){
		to_process_lp.push_back(line_point);
		return to_process_lp.size();
	}

	// Running line point restore by provided data
	bool Run(){
		if( to_process_lp.size() == 0 ) return match_data.matched_left != 0; // nothing to do - no points were added since last time

		auto rSrc  = RManager.getSource();
		if( rSrc.isLocal() ){
			if( match_data.AddBulk( left_lp, first_right_lp_index, to_process_lp, removed_bits_no, rejects, validator ) )
				return true;
		} else {
#ifndef TCOMPERESS_WITH_NETWORK
			throw std::runtime_error( "Network is unsupported" );
#else
			try{
				// // Prepare buffer to send
				uint8_t buf[0xffff];
				LargeBits bits;

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


				int32_t dres = SendGet( rSrc.socket, NET_REQUEST_RESTORE, buf, bits );
				if( dres < 0 ) throw ReconstructorsManager::ClientError ( "send/recieve" );

				BufValuesReader r( buf+3, (uint32_t)dres );
				if( buf[0] == NET_RESTORED ){
					processMatchedResponse( r, first_right_lp_index );
				} else if( buf[0] == NET_NOT_RESTORED ){
					uint16_t num_points = (uint16_t)r.Next(16);
					uint16_t idx = (uint16_t)r.Next(16);
					if( idx != (uint16_t)-1 )
						throw ReconstructorsManager::ClientError( "Incorrect left point index " + std::to_string( idx ) );

					uint128_t lp = r.Next(k_size*2);
					if( match_data.left.size == 0 ){
						match_data.left.Add( idx, checkRestore( left_lp.orig_line_point, lp ), validator.CalculateYs( lp ) );
					}
					else {
						if( match_data.left.points[0].LinePoint() != lp )
							throw ReconstructorsManager::ClientError( "Incorrect left line point");
					}
					for( uint16_t i = 0; i < num_points; i++ ){
						idx = (uint16_t)r.Next(16);
						lp = r.Next( k_size*2 );
						match_data.right.Add( idx + first_right_lp_index, checkRestore( to_process_lp[idx], lp ), validator.CalculateYs(lp) );
					}

				} else
					throw ReconstructorsManager::ClientError( "Incorrect client response type " + std::to_string( (int)buf[0] ) );

				// excption from here and below could broke restore because restore data already processed and in
				// read rejects
				uint16_t rejects_no = (uint16_t)r.Next(12);
				for( uint16_t i = 0; i < rejects_no; i++ )
					rejects.push_back( (uint16_t)r.Next(12) + first_right_lp_index );
			} catch( ReconstructorsManager::ClientError &er) {
				std::cout << er.what() << std::endl;
				rSrc.setBadConnection();
				return Run(); // recursion to try another source
			}
#endif // TCOMPERESS_WITH_NETWORK
		}

		to_process_lp.clear();

		return match_data.matched_left != 0;
	}

	bool RunSecondRound(){
		assert( match_data.left.size == 1 );
		assert( match_data.right.size > 0 );

		// TODO add remote posibility may be with slpit for multiple clients
		auto rSrc  = RManager.getSource();
		if( rSrc.isLocal() )
			return match_data.RunSecondRound( removed_bits_no, validator );

#ifndef TCOMPERESS_WITH_NETWORK
		throw std::runtime_error( "Network is unsupported" );
#else
		try{
			uint8_t buf[0xffff];
			LargeBits bits;

			bits.AppendValue( -1, 16 ); // left point index
			bits.AppendValue( match_data.left.points[0].LinePoint(), k_size*2 );

			bits.AppendValue( match_data.right.size, 16 );
			for( uint16_t i = 0; i < match_data.right.size; i++ ){
				bits.AppendValue( match_data.right.points[i].orig_idx, 16 );
				bits.AppendValue( match_data.right.points[i].LinePoint(), k_size*2 );
			}

			int32_t dres = SendGet( rSrc.socket, NET_REQUEST_SECOND_ROUND, buf, bits );
			if( dres < 0 ) throw ReconstructorsManager::ClientError ( "send/recieve" );

			if( buf[0] == NET_NOT_RESTORED )
				return false; // TODO validation, but it need recompute in this case
			if( buf[0] == NET_RESTORED ) {
				BufValuesReader r( buf+3, (uint32_t)dres );
				processMatchedResponse( r, 0 );
				return true;
			} else
				throw ReconstructorsManager::ClientError ( "Incorrect client response type " + std::to_string( (int)buf[0] ) );
		} catch( ReconstructorsManager::ClientError &er) {
			std::cout << er.what() << std::endl;
			rSrc.setBadConnection();
			return Run(); // recursion to try another source
		}
#endif //TCOMPERESS_WITH_NETWORK
	}

	~Reconstructor(){
		// TODO show rejected
	}

private:
#ifdef TCOMPERESS_WITH_NETWORK

	int32_t SendGet( int32_t socket, uint8_t req_type, uint8_t *buf, LargeBits &bits ){
		buf[0] = req_type; // request for restore
		buf[3] = k_size;
		buf[4] = removed_bits_no;
		memcpy( buf+5, plot_id, 32 );

		uint16_t buf_size = 34 + (bits.GetSize()+7)/8;
		bits.ToBytes( buf + 37 );
		buf[1] = buf_size >> 8;
		buf[2] = buf_size;

		uint16_t rsize;
		uint32_t timeout = CLIENT_TIMEOUT_BASE_MS;
		if( removed_bits_no > 10 ) // TODO: add line points count to this function
			timeout += std::min( CLIENT_TIMEOUT_BASE_MS, CLIENT_TIMEOUT_ADD_PER_CUT_BIT << (removed_bits_no-10) );

		try{
			sendData( socket, buf, buf_size  + 3 );
			do{ getData( socket, buf, 3, timeout ); // get header
			} while( buf[0] == NET_PING_RESPONSE || buf[0] == NET_PING ); // skip all pings

			rsize = (((uint16_t)buf[1])<<8) + buf[2];
			getData( socket, buf + 3, rsize, timeout ); // get packet
		} catch(...){
			return -1;
		}

		return rsize;
	}

	void processMatchedResponse( BufValuesReader &r, uint16_t idx_fix ){
		// validate - check is it really match for reals sent points
		Table2MatchData mdata;
		mdata.AddLeft( checkRestore( left_lp.orig_line_point, r.Next( k_size*2 ) ), validator );
		idx_fix += (uint16_t)r.Next(16);
		mdata.AddRight( idx_fix, r.Next( k_size*2 ), validator ); // TODO check restore here? (need to find point with same index
		if( mdata.matched_left == 0 ) throw ReconstructorsManager::ClientError ( "Incorrect restore" );
		// copy validated
		match_data.matched_left = mdata.matched_left;
		match_data.matched_right_idx = mdata.matched_right_idx ;
		match_data.matched_right = mdata.matched_right;
	}

	inline uint128_t checkRestore( uint128_t cut, uint128_t full ){
		if( cut != (full>>removed_bits_no) )
			throw ReconstructorsManager::ClientError ( "Incorrect line point cut" );
		if( FindNextLinePoint( full, removed_bits_no, validator) != full )
			throw ReconstructorsManager::ClientError ( "Incorrect line point retore" );

		return full;
	}
#endif //TCOMPERESS_WITH_NETWORK
};

}
#endif // SRC_CPP_CMP_RECONSTRUCTOR_HPP_
