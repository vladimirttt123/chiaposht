#ifndef SRC_CPP_CMP_NETWORK_HPP_
#define SRC_CPP_CMP_NETWORK_HPP_

#ifdef TCOMPRESS_SERVER_PORT
#define TCOMPERESS_WITH_NETWORK
#endif

#ifdef TCOMPERESS_WITH_NETWORK
#include <cstring>
#include <iostream>
#include <netinet/in.h>
#include <sys/socket.h>
#include <fcntl.h>
#include <thread>
#include <chrono>
using namespace std::chrono_literals; // for operator""min;

#include "cmp_tools.hpp"
#include "util.hpp"
#include "bits.hpp"

namespace TCompress{

#ifdef TCOMPRESS_SERVER_PORT
const uint32_t DEFAULT_SERVER_PORT = TCOMPRESS_SERVER_PORT;
#endif

// client should be quick enough to provide respoinse in this time
uint32_t CLIENT_TIMEOUT_MS = 2000;

const uint32_t PROTOCOL_VER = 0x00020001;
const inline uint8_t NET_PING = 1, NET_PING_RESPONSE = 2, NET_REQUEST_RESTORE = 3,
		NET_RESTORED = 4, NET_NOT_RESTORED = 5, NET_REQUEST_SECOND_ROUND = 6;

void getData( int32_t socket, uint8_t * buf, int32_t size, int32_t timeoutMs = -1 ){
	assert( size > 0 );

	auto err = [&socket, &size]( int res ){
		try{ close( socket ); } catch(...) {}
		std::cout << "Cannot fully get data " + std::to_string(res) + " of " + std::to_string( size )
										 + " error " + std::to_string(errno) + " - " + ::strerror(errno) << std::endl;
		throw std::runtime_error( "Cannot fully get data " + std::to_string(res) + " of " + std::to_string( size )
															 + " error " + std::to_string(errno) + " - " + ::strerror(errno)  );
		};
	if( timeoutMs < 0 ){
		int res = recv( socket, buf, size, MSG_WAITALL ); // TODO timeout
		if( res != size ) err(res);
		return;
	}

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	for( int32_t got = 0; got != size; ){
		int res = recv( socket, buf + got, size - got, MSG_DONTWAIT );
		if( res < 0 && errno != EAGAIN && errno != EWOULDBLOCK )
			err(res);
		if( res > 0 )	got += res;
		if( got != size ){
			std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
			auto responseTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
			if( responseTimeMs > timeoutMs ){
				std::cout << "Remove slow connection. Response time: " << responseTimeMs << "ms" << std::endl;
				throw std::runtime_error( "Slow connection. get data " + std::to_string(res) + " of " + std::to_string( size )
																 + " error " + std::to_string(errno) + " - " + ::strerror(errno)  );
			}

			std::this_thread::sleep_for(10ms);
		}
	}
}

void sendData( int32_t socket, uint8_t *buf, int32_t size, uint32_t timeoutMs = -1 ){
	assert( size > 0);
	int res = send( socket, buf, size, 0 );
	if( res != size ){
		try{ close( socket ); } catch(...){};
		throw std::runtime_error( "Incorrect size send try to send " + std::to_string(size) + " bytes result " + std::to_string(res)
														 + " error " + std::to_string(errno) + " - " + ::strerror(errno) );
	}
}

struct BufValuesReader{
	const uint8_t *buf;
	const uint32_t buf_size;

	BufValuesReader( const uint8_t *buf, uint32_t buf_zie ) : buf(buf), buf_size(buf_zie){}

	uint128_t Next( uint8_t bits_num ){
		uint128_t val = Util::SliceInt128FromBytes( buf + next_bit/8, next_bit%8, bits_num);
		next_bit += bits_num;
		if( (next_bit+7)/8 > buf_size )
			throw std::runtime_error( "Read data over buffer size" );
		return val;
	}
private:
	uint32_t next_bit = 0;
};

template <typename T>
class Server{
	int serverSocket = -1;
	std::thread *accepting_thread;
	bool accepting = true;
	T& accepter;

public:
	Server( uint16_t port, T& accepter ) :accepter(accepter){
		if( port == 0 ) return;

		// creating socket
		serverSocket = socket( AF_INET, SOCK_STREAM, 0);

		// specifying the address
		sockaddr_in serverAddress;
		serverAddress.sin_family = AF_INET;
		serverAddress.sin_port = htons(port);
		serverAddress.sin_addr.s_addr = INADDR_ANY;

		// binding socket.
		int status = bind( serverSocket, (struct sockaddr*)&serverAddress, sizeof(serverAddress) );
		if( status < 0 ){
			std::cout << "Server: cannot bind socket on port " << port << " with " << status
								<< ", errno " << errno << " - " << ::strerror(errno) << std::endl;
			try{ close(serverSocket); } catch(...){}
			return;
		}
		// listening to the assigned socket
		status = listen( serverSocket, 5 /*backlog size*/ );
		if( status < 0 ){
			std::cout << "Server: cannot start listen on port " << port << " with " << status
								<< ", errno " << errno << " - " << ::strerror(errno) << std::endl;
			try{ close(serverSocket); } catch(...){}
			return;
		}

		status = fcntl( serverSocket, F_SETFL, fcntl(serverSocket, F_GETFL, 0) | O_NONBLOCK);

		if( status < 0 )
			std::cerr << "Server: cannot set socket to non-blocking mode: " << status << std::endl;


		accepting_thread = new std::thread( [this, port](){
			std::cout<< "Server accepting on " << port << std::endl;
			while( accepting ){
				try{
					// accepting connection request
					int clientSocket = accept( serverSocket, nullptr, nullptr );
					if( accepting ) {
						if( clientSocket < 0 ){
							if( !( errno == EAGAIN || errno == EWOULDBLOCK ) )
								std::cerr << "ERROR on accept remote connection " << clientSocket << " errno " << errno << std::endl;
							std::this_thread::sleep_for( 1s ); // possible to change accepting to atomic and wait 1 sec or for change
						}
						else{
							std::cout << "Connection achieved" << std::endl;
							uint8_t buf[4];
							sendData( clientSocket, (uint8_t*)&PROTOCOL_VER, 4 );
							getData( clientSocket, buf, 4, 1000 );
							if( memcmp( &PROTOCOL_VER, buf, 4 ) == 0 ){
								std::cout << "Connection accpeted" << std::endl;
								this->accepter.AddRemoteSource( clientSocket );
							}
							else{
								std::cerr << "ERROR protocol agreement with client";
								close( clientSocket );
							}
						}
					}
				}catch(...){
					std::cerr << "Exception on accept remote conection" << std::endl;
				}
			}
			std::cout << "Server Stopped" << std::endl;
		});
	}

	~Server(){
		accepting = false;
		// closing the socket.
		close(serverSocket);

		if( accepting_thread ){
			std::cout << "Server stop accepting" << std::endl;
			accepting_thread->join();

			delete accepting_thread;
		}
	}
};


void StartClient( uint32_t addr, uint16_t port, uint32_t threads_no = THREADS_PER_LP ){
	int clientSocket = socket( AF_INET, SOCK_STREAM, 0 );

	// specifying address
	sockaddr_in serverAddress;
	serverAddress.sin_family = AF_INET;
	serverAddress.sin_port = htons(port);
	serverAddress.sin_addr.s_addr = addr;

	auto err = [&clientSocket]( const std::string &msg ){
		try{ close(clientSocket); } catch (...) {};
		throw std::runtime_error( msg );
	};

	// sending connection request
	auto cres = connect( clientSocket, (struct sockaddr*)&serverAddress, sizeof(serverAddress));
	if( cres < 0 )
		err( "Cannot connect to server result " + std::to_string( cres )
														 + " errno " + std::to_string(errno) + " - " + ::strerror(errno) );

	std::cout << "Connected to server with result " << cres << std::endl;
	// sending protocol version
	sendData( clientSocket, (uint8_t*)&PROTOCOL_VER, 4 );

	uint8_t buf[0xffff];
	getData( clientSocket, buf, 4 );
	if( memcmp( &PROTOCOL_VER, buf, 4 ) != 0 )
		err( "Incorrect protocol version from server " + Util::HexStr( buf, 4 )
														 + " while expected " + Util::HexStr( (uint8_t*)&PROTOCOL_VER, 4 ) );

	while( true ){
		getData( clientSocket, buf, 3, 40*60*1000 /*40 minutes*/ );

		if( buf[0] == NET_PING ){ // ping request
			buf[0] = NET_PING_RESPONSE; // ping response
			sendData( clientSocket, buf, 3 );
		} else {
			const uint16_t size = (((uint16_t)buf[1])<<8) + buf[2];

			if( size < 34 || size > 65000 )  // k + cut + plot_id = 34 bytes at least
				err( "Incorrect request size " + std::to_string(size) );

			getData( clientSocket, buf+3, size );
			uint8_t k_size = buf[3], removed_bits_no = buf[4], plot_id[32];
			memcpy( plot_id, buf+5, 32 );
			BufValuesReader buf_r( buf + 37, size - 34 );

			Table2MatchData mdata;
			LinePointMatcher validator( k_size, plot_id );
			if( buf[0] == NET_REQUEST_RESTORE ){ // restore
				Timer req_timer;
				std::cout << "Resotre request for k: " << (int)k_size << ", plot_id: " << Util::HexStr( plot_id, 32 )
									<< ", removed_bits: " << (int)removed_bits_no << std::flush;

				const uint32_t lp_size = k_size*2 - removed_bits_no;

				LinePointInfo left_lp;
				left_lp.skip_points = buf_r.Next(8);
				if( left_lp.skip_points == 255 ){
					left_lp.full_line_point = buf_r.Next( k_size*2 );
					mdata.left.Add( -1, left_lp.full_line_point, validator.CalculateYs(left_lp.full_line_point) );
				} else {
					left_lp.orig_line_point = buf_r.Next( lp_size );
				}

				uint16_t line_points_count = (uint16_t)buf_r.Next( 16 );
				std::cout << ", points_no: " << line_points_count << std::endl;
				if( line_points_count > 4095 )
					err( "Too many line points " + std::to_string( line_points_count ) );

				std::vector<uint128_t> line_points(line_points_count);
				for( uint32_t i = 0; i < line_points_count; i++ )
					line_points[i] = buf_r.Next( lp_size );

				std::vector<uint16_t> rejects;
				uint16_t res_size;
				LargeBits bits;

				if( mdata.AddBulk( left_lp, 0, line_points, removed_bits_no, rejects, validator, threads_no ) ){
					std::cout << " -> restored";
					buf[0] = NET_RESTORED; // restored
					bits.AppendValue( mdata.matched_left, k_size*2 );
					bits.AppendValue( mdata.matched_right_idx, 16 );
					bits.AppendValue( mdata.matched_right, k_size*2 );
				}else {
					std::cout << " -> calculated";
					buf[0] = NET_NOT_RESTORED; // not restored

					bits.AppendValue( mdata.right.size + 1, 16 ); // number of points to return
					bits.AppendValue( -1, 16 ); // left line point
					bits.AppendValue( mdata.left.points[0].LinePoint(), k_size*2 );

					for( uint32_t i = 0; i < mdata.right.size; i++ ){
						bits.AppendValue( mdata.right.points[i].orig_idx, 16 );
						bits.AppendValue( mdata.right.points[i].LinePoint(), k_size*2 );
					}
				}

				// append rejects
				bits.AppendValue( rejects.size(), 12 );
				for( uint16_t i = 0; i < rejects.size(); i++ )
					bits.AppendValue( rejects[i], 12 );

				bits.ToBytes( buf + 3 );
				res_size = (bits.GetSize()+7)/8;

				buf[1] = res_size>>8; // size - high
				buf[2] = res_size; // size - low
				sendData( clientSocket, buf, res_size + 3 );

				req_timer.PrintElapsed( " in " );
			} else if( buf[0] == NET_REQUEST_SECOND_ROUND ){ // Second round
				Timer req_timer;
				std::cout << "Second round request k: " << (int)k_size << ", plot_id: " << Util::HexStr( plot_id, 32 )
									<< ", removed_bits: " << (int)removed_bits_no << std::endl;
				uint16_t idx = (uint16_t)buf_r.Next(16);
				if( idx != (uint16_t)-1 )
					throw std::runtime_error( "Incorrect left line point index" );

				Table2MatchData mdata;
				LinePointMatcher validator( k_size, plot_id );

				uint128_t lp = buf_r.Next( k_size * 2);
				mdata.left.Add( idx, lp, validator.CalculateYs(lp) );
				for( int32_t rcount = (uint16_t)buf_r.Next(16); rcount > 0; rcount-- ){
					idx = (uint16_t)buf_r.Next(16);
					lp = buf_r.Next( k_size*2 );
					mdata.right.Add( idx, lp, validator.CalculateYs(lp));
				}

				if( mdata.RunSecondRound( removed_bits_no, validator, threads_no ) ){
					std::cout << " - restored ";
					LargeBits bits;
					buf[0] = NET_RESTORED; // restored
					bits.AppendValue( mdata.matched_left, k_size*2 );
					bits.AppendValue( mdata.matched_right_idx, 16 );
					bits.AppendValue( mdata.matched_right, k_size*2 );
					uint16_t res_size = (bits.GetSize()+7)/8;
					bits.ToBytes( buf + 3 );

					buf[1] = res_size>>8; // size - high
					buf[2] = res_size; // size - low
					sendData( clientSocket, buf, res_size + 3 );

				}else
				{
					buf[0] = NET_NOT_RESTORED;
					buf[1] = 0;
					buf[2] = 0;
					sendData( clientSocket, buf, 3 );
				}
				req_timer.PrintElapsed( " in " );
			}
		}
	}
	// closing socket
	close(clientSocket);
}
} // namespace

#endif // TCOMPERESS_WITH_NETWORK
#endif // SRC_CPP_CMP_NETWORK_HPP_
