#include <cstring>
#include <iostream>
#include <netinet/in.h>
#include <sys/socket.h>
#include <fcntl.h>
#include <thread>

#include "cmp_tools.hpp"
#include "util.hpp"
#include "bits.hpp"

namespace TCompress{
const uint32_t PROTOCOL_VER = 0x00020001;
const inline uint8_t NET_PING = 1, NET_PING_RESPONSE = 2, NET_REQUEST_RESTORE = 3,
		NET_RESTORED = 4, NET_NOT_RESTORED = 5, NET_REQUEST_SECOND_ROUND = 6;

void getData( int32_t socket, uint8_t * buf, int32_t size, uint32_t timeoutMs = -1 ){
	assert( size > 0 );
	int res = recv( socket, buf, size, MSG_WAITALL ); // TODO timeout
	if( res != size ){
		close( socket );
		throw std::runtime_error( "Cannot fully send data " + std::to_string(res) + " of " + std::to_string( size )
														 + " error " + std::to_string(errno) + " - " + ::strerror(errno)  );
	}
}

void sendData( int32_t socket, uint8_t *buf, int32_t size, uint32_t timeoutMs = -1 ){
	assert( size > 0);
	int res = send( socket, buf, size, 0 );
	if( res != size ){
		close( socket );
		throw std::runtime_error( "Incorrect size send try to send " + std::to_string(size) + " bytes result " + std::to_string(res)
														 + " error " + std::to_string(errno) + " - " + ::strerror(errno) );
	}
}

struct BufValuesReader{
	const uint8_t *buf;
	const uint32_t buf_size;

	BufValuesReader( const uint8_t *buf, uint32_t buf_zie ) : buf(buf), buf_size(buf_zie){}

	uint128_t Next( uint8_t bits_num ){
		uint128_t val = Util::SliceInt64FromBytes( buf + next_bit/8, next_bit%8, bits_num);
		next_bit += bits_num;
		if( (next_bit+7)/8 > buf_size )
			throw std::runtime_error( "Read data over buffer size" );
		return val;
	}
private:
	uint32_t next_bit = 0;
};

class Server{
	int serverSocket;
	std::thread *accepting_thread;
	bool accepting = true;

public:
	template <typename F>
	Server( uint16_t port, const F& onAccept ){
		// creating socket
		serverSocket = socket( AF_INET, SOCK_STREAM, 0);

		// specifying the address
		sockaddr_in serverAddress;
		serverAddress.sin_family = AF_INET;
		serverAddress.sin_port = htons(port);
		serverAddress.sin_addr.s_addr = INADDR_ANY;

		// binding socket.
		bind( serverSocket, (struct sockaddr*)&serverAddress, sizeof(serverAddress) );

		// listening to the assigned socket
		listen( serverSocket, 5 /*backlog size*/ );
		int status = fcntl( serverSocket, F_SETFL, fcntl(serverSocket, F_GETFL, 0) | O_NONBLOCK);

		if( status < 0 )
			std::cerr << "Server: cannot set socket to non-blocking mode: " << status << std::endl;


		accepting_thread = new std::thread( [this,&onAccept](){
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
								onAccept( clientSocket );
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
		std::cout << "Server stop accepting" << std::endl;
		accepting_thread->join();


		// closing the socket.
		close(serverSocket);

		delete accepting_thread;
	}
};



void StartClient( uint32_t addr, uint16_t port, uint32_t threads_no = THREADS_PER_LP ){
	int clientSocket = socket(AF_INET, SOCK_STREAM, 0);

	// specifying address
	sockaddr_in serverAddress;
	serverAddress.sin_family = AF_INET;
	serverAddress.sin_port = htons(port);
	serverAddress.sin_addr.s_addr = addr;

	// sending connection request
	auto cres = connect(clientSocket, (struct sockaddr*)&serverAddress, sizeof(serverAddress));
	if( cres < 0 )
		throw std::runtime_error( "Cannot connect to server " + std::to_string( cres ) );

	std::cout << "Connected to server with result " << cres << std::endl;
	// sending protocol version
	sendData( clientSocket, (uint8_t*)&PROTOCOL_VER, 4 );

	uint8_t buf[0xffff];
	getData( clientSocket, buf, 4 );
	if( memcmp( &PROTOCOL_VER, buf, 4 ) != 0 )
		throw std::runtime_error( "Incorrect protocol version from server " + Util::HexStr( buf, 4 )
														 + " while expected " + Util::HexStr( (uint8_t*)&PROTOCOL_VER, 4 ) );

	while( true ){
		getData( clientSocket, buf, 3 );

		if( buf[0] == NET_PING ){ // ping request
			buf[0] = NET_PING_RESPONSE; // ping response
			sendData( clientSocket, buf, 3 );
		} else {
			const uint16_t size = (((uint16_t)buf[1])<<8) + buf[2];

			if( size < 34 || size > 65000 ) { // k + cut + plot_id = 34 bytes at least
				close( clientSocket );
				throw std::runtime_error( "Incorrect request size " + std::to_string(size) );
			}
			getData( clientSocket, buf+3, size );
			uint8_t k_size = buf[3], removed_bits_no = buf[4], plot_id[32];
			memcpy( plot_id, buf+5, 32 );
			BufValuesReader buf_reader( buf + 38, size - 35 );

			if( buf[0] == NET_REQUEST_RESTORE ){ // restore
				std::cout << "Resotre request for k: " << (int)k_size << ", plot_id: " << Util::HexStr( plot_id, 32 )
									<< ", removed_bits: " << (int)removed_bits_no << std::endl;

				const uint32_t lp_size = k_size*2 - removed_bits_no;

				Table2MatchData mdata;
				LinePointMatcher validator( k_size, plot_id );
				LinePointInfo left_lp;
				if( buf[37] == 255 ){
					left_lp.full_line_point = buf_reader.Next( k_size*2 );
					mdata.left.Add( -1, left_lp.full_line_point, validator.CalculateYs(left_lp.full_line_point) );
				}

				if( buf[37] < 255 ){
					left_lp.skip_points = buf[37];
					left_lp.orig_line_point = buf_reader.Next( lp_size );
				}

				uint16_t line_points_count = (uint16_t)buf_reader.Next( 16 );
				std::vector<uint128_t> line_points(line_points_count);
				for( uint32_t i = 0; i < line_points_count; i++ )
					line_points[i] = buf_reader.Next( lp_size );

				std::vector<uint16_t> rejects;
				uint16_t res_size;
				ParkBits bits;

				if( mdata.AddBulk( left_lp, 0, line_points, removed_bits_no, rejects, validator, threads_no ) ){
					res_size = 34;
					buf[0] = NET_RESTORED; // restored
					bits.AppendValue( mdata.matched_left, k_size*2 );
					bits.AppendValue( mdata.matched_right_idx, 16 );
					bits.AppendValue( mdata.matched_right, k_size*2 );
				}else {
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
				bits.AppendValue( rejects.size(), 16 );
				memcpy( buf + res_size + 3, rejects.data(), rejects.size()*sizeof(uint16_t) );

				bits.ToBytes( buf + 3 );
				res_size = (bits.GetSize()+7)/8 + rejects.size()*sizeof(uint16_t);

				buf[1] = res_size>>8; // size - high
				buf[2] = res_size; // size - low
				sendData( clientSocket, buf, res_size + 3 );
			} else if( buf[0] == NET_REQUEST_SECOND_ROUND ){ // Second round

			}
		}
	}
	// closing socket
	close(clientSocket);
}
} // namespace
