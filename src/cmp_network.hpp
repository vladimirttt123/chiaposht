#include <cstring>
#include <iostream>
#include <netinet/in.h>
#include <sys/socket.h>
#include <thread>

#include "cmp_reconstrutor.hpp"

namespace TCompress{
// const uint32_t PROTOCOL_VER = 0x00020001;

// void getData( int32_t socket, uint8_t * buf, uint32_t size, uint32_t timeoutMs = -1 ){
// 	int res = recv( socket, buf, size, MSG_WAITALL ); // TODO timeout
// 	if( res != size ){
// 		close( socket );
// 		throw std::runtime_error( "Cannot fully send data " + std::to_string(res) + " of " + std::to_string( size )  );
// 	}
// }

// void sendData( int32_t socket, uint8_t *buf, size_t size, uint32_t timeoutMs = -1 ){
// 	int res = send( socket, buf, size, 0 );
// 	if( res != size ){
// 		close( socket );
// 		throw std::runtime_error( "Incorrect size send try to send " + std::to_string(size) + " bytes result " + std::to_string(res) );
// 	}
// }

// class Server{
// 	int serverSocket;
// 	std::thread *accepting_thread;
// 	bool accepting = true;

// public:
// 	Server( uint16_t port ){
// 		// creating socket
// 		serverSocket = socket( AF_INET, SOCK_STREAM, 0);

// 		// specifying the address
// 		sockaddr_in serverAddress;
// 		serverAddress.sin_family = AF_INET;
// 		serverAddress.sin_port = htons(port);
// 		serverAddress.sin_addr.s_addr = INADDR_ANY;

// 		// binding socket.
// 		bind( serverSocket, (struct sockaddr*)&serverAddress, sizeof(serverAddress) );

// 		// listening to the assigned socket
// 		listen( serverSocket, 5 /*backlog size*/ );

// 		accepting_thread = new std::thread( [this](){
// 			while( accepting ){
// 				try{
// 					// accepting connection request
// 					int clientSocket = accept(serverSocket, nullptr, nullptr);
// 					std::cout << "Connection achieved" << std::endl;
// 					if( accepting && clientSocket < 0 )
// 						std::cerr << "ERROR on accept remote connection" << clientSocket << std::endl;
// 					else{
// 						uint8_t buf[4];
// 						sendData( clientSocket, buf, 4 );
// 						getData( clientSocket, buf, 4, 1000 );
// 						if( memcmp( &PROTOCOL_VER, buf, 4 ) == 0 ){
// 							std::cout << "Connection accpeted" << std::endl;
// 							RManager.AddRemoteSource( clientSocket );
// 						}
// 						else{
// 							std::cerr << "ERROR protocol agreement with client";
// 							close( clientSocket );
// 						}
// 					}
// 				}catch(...){
// 					std::cerr << "Exception on accept remote conection" << std::endl;
// 				}
// 			}
// 			std::cout << "Server Stopped" << std::endl;
// 		});
// 	}

// 	~Server(){
// 		accepting = false;

// 		// closing the socket.
// 		close(serverSocket); // TODO check if it throws in accept

// 		accepting_thread->join();
// 		delete accepting_thread;
// 	}
// };



// void StartClient( uint32_t addr, uint16_t port ){
// 	int clientSocket = socket(AF_INET, SOCK_STREAM, 0);

// 	// specifying address
// 	sockaddr_in serverAddress;
// 	serverAddress.sin_family = AF_INET;
// 	serverAddress.sin_port = htons(port);
// 	serverAddress.sin_addr.s_addr = addr;

// 	// sending connection request
// 	auto cres = connect(clientSocket, (struct sockaddr*)&serverAddress, sizeof(serverAddress));
// 	if( cres < 0 )
// 		throw std::runtime_error( "Cannot connect to server " + std::to_string( cres ) );

// 	std::cout << "Connected to server with result " << std::endl;
// 	// sending protocol version
// 	sendData( clientSocket, (uint8_t*)&PROTOCOL_VER, 4 );

// 	uint8_t buf[0xffff];
// 	getData( clientSocket, buf, 4 );
// 	if( memcmp( &PROTOCOL_VER, buf, 4 ) != 0 )
// 		throw std::runtime_error( "Incorrect protocol version from server " + Util::HexStr( buf, 4 )
// 														 + " while expected " + Util::HexStr( (uint8_t*)&PROTOCOL_VER, 4 ) );

// 	while( true ){
// 		getData( clientSocket, buf, 3 );

// 		if( buf[0] == 1 ){ // ping request
// 			buf[0] = 2; // ping response
// 			sendData( clientSocket, buf, 3 );
// 		} else {
// 			uint16_t size = (((uint16_t)buf[1])<<8) + buf[2];
// 			if( size < 34 || size > 65000 ) { // k + cut + plot_id = 34 bytes at least
// 				close( clientSocket );
// 				throw std::runtime_error( "Incorrect request size " + std::to_string(size) );
// 			}
// 			getData( clientSocket, buf+3, size );
// 			uint8_t k_size = buf[3], removed_bits_no = buf[4], plot_id[32];
// 			memcpy( plot_id, buf+5, 32 );
// 			uint32_t next_bit = 38*8;

// 			if( buf[0] == 3 ){ // restore
// 				uint128_t left_line_point = 0;
// 				if( buf[37] == 255 ){
// 					left_line_point = Util::SliceInt128FromBytes( buf+ next_bit/8, next_bit%8, k_size*2 );
// 					next_bit += k_size*2;
// 				}

// 				Reconstructor rec( k_size, plot_id, removed_bits_no, left_line_point );

// 				const uint32_t lp_size = k_size*2 - removed_bits_no;
// 				if( buf[37] < 255 ){
// 					LinePointInfo lp;
// 					lp.skip_points = buf[37];
// 					lp.orig_line_point = Util::SliceInt128FromBytes( buf + next_bit/8, next_bit%8, lp_size );
// 					rec.setLeftLinePoint( lp );
// 				}

// 				uint16_t right_idx = (uint16_t)Util::SliceInt64FromBytes( buf + next_bit/8, next_bit%8, 16 );
// 				next_bit += 16;
// 				for( uint32_t bits_in_input = size*8; next_bit+lp_size < bits_in_input; next_bit += lp_size )
// 					rec.addRightLinePoint( Util::SliceInt128FromBytes( buf + next_bit/8, next_bit%8, lp_size ) );

// 				if( rec.Run(0,0) ){
// 					buf[0] = 4; // restored
// 					buf[1] = 0; // size - high
// 					buf[2] = 34; // size - low
// 					((uint128_t*)(buf + 3))[0] = rec.left_LP();
// 					((uint16_t*)(buf + 19))[0] = rec.right_LP_idx() + right_idx;
// 					((uint128_t*)(buf + 21))[0] = rec.left_LP();
// 					sendData( clientSocket, buf, 37 );
// 				}else {
// 					buf[0] = 5; // unrestored
// 				}
// 			}
// 		}
// 	}
// 	// closing socket
// 	close(clientSocket);
// }
} // namespace
