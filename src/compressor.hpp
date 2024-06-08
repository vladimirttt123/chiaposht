#ifndef SRC_CPP_COMPRESSOR_HPP_
#define SRC_CPP_COMPRESSOR_HPP_

#include <cstdint>
#include <cstring>
#ifndef _WIN32
#include <unistd.h>
#endif
#include <string>
#include <iostream>
#include <fstream>

#include "util.hpp"
#include "pos_constants.hpp"
#include "entry_sizes.hpp"

class Compressor {
public:
	explicit Compressor( const std::string& filename ){
		this->filename = filename;
		std::cout << "create compressor " << filename << std::endl;

		disk_file = new std::ifstream( filename, std::ios::in | std::ios::binary );
		if (!disk_file->is_open()) {
			delete disk_file;
			disk_file = nullptr;
			throw std::invalid_argument( "Cannot open file " + filename);
		}
		// 19 bytes  - "Proof of Space Plot" (utf-8)
		// 32 bytes  - unique plot id
		// 1 byte    - k
		// 2 bytes   - format description length
		// x bytes   - format description
		// 2 bytes   - memo length
		// x bytes   - memo
		const int BUF_SIZE = 500;
		uint8_t buf[BUF_SIZE];

		Read( buf, 56 + kFormatDescription.size() );
		if( memcmp( buf, "Proof of Space Plot", 19 ) != 0
				|| Util::TwoBytesToInt( buf + 52) != kFormatDescription.size()
				|| memcmp( buf + 54, kFormatDescription.c_str(), kFormatDescription.size() ) ){
			std::cout << "Incorrect plot format. Support only for original non compressed plots" << std::endl;
			throw std::invalid_argument( "Incorrect plot format" );
		}
		memcpy( plotid, buf+19, 32 );
		k_size = buf[51];
		std::cout << "k: " << (int)k_size << "; id: " << Util::HexStr( plotid, 32 ) <<  std::endl;

		memo_size = Util::TwoBytesToInt( buf+54 + kFormatDescription.size() );
		if( (56 + kFormatDescription.size() + memo_size) > BUF_SIZE ){
			std::cout << "Too big memo " << memo_size;
			throw std::runtime_error( "Too big memo" );
		}

		memo = new uint8_t[memo_size];
		Read( memo, memo_size );

		// now table pointers
		std::cout << "Tables pointers: ";
		Read( buf, 80 );
		for( int i = 0; i < 10; i++ ){
			table_pointers[i] = Util::EightBytesToInt( buf + i*8 );
			std::cout << table_pointers[i] << ", ";
		}
		std::cout << std::endl;
	}

	void CompressTo( const std::string& filename, int level ){
		std::cout << "Start compressing with level " << level << " to " << filename << std::endl;

		uint64_t saved_space = 0;
		for( int i = 1; i <= 7; i++)
			saved_space += AnalizeTable(i);

		std::cout << "saved space: " << saved_space << "( " << (saved_space>>20) << "MiB )"<< std::endl;
	}

	~Compressor(){
		if( disk_file != nullptr ) delete disk_file;
		if( memo != nullptr ) delete [] memo;
	}

private:
	std::string filename;
	std::ifstream *disk_file;
	uint8_t k_size;
	uint8_t plotid[32];
	uint16_t memo_size;
	uint8_t * memo = nullptr;
	uint64_t table_pointers[10];


	uint64_t AnalizeTable( int table_no ){
		const uint32_t line_point_size = EntrySizes::CalculateLinePointSize(k_size);
		const uint32_t stub_size = EntrySizes::CalculateStubsSize(k_size);
		uint32_t delta_size = EntrySizes::CalculateMaxDeltasSize(k_size, table_no );
		uint32_t park_size = EntrySizes::CalculateParkSize( k_size, table_no );

		std::cout << "Table " << table_no << ": line_point_size = " << line_point_size
							<< "; stub size = " << stub_size << "; delta_size = " << delta_size
							<< "; park_size = " << park_size
							<< std::endl;

		uint64_t stats[delta_size];
		memset( stats, 0, delta_size*sizeof(uint64_t) );
		Seek( table_pointers[table_no-1] );

		uint8_t buf[park_size+7];
		int64_t number_of_parks = (table_pointers[table_no]-table_pointers[table_no-1])/park_size;
		uint64_t line_point_sum = 0, line_point_min = -1, line_point_max = 0;

		for( int64_t i = 0; i < number_of_parks; i++ ){
			Read( buf, park_size );

			if( i > 0 ){
				uint32_t line_point_top = bswap_32( ((uint32_t*)buf)[0] )/i;
				line_point_sum += line_point_top;
				if( line_point_max < line_point_top ) line_point_max = line_point_top;
				if( line_point_min > line_point_top ) line_point_min = line_point_top;
			}

			int non_zero_at_end = park_size-1;
			while( buf[non_zero_at_end] == 0 )
				non_zero_at_end--;
			stats[non_zero_at_end-stub_size-line_point_size]++;
		}
		std::cout << "line_point: { " << line_point_min << " : " << line_point_sum/number_of_parks
							<< " : " << line_point_max << "}"<< std::endl;
		std::cout<<"Stats: "<<std::flush;
		for( int s = 0; s < delta_size; s++ )
			if( stats[s] != 0 ) std::cout << " [" << s << "->" << stats[s] << "]" << std::flush;
		std::cout << std::endl;

		// evaluate compaction
		uint32_t new_park_size = park_size-line_point_size-stub_size;
		uint32_t overdraft_parks_no = 0, overdraft_size = 0;
		while( (overdraft_parks_no + stats[new_park_size-1]) < 0xffff ){
			overdraft_parks_no += stats[--new_park_size];
			overdraft_size += overdraft_parks_no;
		}

		new_park_size += line_point_size + stub_size;

		uint64_t new_table_size = (number_of_parks - overdraft_parks_no) * new_park_size + overdraft_size +
															overdraft_parks_no*6 /* 2 bytes index and 4 bites pointer to each one */;
		uint64_t saved_space = park_size*number_of_parks-new_table_size;
		std::cout << "overflowed parks no = " << overdraft_parks_no << "; overflowed size = " << overdraft_size << std::endl;
		std::cout << "park size: " << park_size << " -> " << new_park_size << std::endl;
		std::cout << "Table size: " << park_size*number_of_parks
							<< " - " << new_table_size << " = " << saved_space << std::endl;

		return saved_space;
	}

	void Read( uint8_t * buffer, uint64_t size ){
		int64_t pos = disk_file->tellg();
		disk_file->read(reinterpret_cast<char*>(buffer), size);

		if (disk_file->fail()) {
			std::cout << "Failed to read input disk at position " << pos << " size " << size << std::endl;
			throw std::runtime_error("disk read file " +
															 std::to_string(size) + " at position " + std::to_string(pos));
		}
	}

	void Seek( uint64_t position ){
		disk_file->seekg( position );

		if (disk_file->fail()) {
			std::cout << "Disk seek FAILED to position " << position << std::endl;
			throw std::runtime_error("disk seek failed " + std::to_string(position));
		}
	}
};

#endif  // SRC_CPP_COMPRESSOR_HPP_
