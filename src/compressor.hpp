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
#include "encoding.hpp"
#include "disk.hpp"

#define _DEBUG_NO_WRITE_COMPRESSED

namespace TCompress {

// This is replacement made for kFormatDescription in compress plots
const std::string tFormatDescription = "t0.1";
const char* plotMagicFrase = "Proof of Space Plot";

struct OriginalTableInfo{
public:
	const uint64_t table_pointer, table_size;
	const uint8_t k, line_point_size;
	const uint64_t delta_size, stub_size, park_size, parks_count;

	uint64_t total_deltas_size = 0;
	uint32_t *delta_real_sizes;
	int64_t line_point_min_delta = -1, line_point_max_delta = 0, line_point_sum_deltas = 0;

	uint8_t new_line_point_size = 0;
	uint64_t new_park_size = 0, new_table_size = 0;
	uint32_t overflowed_parks_count = 0, overflowed_parks_additional_size = 0;

	OriginalTableInfo( uint8_t k, uint8_t table_no, uint64_t table_pointer, uint64_t table_size )
			: table_pointer(table_pointer), table_size(table_size), k(k)
			, line_point_size( EntrySizes::CalculateLinePointSize(k) )
			, delta_size( EntrySizes::CalculateMaxDeltasSize( k, table_no ) )
			, stub_size( EntrySizes::CalculateStubsSize( k ) )
			, park_size( EntrySizes::CalculateParkSize( k, table_no ) )
			, parks_count( table_size / park_size )
			, delta_real_sizes( new uint32_t[delta_size+1] )
	{
		memset( delta_real_sizes, 0, (delta_size+1)*sizeof(uint32_t) );
	}


	void addStats( int64_t line_point_delta, uint32_t real_park_size ){
		delta_real_sizes[ real_park_size - line_point_size - stub_size ]++;
		total_deltas_size += real_park_size - line_point_size - stub_size;

		line_point_sum_deltas += line_point_delta;
		if( line_point_max_delta < line_point_delta ) line_point_max_delta = line_point_delta;
		if( line_point_delta > 1 && line_point_min_delta > line_point_delta ) line_point_min_delta = line_point_delta;
	}

	void showStats(){
		std::cout << "line_point: { " << line_point_min_delta << " : " << line_point_sum_deltas/parks_count
							<< " : " << line_point_max_delta << "}"<< std::endl;
		std::cout<<"Stats: " << std::flush;
		for( uint32_t s = 0; s < delta_size; s++ )
			if( delta_real_sizes[s] != 0 ) std::cout << " [" << s << "->" << delta_real_sizes[s] << "]" << std::flush;
		std::cout << std::endl;
	}

	void evaluateCompaction(){
		new_park_size = park_size - line_point_size - stub_size;
		new_table_size = table_size;
		uint64_t next_increase = ( overflowed_parks_count + 4*delta_real_sizes[new_park_size-1] );
		while( (overflowed_parks_additional_size + next_increase ) < 0xffffff && parks_count > next_increase ){
//		while( (overflowed_parks_count + delta_real_sizes[new_park_size-1]) < 0xffff ){
			overflowed_parks_count += delta_real_sizes[--new_park_size];
			overflowed_parks_additional_size += overflowed_parks_count;
			next_increase = ( overflowed_parks_count + 4*delta_real_sizes[new_park_size-1] );
		}

		new_line_point_size = line_point_size;// Util::ByteAlign( (uint32_t)k + k - 63 + std::log2(line_point_max_delta-line_point_min_delta) )/8; // TODO new size;
		new_park_size += new_line_point_size + stub_size;

		//new_table_size = new_park_size * parks_count + overflowed_parks_additional_size + overflowed_parks_count*3/*pointers into addittional part*/;
		new_table_size = total_deltas_size + parks_count*(uint64_t)(stub_size+line_point_size+2/*addition per park to link deltas*/);
		std::cout << "deltas_size: " << (parks_count*(uint64_t)delta_size) << " -> " <<  total_deltas_size << std::endl;
		std::cout << "line_point_size: " << (int)line_point_size << " -> " << (int)new_line_point_size << std::endl;
		std::cout << "overflowed parks no = " << overflowed_parks_count << "; overflowed size = " << overflowed_parks_additional_size << std::endl;
		std::cout << "park size: " << park_size << " -> " << new_park_size << std::endl;
		std::cout << "Table size: " << park_size*parks_count
							<< " - " << new_table_size << " = " << (table_size-new_table_size)
							<< " ( " << ((table_size-new_table_size)>>20) << "MiB )" <<  std::endl;
	}

	~OriginalTableInfo(){
		delete [] delta_real_sizes;
	}
};

struct ReadFileWrapper{
public:
	explicit ReadFileWrapper( const std::string& filename )
			: disk_file( new std::ifstream( filename, std::ios::in | std::ios::binary ) ), isExternal(false)
	{
		if( !disk_file->is_open() )
			throw std::invalid_argument( "Cannot open file " + filename);
	}

	explicit ReadFileWrapper( std::ifstream * file )
			: disk_file( file ), isExternal(true)
	{
		if( !disk_file->is_open() )
			throw std::invalid_argument( "File is not opened" );
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

	void Read( uint64_t position, uint8_t * buffer, uint64_t size ){
		Seek( position );
		Read( buffer, size );
	}

	void Seek( uint64_t position ){
		disk_file->seekg( position );

		if( disk_file->fail() ) {
			std::cout << "Disk seek FAILED to position " << position << std::endl;
			throw std::runtime_error("disk seek failed " + std::to_string(position));
		}
	}

	int64_t Size(){
		int64_t pos = disk_file->tellg();
		disk_file->seekg( 0, std::ios_base::end );
		if (disk_file->fail()) {
			std::cout << "Disk seek FAILED to end of file" << std::endl;
			throw std::runtime_error("disk seek failed to end of file" );
		}
		int64_t size = disk_file->tellg();
		Seek( pos );
		return size;
	}

	~ReadFileWrapper(){
		if( !isExternal ) delete disk_file;
	}
private:
	std::ifstream *disk_file;
	bool isExternal;
};

struct DeltasStorage{
public:
	const uint32_t parks_count;

	uint64_t total_size = 0;
	uint16_t *all_sizes = nullptr;

	DeltasStorage( uint64_t parks_count ) : parks_count(parks_count){
		if( parks_count > 0xffffffffUL )
			throw std::runtime_error("too many parks " + std::to_string(parks_count) );

		all_sizes = new uint16_t[ parks_count ];
	}

	void Add( uint64_t idx, uint16_t size ){
		assert( idx < parks_count );
		total_size += (all_sizes[idx] = size);
	}

	void TotalEndToBuf( uint64_t idx, uint8_t *buf ){
		buf[0] = total_size >> ((idx&1)?24:16);
		buf[1] = total_size>>8;
		buf[2] = total_size;
	}

	// check all partially saved position could be restored
	bool IsDeltasPositionRestorable(){
		uint64_t initial_total_size = total_size;
		uint64_t park_avg_size = total_size/parks_count;
		uint8_t buf[6];
		uint8_t *buf_prev = buf, *buf_cur = buf + 3;
		uint64_t delta_pos;
		uint16_t delta_size;

		total_size = 0; // to work with function TotalEndToBuf

		for( uint64_t i = 0; i < parks_count; i++ ){
			total_size += all_sizes[i];
			TotalEndToBuf( i, buf_cur );
			RestoreParkPositionAndSize( park_avg_size, i, buf_prev, buf_cur, delta_pos, delta_size );
			if( delta_pos != (total_size-all_sizes[i]) || delta_size != all_sizes[i] ){
				std::cout << "Couldn't restore parkd index " << i << " with position " << (total_size-all_sizes[i])
									<< " and size " << all_sizes[i] << ". Predicted position " << delta_pos
									<< " predicted size " << delta_size << std::endl;

				total_size = initial_total_size;
				return false;
			}

			buf_prev = buf_cur;
			buf_cur = buf + 3*(i&1);
		}

		return true;
	}

	static void RestoreParkPositionAndSize( uint64_t park_avg_size, uint64_t park_idx, uint8_t *prev_buf, uint8_t *cur_buf,
																				uint64_t &delta_position, uint16_t &encoded_delta_size ){
		uint64_t cur_pos = (((uint16_t)cur_buf[1])<<8) | cur_buf[2];

		if( park_idx == 0 ){
			delta_position = 0;
			encoded_delta_size = cur_pos;
			return;
		}


		uint64_t prev_pos = (((uint16_t)prev_buf[1])<<8) | prev_buf[2];
		encoded_delta_size = cur_pos - prev_pos;

		if( park_idx&1 ){
			prev_pos |= ((uint64_t)prev_buf[0])<<16;
			bool is_jump = (prev_pos+encoded_delta_size)>>24; // is adding current delta will lead to jump to next value in bits24+
			prev_pos |= (((uint64_t)(cur_buf[0]-(is_jump?1:0)) )&0xffUL)<<24;
		}else{
			cur_pos = (cur_pos | 0x1000000 | (((uint64_t)cur_buf[0])<<16)) - encoded_delta_size;
			prev_pos |= (cur_pos&0xff0000) | (((uint64_t)prev_buf[0])<<24);
		}

		uint64_t expected_by_magic = park_avg_size*park_idx;
		if( (expected_by_magic&0xffffffff) < 0x50000000	&& ( prev_pos&0xffffffff ) > 0xa0000000 )
			expected_by_magic -= 0x100000000UL;
		else if( (expected_by_magic&0xffffffff) > 0xa0000000	&& ( prev_pos&0xffffffff ) < 0x50000000 )
			expected_by_magic += 0x100000000UL;
		prev_pos |= expected_by_magic&~0xffffffffUL;

		delta_position = prev_pos;
	}

	~DeltasStorage(){
		if( all_sizes != nullptr ) delete[]all_sizes;
	}

private:
	static uint64_t compound( uint8_t top, uint8_t high, uint8_t med, uint8_t low )
	{
		return (((uint64_t)top)<<24) | (((uint64_t)high)<<16) | (((uint64_t)med)<<8)  | (((uint64_t)low));
	}
};

struct ProgressCounter{
	ProgressCounter( uint64_t total) : scale( total/20 ){}
	void ShowNext( uint64_t current ){
		for( uint64_t cur = current/scale; shown < cur; shown++ )
			std::cout << (((shown+1)%5)?'-':'*') << std::flush;
	}
	~ProgressCounter() { std::cout<<std::endl; }
private:
	uint64_t scale, shown = 0;
};

class Compressor {
public:
	explicit Compressor( const std::string& filename )
			:filename(filename), disk_file( filename )
	{
		std::cout << "create compressor on " << filename << std::endl;

		// 19 bytes  - "Proof of Space Plot" (utf-8)
		// 32 bytes  - unique plot id
		// 1 byte    - k
		// 2 bytes   - format description length
		// x bytes   - format description
		// 2 bytes   - memo length
		// x bytes   - memo
		const int BUF_SIZE = 500;
		uint8_t buf[BUF_SIZE];

		disk_file.Read( buf, 56 + kFormatDescription.size() );
		if( memcmp( buf, plotMagicFrase, 19 ) != 0
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
		disk_file.Read( memo, memo_size );

		// now table pointers
		std::cout << "Tables pointers: ";
		disk_file.Read( buf, 80 );
		for( int i = 0; i < 10; i++ ){
			table_pointers[i] = Util::EightBytesToInt( buf + i*8 );
			std::cout << table_pointers[i] << ", ";
		}
		table_pointers[10] = disk_file.Size();

		std::cout << std::endl;
	}

	void CompressTo( const std::string& filename, uint8_t level ){
		std::cout << "Start compressing with level " << (int)level << " to " << filename << std::endl;

		// for( uint i = 1; i < 7; i++)
		// 	AnalizeLinePoints( i );

		// return;
		// uint64_t saved_space = AnalizeC3Table();
		// for( int i = 1; i < 7; i++){
		// 	auto tinfo = AnalizeTable(i);
		// 	saved_space += tinfo.table_size - tinfo.new_table_size;
		// }

		// std::cout << "saved space: " << saved_space << "( " << (saved_space>>20) << "MiB )"<< std::endl;


		// writing header
		const uint32_t header_size = 60 + memo_size + 10*8/*tables pointers*/ + 1 /*compression level*/ + 7*4 /*compressed tables parks count*/;
		uint8_t header[header_size];
		memcpy( header, plotMagicFrase, 19 );
		memcpy( header + 19, plotid, 32 );
		header[51] = k_size;
		header[52] = 0;
		header[53] = 4; // size of format
		memcpy( header + 54, tFormatDescription.c_str(), 4 );
		header[58] = memo_size>>8;
		header[59] = memo_size;
		memcpy( header + 60, memo, memo_size );
		// next is new tables pointers
		header[140] = level; // compression level

		uint64_t *new_table_pointers = (uint64_t*)(header+60+memo_size);
		new_table_pointers[0] = header_size;
		uint32_t *parks_count = (uint32_t*)(header+60+memo_size + 10*8 + 1 );
		uint64_t total_saved = 0;

		FileDisk output_file( filename );
		output_file.Write( 0, header, header_size ); // could be done at the end when everythins is full;
		for( uint32_t i = 1; i < 7; i++ ){
			std::cout << "Table " << i << ": compacting " << std::flush;
			auto tinfo = CompactTable( i, &output_file, new_table_pointers[i-1] );
			new_table_pointers[i] = new_table_pointers[i-1] + tinfo.new_table_size;
			if( tinfo.parks_count > 0xffffffffUL )
				throw std::runtime_error("too many parks");
			parks_count[i-1] = tinfo.parks_count;

			uint64_t saved = tinfo.table_size - tinfo.new_table_size;
			total_saved += saved;
			std::cout << "		original size: " << tinfo.table_size << " (" << (tinfo.table_size>>20) << "MiB)"
								<< "; compacted size: " << tinfo.new_table_size << " (" << (tinfo.new_table_size>>20) << "MiB)"
								<< "; saved: " << saved << " (" << (saved>>20) << "MiB)" << std::endl;
		}

		std::cout << "copy tables 7, C1, C2 " << std::flush;
#ifndef _DEBUG_NO_WRITE_COMPRESSED
		CopyData( table_pointers[6], table_pointers[9] - table_pointers[6], &output_file, new_table_pointers[6] );
#endif
		new_table_pointers[7] = new_table_pointers[6] + (table_pointers[7]-table_pointers[6]);
		new_table_pointers[8] = new_table_pointers[7] + (table_pointers[8]-table_pointers[7]);
		new_table_pointers[9] = new_table_pointers[8] + (table_pointers[9]-table_pointers[8]);

		std::cout << "Table C3: compacting " << std::flush;
		parks_count[6] = CompactC3Table( &output_file, new_table_pointers[9], total_saved );


		// invert table pointers - TODO check endians and skip if not needed
		for( uint32_t i = 0; i < 10; i++ ){
			new_table_pointers[i] = bswap_64( new_table_pointers[i] ); //
			if( i < 7 )
				parks_count[i] = bswap_32( parks_count[i] );
		}

		// save final header
		output_file.Write( 0, header, header_size ); // could be done at the end when everythins is full;

		std::cout << "Total decrease size: " << total_saved << " (" << (total_saved>>20) << "MiB)" << std::endl;
	}

	~Compressor(){
		if( memo != nullptr ) delete [] memo;
	}

private:
	const std::string filename;
	ReadFileWrapper disk_file;
	uint8_t k_size;
	uint8_t plotid[32];
	uint16_t memo_size;
	uint8_t * memo = nullptr;
	uint64_t table_pointers[11];


	bool CouldRestoreLPbySaved( uint64_t expected_line_point, uint64_t saved_bit_no, uint64_t line_point ){
		uint64_t exp_top = expected_line_point>>saved_bit_no, lp_top = line_point>>saved_bit_no;

		auto at_begin = [ &saved_bit_no](uint64_t v){ return ((v>>(saved_bit_no-4))&0xf) < 6; };
		auto at_end = [ &saved_bit_no](uint64_t v){ return ((v>>(saved_bit_no-4))&0xf) > 9; };

		if( exp_top > 0 && at_begin(expected_line_point) && at_end(line_point) )
			exp_top--;
		else if( at_end(expected_line_point) && at_begin(line_point) )
			exp_top++;

		if( exp_top == lp_top ) return true;

		return false;
	}
	void AnalizeLinePoints( int table_no ){
		std::cout<<"Analyze line point for table " << table_no << std::endl;

		OriginalTableInfo tinfo( k_size, table_no, table_pointers[table_no-1], table_pointers[table_no] - table_pointers[table_no-1] );

		// read last park
		uint8_t line_point_buf[8+7];
		disk_file.Read( tinfo.table_pointer + tinfo.park_size*(tinfo.parks_count-1), line_point_buf, 8 );
		uint64_t line_point = Util::EightBytesToInt( line_point_buf );
		const uint64_t avg_line_point = line_point / (tinfo.parks_count-1);

		uint64_t bits_stats[65];
		memset( bits_stats, 0, 65*8 );

		for( uint64_t i = 0; i < tinfo.parks_count; i++ ){
			// read park's i checkpoint line point
			disk_file.Read( tinfo.table_pointer + tinfo.park_size*i, line_point_buf, 8 );
			line_point = Util::EightBytesToInt( line_point_buf );


			// evaluate how much we could save;
			uint64_t expected_line_point = avg_line_point*i;
			uint64_t need_bits = 4;
			while( !CouldRestoreLPbySaved(expected_line_point, need_bits, line_point ) )
				need_bits++;

			bits_stats[need_bits]++;
		}

		uint64_t bits_sums[65];
		memset( bits_sums, 0, 65*8 );
		std::cout << "Original line_point size: " << (int)tinfo.line_point_size << "; bits_stas: [";
		for( uint i = 4; i < 65; i++ ){
			bits_sums[i] = bits_sums[i-1] + bits_stats[i];

			if( bits_stats[i] != 0 )
				std::cout << i << ": " << bits_stats[i] << " (" << bits_sums[i] << "), ";
		}
		std::cout << "]" << std::endl;
		std::cout << "Savings on 7 byte: " << ( bits_sums[55] - (/*9bytes*/bits_sums[64] - bits_sums[62]) ) << std::endl
							<< "Saving on 6 bytes: " << ( bits_sums[47]/*6bytes*/ + bits_sums[53]/*6+7bytes*/
																					 - (/*9bytes*/bits_sums[64] - bits_sums[61]) ) <<std::endl
							<< "Saving on 5 bytes: "<< ( bits_sums[39]/*5bytes*/ + bits_sums[45]/*5+6bytes*/
																					 + bits_sums[53]/*5+6+7bytes*/ - (/*9 bytes*/ bits_sums[64] - bits_sums[62] ) ) <<std::endl
							<< "Saving on 4 bytes: "<< ( bits_sums[31] + bits_sums[37] + bits_sums[45] + bits_sums[53]
																					 - ( bits_sums[64] - bits_sums[62] ) ) <<std::endl;

	}

	OriginalTableInfo AnalizeTable( int table_no ){
		OriginalTableInfo tinfo( k_size, table_no, table_pointers[table_no-1], table_pointers[table_no] - table_pointers[table_no-1] );
		uint8_t buf[tinfo.park_size+7];

		// read last line point
		disk_file.Read( table_pointers[table_no-1] + (tinfo.park_size*(tinfo.parks_count-1) ), buf, 8 );
		uint64_t last_line_point = bswap_64( ((uint64_t*)buf)[0] );
		uint64_t avg_line_point = last_line_point / (tinfo.parks_count-1);

		std::cout << "Table " << table_no << ": line_point_size = " << (int)tinfo.line_point_size
							<< "; stub size = " << tinfo.stub_size << "; delta_size = " << tinfo.delta_size
							<< "; park_size = " << tinfo.park_size << "; parks_count = " << tinfo.parks_count << std::endl;

		disk_file.Seek( table_pointers[table_no-1] );

		// std::vector<unsigned char> deltas;
		// auto R = kRValues[table_no - 1];
		uint32_t groups_park_size = 0, max_group_park_size = 0, group_size = 8;

		for( int64_t i = 0; i < tinfo.parks_count; i++ ){
			disk_file.Read( buf, tinfo.park_size );

			int32_t encoded_deltas_size = 0x7fff& ( ((uint32_t)buf[tinfo.line_point_size + tinfo.stub_size])
																			| (((uint32_t)buf[tinfo.line_point_size + tinfo.stub_size + 1])<<8) );

			// deltas = Encoding::ANSDecodeDeltas( buf + tinfo.line_point_size + tinfo.stub_size + 2, encoded_deltas_size, kEntriesPerPark - 1, R );
			// uint8_t max_delta = 0;
			// for( int d = 0; d < kEntriesPerPark-1; d++ )
			// 	if( max_delta < deltas[d] ) max_delta = deltas[d];

			uint32_t real_park_size = tinfo.stub_size + tinfo.line_point_size + encoded_deltas_size + 2; // 2 bytes for deltas size

			if( (i%group_size) == 0 ){
				if( max_group_park_size < groups_park_size ) max_group_park_size = groups_park_size;
				groups_park_size = 0;
			}
			groups_park_size += real_park_size;

			uint64_t line_point = bswap_64( ((uint64_t*)buf)[0] );
			tinfo.addStats( line_point - avg_line_point*i, real_park_size );
		}

		std::cout << "max group size: " << max_group_park_size << ", per group size: " << (max_group_park_size/group_size) << std::endl;

		tinfo.showStats();

		tinfo.evaluateCompaction();

		return tinfo;
	}

	void CopyData( uint64_t src_position, uint64_t size,  FileDisk *output_file, uint64_t output_position ){
		const uint64_t BUF_SIZE = 64*1024;
		uint8_t buf[BUF_SIZE];
		ProgressCounter progress( size );


		disk_file.Seek( src_position );
		for( uint64_t i = 0; i < size; i += BUF_SIZE ){
			progress.ShowNext( i );

			uint64_t to_copy = std::min( BUF_SIZE, i - size );
			disk_file.Read( buf, to_copy );
			output_file->Write( output_position + i, buf, to_copy );
		}
	}
	// Compacted table has following structuer:
	// 1. Part with constant size block i.e. stub, line_point, pointer_to_deltas
	// 2. Deltas - writte one after another without spaces deltas.
	// Each entry of 1st part contains of
	// 1. Line Point
	// 2. Stubs
	// 3. 3 bottom bytes of pointer to deltas location after deltas of this part added.
	//    it is supposed that first pointer is 0 and at length of deltas could be evaluated by substracting prev and next pointer.
	OriginalTableInfo CompactTable( int table_no, FileDisk * output_file, uint64_t output_position ){
		OriginalTableInfo tinfo( k_size, table_no, table_pointers[table_no-1], table_pointers[table_no] - table_pointers[table_no-1] );
		uint8_t buf[tinfo.park_size+7], small_buf[8];

		const uint64_t lps_size =  tinfo.line_point_size + tinfo.stub_size;
		uint64_t deltas_position = output_position + tinfo.parks_count * (3 + lps_size);
		DeltasStorage deltas( tinfo.parks_count );

		ProgressCounter progress( tinfo.parks_count );

		disk_file.Seek( table_pointers[table_no-1] ); // seek to start of the table
		for( uint64_t i = 0; i < tinfo.parks_count; i++ ){
			progress.ShowNext( i );

			disk_file.Read( buf, tinfo.park_size );

#ifndef _DEBUG_NO_WRITE_COMPRESSED
			output_file->Write( output_position + (lps_size + 3)*i, buf, lps_size );
#endif
			deltas.Add( i, 0x7fff & ( ((uint32_t)buf[lps_size])
															| (((uint32_t)buf[lps_size + 1])<<8) ) );

#ifndef _DEBUG_NO_WRITE_COMPRESSED
			output_file->Write( deltas_position, buf+lps_size+2, deltas.all_sizes[i] );
#endif
			deltas_position += deltas.all_sizes[i];

			deltas.TotalEndToBuf( i, small_buf ); // this is buf with pointer already to next postion it is OK
#ifndef _DEBUG_NO_WRITE_COMPRESSED
			output_file->Write( output_position + (lps_size + 3)*i + lps_size, small_buf, 3 );
#endif
		}

		tinfo.new_table_size = (lps_size + 3)*tinfo.parks_count + deltas.total_size;

		// check all partially saved position could be restored
		if( !deltas.IsDeltasPositionRestorable() )
			throw std::runtime_error( "couldn't restore position of deltas" );

		return tinfo;
	}


	uint64_t CompactC3Table( FileDisk * output_file, uint64_t output_position, uint64_t &total_saved ){
		auto park_size = EntrySizes::CalculateC3Size(k_size);
		uint64_t parks_count = (table_pointers[10] - table_pointers[9]) / park_size;
		const uint32_t POS_BUF_SIZE = 30000;
		uint8_t buf[park_size], positions[POS_BUF_SIZE];
		uint32_t cur_pos_buf = 0;
		uint64_t next_buf_pos = output_position;
		uint64_t deltas_position = output_position + parks_count*3;

		DeltasStorage deltas(parks_count);
		ProgressCounter progress( parks_count );

		disk_file.Seek( table_pointers[9] );
		for( uint64_t i = 0; i < parks_count; i++ ){
			progress.ShowNext( i );

			disk_file.Read( buf, park_size );
			deltas.Add( i,  Bits(buf, 2, 16).GetValue() );

#ifndef _DEBUG_NO_WRITE_COMPRESSED
			output_file->Write( deltas_position, buf + 2, deltas.all_sizes[i] );
#endif
			deltas_position += deltas.all_sizes[i];

			deltas.TotalEndToBuf( i, positions + cur_pos_buf );
			cur_pos_buf += 3;
			if( cur_pos_buf >= POS_BUF_SIZE ){
				output_file->Write( next_buf_pos, positions, cur_pos_buf );
				next_buf_pos += POS_BUF_SIZE;
				cur_pos_buf = 0;
			}
		}

#ifndef _DEBUG_NO_WRITE_COMPRESSED
		output_file->Write( next_buf_pos, positions, cur_pos_buf ); // write last positions
#endif
		// check all partially saved position could be restored
		if( !deltas.IsDeltasPositionRestorable() )
			throw std::runtime_error( "couldn't restore position of deltas" );

		uint64_t original_size = park_size * parks_count,
				compacted_size = (3*parks_count + deltas.total_size);
		uint64_t saved = original_size - compacted_size ;
		total_saved += saved;
		std::cout << "		original size " << original_size
							<< "; compressed size: " << compacted_size << std::endl;

		return parks_count;
	}
};

} // end of namespace
#endif  // SRC_CPP_COMPRESSOR_HPP_
