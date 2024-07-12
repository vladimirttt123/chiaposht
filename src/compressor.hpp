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
#include "cmp_tools.hpp"
#include "decompressor.hpp"


namespace TCompress {

struct TooSmallMinDeltasException : public std::exception {
	const uint16_t used_min_value, found_value;
	TooSmallMinDeltasException( uint16_t used_value,  uint16_t found_value )
			: used_min_value(used_value), found_value(found_value) {}
	~TooSmallMinDeltasException() throw() {}  // Updated
	const char* what() const throw() {	return "used min deltas sizes value is too small"; }
};

struct OriginalTableInfo{
public:
	const uint64_t table_pointer, table_size;
	const uint8_t k, line_point_size;
	const uint64_t delta_size, single_stub_size_bits, stubs_size, park_size, parks_count;
	const uint64_t new_line_point_size, new_single_stub_size_bits, new_stubs_size;
	const uint64_t single_stub_mask;

	uint64_t new_table_size = 0;
	std::unique_ptr<DeltasStorage> new_deltas;

	OriginalTableInfo( uint8_t k, uint8_t table_no, uint64_t table_pointer, uint64_t table_size, uint32_t stubs_bits_to_remove = 0 )
			: table_pointer(table_pointer), table_size(table_size), k(k)
			, line_point_size( EntrySizes::CalculateLinePointSize(k) )
			, delta_size( EntrySizes::CalculateMaxDeltasSize( k, table_no ) )
			, single_stub_size_bits( k - kStubMinusBits )
			, stubs_size( EntrySizes::CalculateStubsSize( k ) )
			, park_size( EntrySizes::CalculateParkSize( k, table_no ) )
			, parks_count( table_size / park_size )
			, new_line_point_size( (2*k-stubs_bits_to_remove+7)/8 )
			, new_single_stub_size_bits( single_stub_size_bits - stubs_bits_to_remove )
			, new_stubs_size( Util::ByteAlign((kEntriesPerPark - 1) * ( new_single_stub_size_bits )) / 8 )
			, single_stub_mask( ((uint64_t)-1)>>(64-single_stub_size_bits))
			, new_deltas( new DeltasStorage( parks_count ) )
	{}

	OriginalTableInfo( uint8_t k, uint8_t table_no, uint64_t table_pointer, uint64_t table_size, const Decompressor * decompressor )
			: table_pointer(table_pointer), table_size(table_size), k(k)
			, line_point_size( (k*2 - RemovedBitsNumber(table_no, decompressor)  + 7 )/8 )
			, delta_size( EntrySizes::CalculateMaxDeltasSize( k, table_no ) )
			, single_stub_size_bits( k - kStubMinusBits - RemovedBitsNumber(table_no, decompressor) )
			, stubs_size( (single_stub_size_bits*(kEntriesPerPark-1)+7)/8 )
			, park_size( EntrySizes::CalculateParkSize( k, table_no ) )
			, parks_count( decompressor->getParksCount( table_no - 1 ) )
			, new_line_point_size( line_point_size )
			, new_single_stub_size_bits( single_stub_size_bits )
			, new_stubs_size( stubs_size )
			, single_stub_mask( ((uint64_t)-1)>>(64-single_stub_size_bits))
			, new_deltas( new DeltasStorage( parks_count ) )
	{}

	// check all partially saved position could be restored
	std::thread CheckNewDeltas(){
		return std::thread( []( DeltasStorage * deltas){
			if( !deltas->IsDeltasPositionRestorable() ){
				delete deltas;
				throw std::runtime_error( "couldn't restore position of deltas" );
			}
			delete deltas;
		}, new_deltas.release() );
	}
private:
	static uint8_t RemovedBitsNumber( uint8_t table_no, const Decompressor * decompressor ){
		switch(table_no){
		case 1: return decompressor->GetNumberOfRemovedBits();
		case 2: return decompressor->isTable2Cutted()?11:0;
		}
		return 0;
	}
};

struct BufWriter{
	static const uint32_t SIZE = 1024*1024; // 1Mb
	const uint64_t start_write_pos;
	BufWriter( FileDisk *output_file, uint64_t output_pos )
			: start_write_pos(output_pos), output_file(output_file), write_pos( output_pos ){}

	void Write( const uint8_t *next_buf, uint32_t length ){
		assert( length < SIZE );
		if( (buf_length + length) > SIZE )
			Flush();

		memcpy( buf + buf_length, next_buf, length );
		buf_length += length;
	}

	uint64_t Flush() {
		output_file->Write( write_pos, buf, buf_length );
		write_pos += buf_length;
		buf_length = 0;
		return write_pos;
	}

	inline uint64_t GetPosition() const { return write_pos + buf_length; }
	inline uint64_t GetWirttenAmount() const { return GetPosition() - start_write_pos; }

	~BufWriter() {Flush();}
private:
	FileDisk *output_file;
	uint64_t write_pos;
	uint32_t buf_length = 0;
	uint8_t buf[SIZE];
};

struct TableWriter {
public:
	const uint32_t min_deltas_size, stubs_size, park_size, line_point_size;
	const uint64_t output_position, parks_count;

	// line_point_size is in bytes and it transfered as first part of stubs when writing
	// stubs_size is a size of stubs excluding line point size.
	TableWriter( FileDisk * output_file, uint64_t output_position, uint16_t min_deltas_size, uint64_t parks_count, uint32_t line_point_size,
							uint32_t stubs_size, DeltasStorage *deltas_sizes_storage)
			: min_deltas_size( min_deltas_size )
			, stubs_size(stubs_size), park_size( line_point_size + stubs_size + min_deltas_size + overdraftPointerSize)
			, line_point_size(line_point_size) // TODO move park size evaluation to better place
			, output_position(output_position), parks_count(parks_count), output_file(output_file)
			, deltas_sizes_storage(deltas_sizes_storage), stubs_writer( output_file, output_position )
			, overdrafts_writer( output_file, output_position + parks_count*park_size )
	{	}

	void WriteNext( uint8_t * deltas, uint16_t deltas_size ){ // for c3 table
		assert( line_point_size == 0 );
		assert( stubs_size == 0 );

		if( output_file != nullptr ){

			deltas_size = fitDeltasToMin( deltas, deltas_size );

			stubs_writer.Write( deltas, min_deltas_size );
			uint8_t end_buf[overdraftPointerSize];
			deltas_sizes_storage->Add( park_idx, deltas_size - min_deltas_size );
			deltas_sizes_storage->TotalEndToBuf( park_idx, end_buf );
			stubs_writer.Write( end_buf, overdraftPointerSize );

			if( deltas_size > min_deltas_size )
				overdrafts_writer.Write( deltas + min_deltas_size, deltas_size - min_deltas_size );

			assert( overdrafts_writer.GetWirttenAmount() == deltas_sizes_storage->total_size );
			park_idx++;
		}
	}

	void WriteNext( const uint8_t *stubs, uint8_t * deltas, uint64_t deltas_size ){
		WriteNext( stubs, stubs + line_point_size, deltas, deltas_size );
	}

	void WriteNext( const uint8_t * line_point, const uint8_t *stubs, uint8_t * deltas, uint64_t deltas_size ){
		assert( stubs_writer.GetPosition() < output_position + (uint64_t)parks_count*park_size ); // check for too many parks
		assert( park_size >= deltas_size + line_point_size + overdraftPointerSize ); // is deltas fit to park

		if( output_file != nullptr ){
			stubs_writer.Write( line_point, line_point_size ); // first part of stub is the check point line point
			deltas_size = fitDeltasToMin( deltas, deltas_size );
			stubs_writer.Write( deltas, deltas_size ); // write deltas
			uint32_t overdraft_size = deltas_size - min_deltas_size;
			if( overdraft_size > 254 ) throw TooSmallMinDeltasException( min_deltas_size, overdraft_size );
			stubs_writer.Write( stubs, stubs_size - overdraft_size ); // fill free space with stub
			// define what to write to end
			deltas_sizes_storage->Add( park_idx, overdraft_size );
			uint8_t end_buf[overdraftPointerSize];
			deltas_sizes_storage->TotalEndToBuf( park_idx, end_buf );
			stubs_writer.Write( end_buf, overdraftPointerSize );

			if( overdraft_size > 0 )
				overdrafts_writer.Write( stubs + stubs_size - overdraft_size, overdraft_size );

			park_idx++;

			assert( park_idx * (uint64_t)park_size == stubs_writer.GetWirttenAmount() );
		}
	}
	uint64_t getWrittenSize() const { return overdrafts_writer.GetPosition() - output_position; }

	void Flush(){
		stubs_writer.Flush();
		overdrafts_writer.Flush();
	}
private:
	FileDisk *output_file;
	DeltasStorage *deltas_sizes_storage;
	BufWriter stubs_writer, overdrafts_writer;
	uint64_t park_idx = 0;

	inline uint16_t fitDeltasToMin( uint8_t * deltas, uint16_t deltas_size ){
		if( deltas_size >= min_deltas_size ) return deltas_size;

		// WARNING supposed deltas buffer has at least  min_deltas_size bytes
		memset( deltas + deltas_size, 0, min_deltas_size - deltas_size ); // fill deltas with 0 to recieve identical file each time
		return min_deltas_size;  // fix delta size to be at least min_delta_size
	}
};

struct ProgressCounter{
	ProgressCounter( uint64_t total ) : scale( total/20 ){}
	void ShowNext( uint64_t current ){
		for( uint64_t cur = current/scale; shown < cur; shown++ )
			std::cout << (((shown+1)%5)?'-':'*') << std::flush;
	}
	~ProgressCounter() { progres_time.PrintElapsed( "  " ); }
private:
	uint64_t scale, shown = 0;
	Timer progres_time;
};

class Compressor {
public:
	explicit Compressor( const std::string& filename )
			:filename(filename), disk_file( filename )
	{
		std::cout << program_header << std::endl
							<< "Create compressor on " << filename << std::endl;

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
		bool is_correct = memcmp( buf, plotMagicFrase, 19 ) == 0 &&
				 Util::TwoBytesToInt( buf + 52) == kFormatDescription.size();
		if( is_correct && memcmp( buf + 54, tFormatDescription.c_str(), tFormatDescription.size() ) == 0 )
			decompressor.reset( new Decompressor( filename ) );
		else is_correct &= memcmp( buf + 54, kFormatDescription.c_str(), kFormatDescription.size() ) == 0;

		if( !is_correct ) {
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
		disk_file.Read( buf, 80 );
		for( int i = 0; i < 10; i++ ){
			table_pointers[i] = Util::EightBytesToInt( buf + i*8 );
		}
		table_pointers[10] = disk_file.Size();

		if( decompressor ){
			decompressor->init( disk_file, memo_size, k_size, plotid );
			std::cout << "Input file is a compressed file." << std::endl;
			decompressor->ShowInfo( false );
		}
	}

	void CompressTo( const std::string& filename, uint8_t level, uint8_t io_optimized = 0 ){
		uint8_t bits_to_cut = cmp_level_to_bits_cut( level );
		bool cut_table2 = bits_to_cut > 18;
		if( cut_table2 ) bits_to_cut -= 11;
		const int16_t io_p = std::min( 10, (int32_t) io_optimized );

		std::cout << "Start compressing with level " << (int)level << " (" << (cut_table2?"table2 + ":"") << (int)bits_to_cut
							<< "bits) and io_optimization " << io_p << std::endl << "Output file: " << filename
							<< (filename == "plot.dat" ? "\n *** !!! SEEMS YOU MISSING OUTPUT FILE PARAMETER !!! *** " : "" ) << std::endl;

		if( std::filesystem::exists(filename) )
			throw std::invalid_argument( "output file exists please delete manually: " + filename );

		// writing header
		const uint32_t header_size = 60 + memo_size + 10*8/*tables pointers*/ + 1 /*compression level*/ + 7*4 /*compressed tables parks count*/ + 7*2 /*min deltas sizes*/;
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
		// next 80 bytes is new tables pointers
		header[140+memo_size] = bits_to_cut | (cut_table2?0x80:0) | 0x40/*flag of improved align*/; // compression level

		uint64_t *new_table_pointers = (uint64_t*)(header+60+memo_size);
		new_table_pointers[0] = header_size;
		uint32_t *parks_count = (uint32_t*)(header+60+memo_size + 10*8/*table_pointers*/ + 1 /*compression level*/ );
		uint16_t *min_deltas_sizes = (uint16_t*)(header+60+memo_size + 10*8/*table_pointers*/ + 1 /*compression level*/ + 7*4 /*compressed tables parks count*/ );
		// when compressing table2 than in most cases 1 table1 park read up to end that increases IO requests
		min_deltas_sizes[0] = 900 + std::max( 0, io_p-2 )*12;
		// the tables 2-6 usually reads for 1 line point than it has very little impact on IO reqeusts
		min_deltas_sizes[1] = min_deltas_sizes[2] = min_deltas_sizes[3] = min_deltas_sizes[4]
				= 770 + std::max( 0, io_p-4 ) * 7;
		min_deltas_sizes[5] = min_deltas_sizes[1] - 20;
		// by my measures C3 sizes for diffrerent plots varies in range 2300-2460.
		// but c3 park alway read fully than small values bring 2 io request instead of 1
		// Original size for this is 3000 bytes...
		min_deltas_sizes[6] = 2300 + std::min( 160, io_p*22 );

		uint64_t total_saved = 0;
		Timer total_timer;
		std::vector<std::thread> delta_check_threads;
		FileDisk output_file( filename );

		for( uint32_t i = 1; i < 7; i++ ){
			std::cout << "Table " << i << ": " << std::flush;

			auto tinfo = decompressor ? RealignTable( i, &output_file, new_table_pointers[i-1], min_deltas_sizes[i-1] ) :
											 ( ((i == 1 && bits_to_cut > 0 ) || (i==2 && cut_table2) ) ?
															CompressTable( i, &output_file, new_table_pointers[i-1], i==1 ? bits_to_cut : 11, min_deltas_sizes[i-1] )
														: CompactTable( i, &output_file, new_table_pointers[i-1], min_deltas_sizes[i-1] ) );
			new_table_pointers[i] = new_table_pointers[i-1] + tinfo.new_table_size;

			if( tinfo.parks_count > 0xffffffffUL ) throw std::runtime_error("too many parks");
			parks_count[i-1] = tinfo.parks_count;

			tinfo.new_deltas->showStats(); // DEBUG
			total_saved += ShowTableSaved( tinfo.table_size, tinfo.new_table_size );

			tinfo.new_deltas->IsDeltasPositionRestorable();
			delta_check_threads.push_back( tinfo.CheckNewDeltas() );
		}

		std::cout << "Tables 7, C1, C2: copy     " << std::flush;
		CopyData( table_pointers[6], table_pointers[9] - table_pointers[6], &output_file, new_table_pointers[6] );

		new_table_pointers[7] = new_table_pointers[6] + (table_pointers[7]-table_pointers[6]);
		new_table_pointers[8] = new_table_pointers[7] + (table_pointers[8]-table_pointers[7]);
		new_table_pointers[9] = new_table_pointers[8] + (table_pointers[9]-table_pointers[8]);


		std::cout << "Table C3: " << (decompressor? "realign  " : "compacting" ) << "       " << std::flush;
		parks_count[6] = CompactC3Table( &output_file, new_table_pointers[9], total_saved, min_deltas_sizes[6] );


		for( uint32_t i = 0; i < delta_check_threads.size(); i++ ) delta_check_threads[i].join();

		// invert table pointers - TODO check endians and skip if not needed
		for( uint32_t i = 0; i < 10; i++ ){
			new_table_pointers[i] = bswap_64( new_table_pointers[i] ); //
			if( i < 7 )
				parks_count[i] = bswap_32( parks_count[i] );
		}

		// save final header
		output_file.Write( 0, header, header_size ); // done at the end when everythins is full;
		total_timer.PrintElapsed( "Total time: " );
		std::cout << "Original file size: " << (table_pointers[10]>>20) << "MiB; compacted file size: "
							<< ((table_pointers[10]-total_saved)>>20) << "MiB ("
							<< std::fixed << std::setprecision(2)
							<< ((table_pointers[10]-total_saved)*100.0/table_pointers[10]) << "%); saved: "
							<< (total_saved>>20) << "MiB" << std::endl;
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
	std::unique_ptr<Decompressor> decompressor;


	uint64_t ShowTableSaved( uint64_t original_size, uint64_t compacted_size ){
		uint64_t saved = original_size - compacted_size;
		std::cout << std::fixed << std::setprecision(2)
							<< "        original size: " << (original_size>>20) << "MiB, compacted size: "
							<< (compacted_size>>20) << "MiB (" << (compacted_size*100.0/original_size)
							<< "%), saved: " << (saved>>20) << "MiB" << std::endl;
		return saved;
	}

	void CopyData( uint64_t src_position, uint64_t size,  FileDisk *output_file, uint64_t output_position ){
		const uint64_t BUF_SIZE = 128*1024;
		uint8_t buf[BUF_SIZE];
		ProgressCounter progress( size );


		disk_file.Seek( src_position );
		for( uint64_t i = 0; i < size; i += BUF_SIZE ){
			progress.ShowNext( i );

			uint64_t to_copy = std::min( BUF_SIZE, i - size );
			disk_file.Read( buf, to_copy );
			output_file->Write( output_position + i, buf, to_copy );
		}
		progress.ShowNext( size );
	}

	// Compacted table has following structuer:
	// 1. Part with constant size block i.e. line_point, deltas, first_part_of_stubs, pointer_to_overdraft
	// 2. Overdraft - end of tabs that didn't fit in main part (because varies in size of deltas). Writte one after another without spaces deltas.
	// Each entry of 1st part contains of
	// 1. Line Point
	// 2. Deltas
	// 3. Begginig of stubs
	// 4. 2 byes from with restored position and size of overdraft.
	OriginalTableInfo CompactTable( int table_no, FileDisk * output_file, uint64_t output_position, uint16_t min_deltas_sizes ){
		std::cout << " compacting       " << std::flush;

		OriginalTableInfo tinfo( k_size, table_no, table_pointers[table_no-1], table_pointers[table_no] - table_pointers[table_no-1] );
		const uint64_t lps_size =  tinfo.line_point_size + tinfo.stubs_size;
		TableWriter writer( output_file, output_position, min_deltas_sizes, tinfo.parks_count, tinfo.line_point_size, tinfo.stubs_size, tinfo.new_deltas.get() );

		uint8_t buf[tinfo.park_size+7];
		ProgressCounter progress( tinfo.parks_count );

		disk_file.Seek( table_pointers[table_no-1] ); // seek to start of the table
		for( uint64_t i = 0; i < tinfo.parks_count; i++ ){
			progress.ShowNext( i );

			disk_file.Read( buf, tinfo.park_size );
			writer.WriteNext( buf/*line_point+stubs*/, buf+lps_size+2/*deltas_sizes*/, parseDeltasSize( buf+lps_size ) );
		}

		tinfo.new_table_size = writer.getWrittenSize();

		return tinfo;
	}

	// This used instead of CompactTable when "compacting" already compressed or compacted table.
	// but seems it could be joined to one function with CompactTable.
	OriginalTableInfo RealignTable( int table_no, FileDisk * output_file, uint64_t output_position, uint16_t min_deltas_sizes ){
		assert( decompressor );
		std::cout << " realign " << std::flush;

		OriginalTableInfo tinfo( k_size, table_no, table_pointers[table_no-1], table_pointers[table_no] - table_pointers[table_no-1], decompressor.get() );
		TableWriter writer( output_file, output_position, min_deltas_sizes, tinfo.parks_count, tinfo.line_point_size, tinfo.stubs_size, tinfo.new_deltas.get() );

		ProgressCounter progress( tinfo.parks_count );

		for( uint64_t i = 0; i < tinfo.parks_count; i++ ){
			progress.ShowNext( i );

			auto park = decompressor->GetParkReader( disk_file, table_no-1, i*kEntriesPerPark, kEntriesPerPark );
			writer.WriteNext( park.check_point_buf(), park.stubs_buf(), park.deltas_buf(), park.deltas_size );
		}

		tinfo.new_table_size = writer.getWrittenSize();

		return tinfo;
	}

	inline uint16_t parseDeltasSize( uint8_t *buf ){
		uint16_t deltas_size = ( ((uint16_t)buf[0]) | (((uint16_t)buf[1])<<8) );
		if( deltas_size&0x8000 && deltas_size&0x7fff )
			throw std::runtime_error( "Uncompressed deltas does not supported" );
		return deltas_size & 0x7fff;
	}

	uint64_t CompactC3Table( FileDisk * output_file, uint64_t output_position, uint64_t &total_saved, uint16_t min_deltas_sizes ){
		const uint32_t park_size = EntrySizes::CalculateC3Size(k_size);
		const uint64_t parks_count = decompressor ? decompressor->getParksCount(10) : ((table_pointers[10] - table_pointers[9]) / park_size);
		DeltasStorage deltas(parks_count);
		TableWriter writer( output_file, output_position, min_deltas_sizes, parks_count, 0, 0, &deltas );

		uint8_t buf[park_size];

		{ // to show timing at the end of writing
			ProgressCounter progress( parks_count );

			disk_file.Seek( table_pointers[9] );
			for( uint64_t i = 0; i < parks_count; i++ ){
				progress.ShowNext( i );

				if( decompressor )
					writer.WriteNext( buf, decompressor->ReadC3Park( disk_file, i, buf, park_size ) );
				else{
					disk_file.Read( buf, park_size );
					writer.WriteNext( buf + 2, bswap_16( ((uint16_t*)buf)[0] ) );
				}
			}
		}

		// check all partially saved position could be restored
		if( !deltas.IsDeltasPositionRestorable() )
			throw std::runtime_error( "couldn't restore position of deltas" );

		deltas.showStats(); // DEBUG

		uint64_t original_size = table_pointers[10] - table_pointers[9], compacted_size = writer.getWrittenSize();
		total_saved += ShowTableSaved( original_size, compacted_size );

		return parks_count;
	}

	// The compress table is repacked.
	// The small deltas still is deltas but stubs now is not a deltas but real low bits of line point.
	// Than we need to restore all line points and pack them new way
	OriginalTableInfo CompressTable( int table_no, FileDisk * output_file, uint64_t output_position, uint32_t bits_to_remove, uint16_t min_deltas_size ){
		std::cout << " compressing (" << bits_to_remove << ") " << std::flush;

		OriginalTableInfo tinfo( k_size, table_no, table_pointers[table_no-1],
														table_pointers[table_no] - table_pointers[table_no-1], bits_to_remove );
		if( bits_to_remove > tinfo.single_stub_size_bits )
			throw std::runtime_error( "Compression is too high Table 1 stub has " + std::to_string(tinfo.single_stub_size_bits)
															 + " and compression need to remove " + std::to_string(bits_to_remove ) );

		const double R = kRValues[table_no-1];
		const auto max_deltas_size = EntrySizes::CalculateMaxDeltasSize( k_size, table_no );
		TableWriter writer( output_file, output_position, min_deltas_size, tinfo.parks_count, tinfo.new_line_point_size, tinfo.new_stubs_size, tinfo.new_deltas.get() );

		uint8_t read_buf[tinfo.park_size+7],
				new_stubs_buf[ tinfo.new_stubs_size + tinfo.line_point_size + 7/* safe distance*/],
				new_compressed_deltas_buf[kEntriesPerPark];
		ProgressCounter progress( tinfo.parks_count );

		disk_file.Seek( table_pointers[table_no-1] ); // seek to start of the table
		for( uint64_t i = 0; i < tinfo.parks_count; i++ ){
			progress.ShowNext( i );

			disk_file.Read( read_buf, tinfo.park_size ); // read the original park

			ParkReader park( read_buf, k_size*2, tinfo.single_stub_size_bits, max_deltas_size, table_no );


			ParkBits park_stubs_bits;
			std::vector<uint8_t> park_deltas;
			bool deltas_changed = false;
			uint128_t line_point = park.NextLinePoint() >> bits_to_remove;
			Bits line_point_bits( line_point, k_size*2-bits_to_remove );
			line_point_bits.ToBytes( new_stubs_buf ); // save cutted check point line point

			// Restore each line point and repack it new type of stubs
			for( uint32_t j = 0; j < park.DeltasSize(); j++ ){
				assert( park.GetNextIdx() == j+1 );
				uint128_t next_line_point = park.NextLinePoint();
				uint128_t new_line_point = next_line_point>>bits_to_remove;

				uint128_t new_big_delta = new_line_point - line_point;
				uint32_t new_small_delta = new_big_delta >> tinfo.new_single_stub_size_bits;
				if( new_small_delta > 255 )
					throw std::runtime_error( "new delta is too big " + std::to_string(new_small_delta) );
				park_deltas.push_back( new_small_delta );
				if( !deltas_changed && new_small_delta != park.deltas[j] )
					deltas_changed = true;

				uint64_t new_stub = new_big_delta & (tinfo.single_stub_mask >> bits_to_remove);
				park_stubs_bits.AppendValue( new_stub, tinfo.new_single_stub_size_bits );

				line_point = new_line_point;
			}

			// The park repacked now need to save it new way
			size_t deltas_size = park.deltas_size;
			uint8_t * deltas_buf = park.deltas_buf();
			if( deltas_changed ){
				deltas_size = Encoding::ANSEncodeDeltas(park_deltas, R, new_compressed_deltas_buf);
				if( deltas_size <= 0 || deltas_size >= kEntriesPerPark-1 ){
					// uncompressed deltas -> need to think what to do
					throw std::runtime_error( "Unsupported case of uncompressed deltas at park " + std::to_string(i) + " of " + std::to_string( tinfo.parks_count ) );
				}
				deltas_buf = new_compressed_deltas_buf;
			}

			assert( park_stubs_bits.GetSize() ==  park.DeltasSize()*tinfo.new_single_stub_size_bits );

			park_stubs_bits.ToBytes( new_stubs_buf + tinfo.new_line_point_size );
			writer.WriteNext( new_stubs_buf, deltas_buf, deltas_size );
		}

		tinfo.new_table_size = writer.getWrittenSize();

		return tinfo;
	}


	static uint8_t cmp_level_to_bits_cut( uint8_t compression_level ){
		switch( compression_level ){
			case 0: return 0;
			case 1: return 8;
			case 2: return 11;
			case 3: return 14;
			default: return 11 + compression_level;
		}
	}
};

} // end of namespace
#endif  // SRC_CPP_COMPRESSOR_HPP_
