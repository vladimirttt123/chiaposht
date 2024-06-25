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


namespace TCompress {


struct OriginalTableInfo{
public:
	const uint64_t table_pointer, table_size;
	const uint8_t k, line_point_size;
	const uint64_t delta_size, single_stub_size_bits, stubs_size, park_size, parks_count;
	const uint64_t new_single_stub_size_bits, new_stubs_size, new_line_point_size;
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
			, new_single_stub_size_bits( single_stub_size_bits - stubs_bits_to_remove )
			, new_stubs_size( Util::ByteAlign((kEntriesPerPark - 1) * ( new_single_stub_size_bits )) / 8 )
			, new_line_point_size( Util::ByteAlign( 2*k-stubs_bits_to_remove )/8 )
			, single_stub_mask( ((uint64_t)-1)>>(64-single_stub_size_bits))
			, new_deltas( new DeltasStorage( parks_count ) )
	{
	}

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
};

struct BufWriter{
	static const uint32_t SIZE = 1024*1024; // 1Mb
	BufWriter( FileDisk *output_file, uint64_t output_pos )
			: output_file(output_file), write_pos( output_pos ){}

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

	uint64_t GetPosition() const { return write_pos + buf_length; }

	~BufWriter() {Flush();}
private:
	FileDisk *output_file;
	uint64_t write_pos;
	uint32_t buf_length = 0;
	uint8_t buf[SIZE];
};

struct TableWriter {
public:
	const uint64_t output_position, stub_size, parks_count;
	TableWriter( FileDisk * output_file, uint64_t output_position, uint64_t stub_size, uint64_t parks_count )
			:output_position(output_position), stub_size(stub_size), parks_count(parks_count)
			, output_file(output_file), stubs_writer( output_file, output_position )
			, deltas_writer( output_file, output_position + parks_count*stub_size)
			//, stubs_position( output_position ), deltas_position( output_position + parks_count*stub_size )
	{	}

	void WriteNext( uint8_t *stubs, uint8_t * deltas, uint64_t deltas_size ){
		assert( stubs_writer.GetPosition() < output_position + parks_count*stub_size ); // check for too many parks
		if( output_file != nullptr ){
			stubs_writer.Write( stubs, stub_size );
			deltas_writer.Write( deltas, deltas_size );
			// output_file->Write( stubs_position, stubs, stub_size );
			// output_file->Write( deltas_position, deltas, deltas_size );
		}
		// stubs_position += stub_size;
		// deltas_position += deltas_size;
	}

	void WriteNext( uint8_t *stubs, uint8_t *stubs_end, uint8_t * deltas, uint64_t deltas_size ){
		assert( stubs_writer.GetPosition() < output_position + parks_count*stub_size ); // check for too many parks
		if( output_file != nullptr ){
			stubs_writer.Write( stubs, stub_size - 3 );
			stubs_writer.Write( stubs_end, 3 );
			deltas_writer.Write( deltas, deltas_size );

			// output_file->Write( deltas_position, deltas, deltas_size );
			// memcpy( stubs + stub_size-3, stubs_end, 3 );
			// output_file->Write( stubs_position, stubs, stub_size );
		}
		// stubs_position += stub_size;
		// deltas_position += deltas_size;
	}
	uint64_t getWrittenSize() const { return deltas_writer.GetPosition() - output_position; }

	void Flush(){
		stubs_writer.Flush();
		deltas_writer.Flush();
	}
private:
	FileDisk *output_file;
	BufWriter stubs_writer, deltas_writer;
	//uint64_t stubs_position, deltas_position;
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
		// std::cout << "Compress to " << filename << std::endl;
		// std::cout << "Tables pointers: ";
		disk_file.Read( buf, 80 );
		for( int i = 0; i < 10; i++ ){
			table_pointers[i] = Util::EightBytesToInt( buf + i*8 );
			// std::cout << table_pointers[i] << ", ";
		}
		table_pointers[10] = disk_file.Size();

		// std::cout << std::endl;
	}

	void CompressTo( const std::string& filename, uint8_t level ){
		uint8_t bits_to_cut = cmp_level_to_bits_cut( level );
		bool cut_table2 = bits_to_cut > 18;
		if( cut_table2 ) bits_to_cut -= 11;
		std::cout << "Start compressing with level " << (int)level << " (" << (cut_table2?"table2 + ":"") << (int)bits_to_cut << "bits) to " << filename
							<< (filename == "plot.dat" ? "\n *** !!! SEEMS YOU MISSING OUTPUT FILE PARAMETER !!! *** " : "" ) << std::endl;

		if( std::filesystem::exists(filename) )
			throw std::invalid_argument( "output file exists please delete manually: " + filename );

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
		// next 80 bytes is new tables pointers
		header[140+memo_size] = bits_to_cut | (cut_table2?0x80:0); // compression level

		uint64_t *new_table_pointers = (uint64_t*)(header+60+memo_size);
		new_table_pointers[0] = header_size;
		uint32_t *parks_count = (uint32_t*)(header+60+memo_size + 10*8 + 1 );
		uint64_t total_saved = 0;
		Timer total_timer;
		std::vector<std::thread> delta_check_threads;

		FileDisk output_file( filename );
		output_file.Write( 0, header, header_size ); // could be done at the end when everythins is full;
		for( uint32_t i = 1; i < 7; i++ ){
			std::cout << "Table " << i << ": " << std::flush;

			auto tinfo = ((i == 1 && bits_to_cut > 0 ) || (i==2 && cut_table2) ) ?
											CompressTable( i, &output_file, new_table_pointers[i-1], i==1 ? bits_to_cut : 11 )
										: CompactTable( i, &output_file, new_table_pointers[i-1] );
			new_table_pointers[i] = new_table_pointers[i-1] + tinfo.new_table_size;

			if( tinfo.parks_count > 0xffffffffUL ) throw std::runtime_error("too many parks");
			parks_count[i-1] = tinfo.parks_count;

			total_saved += ShowTableSaved( tinfo.table_size, tinfo.new_table_size );

			delta_check_threads.push_back( tinfo.CheckNewDeltas() );
		}

		std::cout << "Tables 7, C1, C2: copy     " << std::flush;
		CopyData( table_pointers[6], table_pointers[9] - table_pointers[6], &output_file, new_table_pointers[6] );

		new_table_pointers[7] = new_table_pointers[6] + (table_pointers[7]-table_pointers[6]);
		new_table_pointers[8] = new_table_pointers[7] + (table_pointers[8]-table_pointers[7]);
		new_table_pointers[9] = new_table_pointers[8] + (table_pointers[9]-table_pointers[8]);

		std::cout << "Table C3: compacting       " << std::flush;
		parks_count[6] = CompactC3Table( &output_file, new_table_pointers[9], total_saved );


		for( uint32_t i = 0; i < delta_check_threads.size(); i++ ) delta_check_threads[i].join();

		// invert table pointers - TODO check endians and skip if not needed
		for( uint32_t i = 0; i < 10; i++ ){
			new_table_pointers[i] = bswap_64( new_table_pointers[i] ); //
			if( i < 7 )
				parks_count[i] = bswap_32( parks_count[i] );
		}

		// save final header
		output_file.Write( 0, header, header_size ); // could be done at the end when everythins is full;
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
	// 1. Part with constant size block i.e. stub, line_point, pointer_to_deltas
	// 2. Deltas - writte one after another without spaces deltas.
	// Each entry of 1st part contains of
	// 1. Line Point
	// 2. Stubs
	// 3. 3 bottom bytes of pointer to deltas location after deltas of this part added.
	//    it is supposed that first pointer is 0 and at length of deltas could be evaluated by substracting prev and next pointer.
	OriginalTableInfo CompactTable( int table_no, FileDisk * output_file, uint64_t output_position ){
		std::cout << " compacting       " << std::flush;

		OriginalTableInfo tinfo( k_size, table_no, table_pointers[table_no-1], table_pointers[table_no] - table_pointers[table_no-1] );
		const uint64_t lps_size =  tinfo.line_point_size + tinfo.stubs_size;
		TableWriter writer( output_file, output_position, 3+lps_size, tinfo.parks_count );

		uint8_t buf[tinfo.park_size+7], small_buf[8];
		ProgressCounter progress( tinfo.parks_count );

		disk_file.Seek( table_pointers[table_no-1] ); // seek to start of the table
		for( uint64_t i = 0; i < tinfo.parks_count; i++ ){
			progress.ShowNext( i );

			disk_file.Read( buf, tinfo.park_size );

			tinfo.new_deltas->Add( i, 0x7fff & ( ((uint32_t)buf[lps_size])
															| (((uint32_t)buf[lps_size + 1])<<8) ) );
			tinfo.new_deltas->TotalEndToBuf( i, small_buf ); // this is buf with pointer already to next postion it is OK

			writer.WriteNext( buf, small_buf, buf+lps_size+2, tinfo.new_deltas->all_sizes[i] );
		}

		tinfo.new_table_size = (lps_size + 3)*tinfo.parks_count + tinfo.new_deltas->total_size;

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
		{ // to show timing at the end of writing
			ProgressCounter progress( parks_count );

			disk_file.Seek( table_pointers[9] );
			for( uint64_t i = 0; i < parks_count; i++ ){
				progress.ShowNext( i );

				disk_file.Read( buf, park_size );
				deltas.Add( i,  Bits(buf, 2, 16).GetValue() );

				output_file->Write( deltas_position, buf + 2, deltas.all_sizes[i] );
				deltas_position += deltas.all_sizes[i];

				deltas.TotalEndToBuf( i, positions + cur_pos_buf );
				cur_pos_buf += 3;
				if( cur_pos_buf >= POS_BUF_SIZE ){
					output_file->Write( next_buf_pos, positions, cur_pos_buf );
					next_buf_pos += POS_BUF_SIZE;
					cur_pos_buf = 0;
				}
			}

			output_file->Write( next_buf_pos, positions, cur_pos_buf ); // write last positions
		}

		// check all partially saved position could be restored
		if( !deltas.IsDeltasPositionRestorable() )
			throw std::runtime_error( "couldn't restore position of deltas" );

		uint64_t original_size = park_size * parks_count,
				compacted_size = (3*parks_count + deltas.total_size);
		total_saved += ShowTableSaved( original_size, compacted_size );

		return parks_count;
	}

	// The compress table is repacked.
	// The small deltas still is deltas but stubs now is not a deltas but real low bits of line point.
	// Than we need to restore all line points and pack them new way
	OriginalTableInfo CompressTable( int table_no, FileDisk * output_file, uint64_t output_position, uint32_t bits_to_remove = 0 ){
		std::cout << " compressing (" << bits_to_remove << ") " << std::flush;

		OriginalTableInfo tinfo( k_size, table_no, table_pointers[table_no-1],
														table_pointers[table_no] - table_pointers[table_no-1], bits_to_remove );
		if( bits_to_remove > tinfo.single_stub_size_bits )
			throw std::runtime_error( "Compression is too high Table 1 stub has " + std::to_string(tinfo.single_stub_size_bits)
															 + " and compression need to remove " + std::to_string(bits_to_remove ) );

		const uint64_t new_lps_size = tinfo.new_line_point_size + tinfo.new_stubs_size;
		const double R = kRValues[table_no-1];

		TableWriter writer( output_file, output_position, 3+new_lps_size, tinfo.parks_count );

		uint8_t buf[tinfo.park_size+7],
				new_stubs_buf[new_lps_size + 3 + 7/* safe distance*/],
				new_compressed_deltas_buf[kEntriesPerPark];
		ProgressCounter progress( tinfo.parks_count );

		disk_file.Seek( table_pointers[table_no-1] ); // seek to start of the table
		for( uint64_t i = 0; i < tinfo.parks_count; i++ ){
			progress.ShowNext( i );

			disk_file.Read( buf, tinfo.park_size ); // read the original park

			ParkReader park( buf, k_size * 2, tinfo.single_stub_size_bits, EntrySizes::CalculateMaxDeltasSize( k_size, table_no ), table_no );

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
					throw std::runtime_error( "uncopressed deltas is not support yet" );
				}
				deltas_buf = new_compressed_deltas_buf;
			}

			assert( park_stubs_bits.GetSize() ==  park.DeltasSize()*tinfo.new_single_stub_size_bits );


			park_stubs_bits.ToBytes( new_stubs_buf + tinfo.new_line_point_size );
			if( park.DeltasSize() < (kEntriesPerPark-1) ){ // need to fill end of stabus by zeros to create same file each time
				uint32_t new_stubs_size_bytes = (park.DeltasSize()*tinfo.new_single_stub_size_bits + 7)/8;
				memset( new_stubs_buf + tinfo.line_point_size + new_stubs_size_bytes, 0, tinfo.new_stubs_size-new_stubs_size_bytes );
			}

			tinfo.new_deltas->Add( i, deltas_size );
			tinfo.new_deltas->TotalEndToBuf( i, new_stubs_buf + new_lps_size );

			writer.WriteNext( new_stubs_buf, deltas_buf, deltas_size );
		}

		tinfo.new_table_size = (new_lps_size + 3)*tinfo.parks_count + tinfo.new_deltas->total_size;
		assert( tinfo.new_table_size == writer.getWrittenSize() );

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
