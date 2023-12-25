// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//    http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef PARK_WRITER_HPP
#define PARK_WRITER_HPP

#include "disk.hpp"
#include "encoding.hpp"
#include "pos_constants.hpp"

// This writes a number of entries into a file, in the final, optimized format. The park
// contains a checkpoint value (which is a 2k bits line point), as well as EPP (entries per
// park) entries. These entries are each divided into stub and delta section. The stub bits are
// encoded as is, but the delta bits are optimized into a variable encoding scheme. Since we
// have many entries in each park, we can approximate how much space each park with take. Format
// is: [2k bits of first_line_point]  [EPP-1 stubs] [Deltas size] [EPP-1 deltas]....
// [first_line_point] ...
void PrepareParkBuffer(
		uint32_t park_size_bytes,
		uint128_t first_line_point,
		const std::vector<uint8_t> &park_deltas,
		const std::vector<uint64_t> &park_stubs,
		uint8_t k,
		uint8_t table_index,
		uint8_t *park_buffer,
		uint64_t const park_buffer_size)
{
	uint8_t *index = park_buffer;

	first_line_point <<= 128 - 2 * k;
	Util::IntTo16Bytes(index, first_line_point);
	index += EntrySizes::CalculateLinePointSize(k);

	// We use ParkBits instead of Bits since it allows storing more data
	ParkBits park_stubs_bits;
	for (uint64_t stub : park_stubs) {
		park_stubs_bits.AppendValue(stub, (k - kStubMinusBits));
	}
	uint32_t stubs_size = EntrySizes::CalculateStubsSize(k);
	uint32_t stubs_valid_size = cdiv(park_stubs_bits.GetSize(), 8);
	park_stubs_bits.ToBytes(index);
	memset(index + stubs_valid_size, 0, stubs_size - stubs_valid_size);
	index += stubs_size;

	// The stubs are random so they don't need encoding. But deltas are more likely to
	// be small, so we can compress them
	double R = kRValues[table_index - 1];
	uint8_t *deltas_start = index + 2;
	size_t deltas_size = Encoding::ANSEncodeDeltas(park_deltas, R, deltas_start);

	if (!deltas_size) {
		// Uncompressed
		deltas_size = park_deltas.size();
		Util::IntToTwoBytesLE(index, deltas_size | 0x8000);
		memcpy(deltas_start, park_deltas.data(), deltas_size);
	} else {
		// Compressed
		Util::IntToTwoBytesLE(index, deltas_size);
	}

	index += 2 + deltas_size;

	if ((uint32_t)(index - park_buffer) > park_buffer_size) {
		std::cout << "index-park_buffer " << index - park_buffer << " park_buffer_size "
							<< park_buffer_size << std::endl;
		throw InvalidStateException(
				"Overflowed park buffer, writing " + std::to_string(index - park_buffer) +
				" bytes. Space: " + std::to_string(park_buffer_size));
	}
	memset(index, 0x00, park_size_bytes - (index - park_buffer));
}

void WriteParkToFile(
		FileDisk &final_disk,
		uint64_t table_start,
		uint64_t park_index,
		uint32_t park_size_bytes,
		uint128_t first_line_point,
		const std::vector<uint8_t> &park_deltas,
		const std::vector<uint64_t> &park_stubs,
		uint8_t k,
		uint8_t table_index,
		uint8_t *park_buffer,
		uint64_t const park_buffer_size)
{
	PrepareParkBuffer( park_size_bytes, first_line_point, park_deltas,
										park_stubs, k, table_index, park_buffer, park_buffer_size );
	// Parks are fixed size, so we know where to start writing. The deltas will not go over
	// into the next park.
	uint64_t writer = table_start + park_index * park_size_bytes;
	final_disk.Write(writer, (uint8_t *)park_buffer, park_size_bytes);
}


struct ParksBuffer{
	uint32_t park_first_index = 0;
	std::unique_ptr<uint8_t,Util::Deleter<uint8_t>> buf;
	std::atomic_uint32_t written_parks = 0;

	ParksBuffer() : buf(Util::allocate<uint8_t>(HUGE_MEM_PAGE_SIZE) ) { }
};

struct ParkAggregator{
	static const uint32_t num_buffers = 2;
	const uint32_t park_size_bytes;
	const uint32_t parks_per_buffer;
	const uint64_t table_start;

	ParkAggregator( FileDisk *final_disk, const uint8_t &k, const uint8_t &table_index, uint64_t table_start, uint32_t park_size_bytes )
			: park_size_bytes(park_size_bytes)
			, parks_per_buffer( (HUGE_MEM_PAGE_SIZE-MEM_SAFE_BUF_SIZE)/park_size_bytes )
			, table_start(table_start)
			, final_disk(final_disk), max_index( parks_per_buffer*num_buffers ){
		for( uint32_t i = 1; i < num_buffers; i++)
			buffers[i].park_first_index = buffers[i-1].park_first_index + parks_per_buffer;
	}

	inline uint8_t * get_for( uint64_t index ){
		uint64_t old;
		while( (old = max_index.load(std::memory_order::relaxed) ) <= index )
			max_index.wait( old, std::memory_order::relaxed );

		uint32_t i = i_for_index( index );
		return buffers[i].buf.get() + (index-buffers[i].park_first_index)*park_size_bytes;
	}

	inline void finish( uint64_t index ){
		uint32_t i = i_for_index( index );
		buffers[i].written_parks.fetch_add(1,std::memory_order::relaxed);
		if( buffers[i].written_parks == parks_per_buffer ){
			{
				std::lock_guard<std::mutex> lk(write_mutex);
				final_disk->Write( table_start + buffers[i].park_first_index * park_size_bytes, buffers[i].buf.get(), parks_per_buffer*park_size_bytes );
			}
			buffers[i].written_parks.store( 0, std::memory_order::relaxed );
			buffers[i].park_first_index = max_index.fetch_add( parks_per_buffer, std::memory_order::relaxed );
			max_index.notify_all();
		}
	}
	~ParkAggregator(){
		bool chk = true;
		for( uint32_t i = 0; i < num_buffers; i++ )
			if( buffers[i].written_parks > 0 ){
				assert( chk ); chk = false;
				final_disk->Write( table_start + buffers[i].park_first_index * park_size_bytes, buffers[i].buf.get(), parks_per_buffer*park_size_bytes );
			}
	}
private:
	FileDisk *final_disk;
	ParksBuffer buffers[num_buffers];
	std::atomic_uint64_t max_index;
	std::mutex write_mutex;

	inline uint32_t i_for_index( uint64_t index ){
		for( uint32_t i = 0; i < num_buffers; i++ ){
			if( buffers[i].park_first_index <= index && index < (buffers[i].park_first_index + parks_per_buffer) )
				return i;
		}

		throw InvalidStateException( "Park aggregator for processed index" );
	}
};

// Thread safe version of park writer
struct ParkWriterTS{

	ParkWriterTS( ParkAggregator&aggregator, const uint8_t &k, const uint8_t &table_index )
			: aggregator(aggregator), k(k), table_index(table_index)	{}

	void Write( uint64_t park_index, uint128_t first_line_point,
						 std::unique_ptr<std::vector<uint8_t>> &park_deltas,
						 std::unique_ptr<std::vector<uint64_t>> &park_stubs ){

		// for( uint32_t k = 15; k < 55; k++ ){
		// 	for( uint32_t table_index = 1; table_index < 8; table_index ++ ){
		// 		auto was = EntrySizes::CalculateLinePointSize(k)
		// 										 + EntrySizes::CalculateStubsSize(k) + 2
		// 											 + EntrySizes::CalculateMaxDeltasSize(k, 1);
		// 		auto now = EntrySizes::CalculateParkSize(k, table_index );
		// 		std::cout << "k=" << k << "; t=" << table_index << "; was=" << was << "; now=" << now << "; was-now=" << (was - now) << std::endl;
		// 		//assert( was == now );
		// 	}
		// }

		PrepareParkBuffer( aggregator.park_size_bytes,
											first_line_point, *park_deltas.get(), *park_stubs.get(),
											k, table_index,
											aggregator.get_for( park_index ), aggregator.park_size_bytes );

		aggregator.finish( park_index );

		park_deltas->clear();
		park_stubs->clear();
	}

private:
	ParkAggregator &aggregator;
	const uint8_t k;
	const uint8_t table_index;
};


#endif // PARK_WRITER_HPP
