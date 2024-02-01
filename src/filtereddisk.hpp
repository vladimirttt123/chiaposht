#ifndef FILTERED_DISK_HPP_
#define FILTERED_DISK_HPP_

#include <utility>
#include <atomic>
#include <thread>

#include "bitfield.hpp"
#include "disk.hpp"
#include "last_table_file.hpp"
#include "util.hpp"

struct FilteredDisk
{
	struct Chunk{
		std::unique_ptr<uint8_t,Util::Deleter<uint8_t>> buf;
		std::atomic_uint_fast8_t state = 0; // 0 - empty, 1 - processing, 2 - ready
		uint64_t file_read_pos, buf_used = 0;
	};

	FilteredDisk( FileDisk* underlying, MemoryManager &memory_manager, bitfield *filter, uint16_t entry_size_bits,
							 uint64_t file_size, uint32_t num_threads, bool parallel_read )
		: entry_size_((entry_size_bits+7)/8), file_size(file_size), num_threads(num_threads)
		, chunk_size( (HUGE_MEM_PAGE_SIZE - MEM_SAFE_BUF_SIZE)/entry_size_/8*entry_size_*8 ), num_chunks(num_threads+1)
		, parallel_read(parallel_read), chunks( new Chunk[num_chunks] ), num_empty_chunks(num_chunks)
		, memory_manager(memory_manager), filter_( filter ), underlying_(underlying)
		, table_reader( *underlying, entry_size_bits )
	{
		assert( entry_size_ > 0 );
		assert( (file_size%entry_size_) == 0 );

		underlying_->setClearAfterRead();
		for( uint64_t i = 0; i < num_chunks; i++ ){
			chunks[i].file_read_pos = i * chunk_size;
			chunks[i].buf = Util::allocate<uint8_t>( chunk_size + MEM_SAFE_BUF_SIZE );
		}
	}

	inline uint8_t EntrySize() const { return entry_size_; }
	inline uint64_t CurrentBucketStart() const { return cur_chunk_start_pos; }
	inline uint64_t CurrentBucketEnd() const { return cur_chunk_start_pos + chunks[cur_chunk_no].buf_used; }
	inline uint64_t CurrentBucketSize() const { return chunks[cur_chunk_no].buf_used; }
	inline bool CurrentBucketIsLast() const { return chunks[(cur_chunk_no+1)%num_chunks].file_read_pos >= file_size; }
	inline const uint8_t * CurrentBucketBuffer() const { return chunks[cur_chunk_no].buf.get(); }
	inline void EnsureSortingStarted() {
		if( running_threads.size() == 0 ){
			RunChunkThread();
			waitForValue( chunks[0].state, (uint8_t)2 ); // wait for first bucket be ready
			assert( chunks[0].state == 2 );
		}
	}

	inline void SwitchNextBucket(){
		if( CurrentBucketIsLast() )
			throw InvalidValueException( "Switching after last bucket" );
		cur_chunk_start_pos += chunks[cur_chunk_no].buf_used;
		chunks[cur_chunk_no].file_read_pos += chunk_size*(uint64_t)num_chunks;
		chunks[cur_chunk_no].state.store( 0, std::memory_order::relaxed ); // free current
		num_empty_chunks.fetch_add( 1, std::memory_order::relaxed );
		num_empty_chunks.notify_all();

		// wait for ready
		cur_chunk_no = (cur_chunk_no+1)%num_chunks; // switch index
		if( waitForValue( chunks[cur_chunk_no].state, (uint8_t)2 ) )
			RunChunkThread();
	}

	inline uint8_t const* ReadNext()
	{
		if( read_next_last_pos >= CurrentBucketEnd() )
			SwitchNextBucket();

		auto res = CurrentBucketBuffer() + read_next_last_pos - CurrentBucketStart();

		read_next_last_pos += entry_size_;

		return res;
	}

	inline void RunChunkThread(){
		if( running_threads.size() < num_threads ){
			running_threads.push_back( new std::thread([this](){ReadChunkThread();}) );
			std::cout << "\rtable1 adding read thread " << running_threads.size() << " of " << num_threads
#ifndef NDEBUG
								<< std::endl
#endif
					;
		}
	}

	inline void ReadChunkThread(){
		auto read_buf = Util::allocate<uint8_t>( chunk_size + MEM_SAFE_BUF_SIZE );
		bitfieldReader filter_reader( *filter_ );
		std::unique_ptr<FileDisk> local_disk;
		std::unique_ptr<TableFileReader> tbl_rd( &table_reader );
		if( parallel_read ){
			local_disk.reset( new FileDisk( underlying_->GetFileName(), false ) );
			local_disk->setClearAfterRead();
			tbl_rd.release();// release previous to not free it
			tbl_rd.reset( new TableFileReader( *local_disk.get(), table_reader.block_size_bytes ) );
		}

		for( bool done = false; !done; ){
			num_empty_chunks.wait( 0, std::memory_order::relaxed );
			for( uint32_t i = 0; i < num_chunks; i++ ){
				if( chunks[i].file_read_pos >= file_size )
					done = true;
				else {
					uint8_t expected_state = 0;
					if( chunks[i].state.compare_exchange_weak( expected_state, 1, std::memory_order::relaxed ) ){
						num_empty_chunks.fetch_sub( 1, std::memory_order_relaxed );
						ReadChunk( chunks[i], read_buf.get(), filter_reader, tbl_rd.get() );
						chunks[i].state.store( 2, std::memory_order::relaxed );
						chunks[i].state.notify_all();
					}
				}
			}
		}
		if( !parallel_read ) tbl_rd.release(); // shouldn't free
	}

	inline void ReadChunk( Chunk & chunk, uint8_t * read_buf, bitfieldReader &filter_reader, TableFileReader *tbl_rd ){
		const uint64_t read_size = std::min( (uint64_t)chunk_size, file_size - chunk.file_read_pos );

		if( filter_->is_readonly() ) filter_read_mutex.lock();// for readonly filter set limits should be done under mutex!
		filter_reader.setLimits( chunk.file_read_pos/entry_size_, read_size/entry_size_ );
		if( filter_->is_readonly() ) filter_read_mutex.unlock();

		if( !parallel_read ) table_read_mutex.lock();
		tbl_rd->ReadBuffer( chunk.file_read_pos, read_buf, read_size );
		if( !parallel_read ) table_read_mutex.unlock();


		if( !filter_->is_readonly() )
			filter_reader.setLimits( chunk.file_read_pos/entry_size_, read_size/entry_size_ );

		chunk.buf_used = 0; // clear the chunk buffer
		for( uint64_t pos_idx = 0, up_to = read_size/entry_size_; pos_idx < up_to; pos_idx++ ){
			if( filter_reader.get(pos_idx) ){
				table_reader.restore_to( read_buf, pos_idx, chunk.buf.get() + chunk.buf_used );
				chunk.buf_used += entry_size_;
			}
		}
	}

	void FreeMemory()
	{
		if( filter_ != nullptr ){
			while( running_threads.size() > 0 ){
				chunks[0].file_read_pos = file_size + 1; // set stopping condition
				num_empty_chunks = num_threads + 2;
				num_empty_chunks.notify_all();// send stopping signal

				running_threads.back()->join();
				delete running_threads.back();
				running_threads.pop_back();
			}

			memory_manager.release( filter_->memSize() );
			filter_->RemoveFile();
			delete filter_;
			filter_ = nullptr;

			underlying_->Remove();

			chunks.reset();
		}
	}

private:
	const uint16_t entry_size_;
	const uint64_t file_size;
	const uint32_t num_threads;
	const uint32_t chunk_size, num_chunks;
	const bool parallel_read;

	std::unique_ptr<Chunk[]> chunks;
	uint64_t cur_chunk_start_pos = 0;
	uint32_t cur_chunk_no = 0;
	std::atomic_int_fast32_t num_empty_chunks;
	std::mutex table_read_mutex, filter_read_mutex;
	std::vector<std::thread*> running_threads;

	uint64_t read_next_last_pos = 0;

	MemoryManager &memory_manager;
	// only entries whose bit is set should be read
	bitfield *filter_;
	FileDisk *underlying_;
	TableFileReader table_reader;
};

#endif // FILTERED_DISK_HPP_

