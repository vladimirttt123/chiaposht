#ifndef FILTERED_DISK_HPP_
#define FILTERED_DISK_HPP_

#include <utility>
#include <atomic>
#include <thread>

#include "bitfield.hpp"
#include "disk.hpp"
#include "util.hpp"

struct FilteredDisk
{
		FilteredDisk( FileDisk* underlying, MemoryManager &memory_manager, bitfield *filter, int entry_size, uint64_t file_size )
			: buf_size( (HUGE_MEM_PAGE_SIZE - MEM_SAFE_BUF_SIZE)/2/entry_size*entry_size), entry_size_(entry_size), file_size(file_size)
			, buf( Util::allocate<uint8_t>( buf_size*2 + MEM_SAFE_BUF_SIZE ) ), buf_start( buf_size )
			, memory_manager(memory_manager), filter_( filter ), underlying_(underlying)
			, pos_in_buffer( -entry_size )
		{
			assert( entry_size_ > 0 );
			assert( (file_size%entry_size) == 0 );
			read_thread = new std::thread( [this](){
				for( uint64_t read_pos = 0, next_read_size, start_pos = 0;
						 !stop_read_thread && read_pos < this->file_size;
						 read_pos += next_read_size, start_pos ^= buf_size ){

					while( next_buf_used.load(std::memory_order::relaxed) > 0 )
						next_buf_used.wait( next_read_size, std::memory_order::relaxed ); // whait for free buffer

					// read next buffer
					next_read_size = std::min( (uint64_t)buf_size, this->file_size - read_pos );
					underlying_->Read( read_pos, buf.get() + start_pos, next_read_size );
					next_buf_used.store( (uint32_t)next_read_size, std::memory_order::relaxed );
					next_buf_used.notify_all();
				}
			} );
		}

		uint8_t const* ReadNext()
		{
				// find next entry
				do {
					pos_in_buffer += entry_size_;
				} while( !filter_->get( ++next_entry_no ) );

				// read up to next entry
				while( buf_used <= pos_in_buffer ){
					pos_in_buffer -= buf_used;
					while( next_buf_used.load(std::memory_order::relaxed) == 0 )
						next_buf_used.wait( 0, std::memory_order::relaxed ); // whait for next buffer

					buf_used = next_buf_used.load(std::memory_order::relaxed);
					buf_start ^= buf_size;
					next_buf_used.store( 0, std::memory_order::relaxed );
					next_buf_used.notify_all();
				}

				return buf.get() + buf_start + pos_in_buffer;
		}

		void FreeMemory()
		{
			if( filter_ != nullptr ){
				stop_read_thread = true;
				next_buf_used = 0; // stopping read thread

				memory_manager.release( filter_->memSize() );
				filter_->RemoveFile();
				delete filter_;
				filter_ = nullptr;


				underlying_->Remove();

				read_thread->join();
				delete read_thread;

				buf.reset();
			}
		}

private:
	const uint32_t buf_size;
	const uint16_t entry_size_;
	const uint64_t file_size;

	std::unique_ptr<uint8_t, Util::Deleter<uint8_t>> buf;
	uint32_t buf_used = 0, buf_start;
	std::atomic_uint32_t next_buf_used = 0;

	std::thread *read_thread;
	bool stop_read_thread = false;

	MemoryManager &memory_manager;
	// only entries whose bit is set should be read
	bitfield *filter_;
	FileDisk *underlying_;

	// position in file where next read from file should be done
	//uint64_t pos_in_file = 0;

	// position in buffer where next entry starts
	int64_t pos_in_buffer;

	// no of entry to be read
	int64_t next_entry_no = -1;
};

#endif // FILTERED_DISK_HPP_
