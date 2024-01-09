#ifndef FILTERED_DISK_HPP_
#define FILTERED_DISK_HPP_

#include <utility>
#include <atomic>
#include <thread>

#include "bitfield.hpp"
#include "disk.hpp"
#include "stream_buffer.hpp"
#include "util.hpp"

struct FilteredDisk
{
		FilteredDisk( FileDisk* underlying, MemoryManager &memory_manager, bitfield *filter, int entry_size, uint64_t file_size, uint32_t num_threads )
			: bucket_size( (HUGE_MEM_PAGE_SIZE - MEM_SAFE_BUF_SIZE)/2/entry_size*entry_size*std::max(1U,num_threads)), entry_size_(entry_size), file_size(file_size)
			, buckets_buf( Util::allocate<uint8_t>( bucket_size*2 + MEM_SAFE_BUF_SIZE ) ), bucket_buf_start( bucket_size )
			, last_pos_in_buffer(-entry_size), memory_manager(memory_manager), filter_( filter ), underlying_(underlying)
		{
			assert( entry_size_ > 0 );
			assert( (file_size%entry_size) == 0 );
		}

		inline uint8_t EntrySize() const { return entry_size_; }
		inline uint64_t CurrentBucketStart() const { return bucket_end_pos - bucket_buf_used; }
		inline uint64_t CurrentBucketEnd() const { return bucket_end_pos; }
		inline uint64_t CurrentBucketSize() const { return bucket_buf_used; }
		inline bool CurrentBucketIsLast() const { return is_bucket_last; }
		inline const uint8_t * CurrentBucketBuffer() const { return buckets_buf.get() + bucket_buf_start; }
		inline void EnsureSortingStarted() {
			if( read_thread ) return;
			read_thread.reset( new std::thread( [this](){
				StreamBuffer file_read_buf( BUF_SIZE/entry_size_*entry_size_ );
				// position in buffer where next entry starts
				int64_t pos_in_buffer = -entry_size_;
				// no of entry to be read
				int64_t next_entry_no = -1;

				uint64_t file_read_pos = 0, basket_buf_start_pos = 0;

				while( !stop_read_thread ) {

					next_bucket_ready.wait( true, std::memory_order::relaxed );

					uint64_t i = 0;
					for( ; i < bucket_size; i+=entry_size_ ){
						// find next entry
						do {
							pos_in_buffer += entry_size_;

							if( file_read_buf.used() <= pos_in_buffer ) {// read next buffer
								pos_in_buffer -= file_read_buf.used();

								file_read_buf.setUsed( std::min( (uint64_t)file_read_buf.size(), this->file_size - file_read_pos ) );
								if( file_read_buf.used() == 0 ){
									is_next_bucket_last = true;
									next_bucket_size = i;
									next_bucket_ready.store( true, std::memory_order::relaxed );
									next_bucket_ready.notify_all();
									return;
								}

								underlying_->Read( file_read_pos, file_read_buf.get(), file_read_buf.used() );
								file_read_pos += file_read_buf.used();
							} // reading from file
						} while( !filter_->get( ++next_entry_no ) );


						memcpy( buckets_buf.get() + basket_buf_start_pos + i, file_read_buf.get() + pos_in_buffer, entry_size_ );
					} // end of loop filling global buffer

					basket_buf_start_pos ^= bucket_size; // switch to next global buffer for next read
					is_next_bucket_last = pos_in_buffer >= file_read_buf.used() && this->file_size == file_read_pos;
					next_bucket_size = i;

					next_bucket_ready.store( true, std::memory_order::relaxed );
					next_bucket_ready.notify_all();
				}
			} ) );

			SwitchNextBucket();
		}

		inline void SwitchNextBucket(){
			assert( read_thread );
			next_bucket_ready.wait( false, std::memory_order::relaxed );
			bucket_buf_start ^= bucket_size;
			is_bucket_last = is_next_bucket_last;
			bucket_buf_used = next_bucket_size;
			bucket_end_pos += next_bucket_size;
			next_bucket_ready.store( false, std::memory_order::relaxed );
			next_bucket_ready.notify_all();
		}

		inline uint8_t const* ReadNext()
		{
			assert( !read_thread );
			assert( last_pos_in_buffer < bucket_buf_used || position_in_file < file_size );

			do {
				last_pos_in_buffer += entry_size_;
			}while( !filter_->get( ++last_entry_idx ) );


			while( last_pos_in_buffer >= bucket_buf_used ) {// read next buffer
				last_pos_in_buffer -= bucket_buf_used;

				bucket_buf_used = std::min( (uint64_t)bucket_size, file_size - position_in_file );
				assert( bucket_buf_used > 0 );
				underlying_->Read( position_in_file, buckets_buf.get(), bucket_buf_used );
				position_in_file += bucket_buf_used;
			} // reading from file

			assert( last_pos_in_buffer < bucket_buf_used );
			return buckets_buf.get() + last_pos_in_buffer;
		}

		void FreeMemory()
		{
			if( filter_ != nullptr ){
				if( read_thread ){
					stop_read_thread = true;
					read_thread->join();
					read_thread.reset();
				}

				memory_manager.release( filter_->memSize() );
				filter_->RemoveFile();
				delete filter_;
				filter_ = nullptr;

				underlying_->Remove();

				buckets_buf.reset();
			}
		}

private:
	const uint32_t bucket_size;
	const uint16_t entry_size_;
	const uint64_t file_size;

	std::unique_ptr<uint8_t, Util::Deleter<uint8_t>> buckets_buf;
	uint32_t bucket_buf_used = 0, bucket_buf_start;
	std::atomic_bool next_bucket_ready = false;
	uint64_t bucket_end_pos = 0, position_in_file = 0, next_bucket_size = 0;
	uint64_t last_pos_in_buffer, last_entry_idx = -1;

	std::unique_ptr<std::thread>read_thread;
	bool is_bucket_last = false, is_next_bucket_last = false;
	bool stop_read_thread = false;

	MemoryManager &memory_manager;
	// only entries whose bit is set should be read
	bitfield *filter_;
	FileDisk *underlying_;
};

