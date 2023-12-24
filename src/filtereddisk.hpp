#ifndef FILTERED_DISK_HPP_
#define FILTERED_DISK_HPP_

#include <utility>
#include <atomic>
#include <thread>

#include "bitfield.hpp"
#include "disk.hpp"
#include "stream_buffer.hpp"
#include "util.hpp"

#define NEW_FILTERED
#ifdef NEW_FILTERED
struct FilteredDisk
{
		FilteredDisk( FileDisk* underlying, MemoryManager &memory_manager, bitfield *filter, int entry_size, uint64_t file_size )
			: bucket_size( (HUGE_MEM_PAGE_SIZE - MEM_SAFE_BUF_SIZE)/2/entry_size*entry_size), entry_size_(entry_size), file_size(file_size)
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
#else //NEW_FILTERED
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
#endif // NEW_FILTERED
#endif // FILTERED_DISK_HPP_
