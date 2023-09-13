#ifndef FILTERED_DISK_HPP_
#define FILTERED_DISK_HPP_

#include <utility>
#include "bitfield.hpp"
#include "disk.hpp"
#include "stream_buffer.hpp"


struct FilteredDisk
{
		FilteredDisk( FileDisk* underlying, MemoryManager &memory_manager, bitfield *filter, int entry_size, uint64_t file_size )
			: buf_size(BUF_SIZE/entry_size*entry_size), buf(buf_size)
			, memory_manager(memory_manager), filter_( filter )
			, underlying_(underlying), entry_size_(entry_size), file_size(file_size)
			, pos_in_buffer( -entry_size )
		{
				assert(entry_size_ > 0);
				buf.ensureSize(buf.size());
		}

		uint8_t const* ReadNext()
		{
				// find next entry
				do {
					pos_in_buffer += entry_size_;
				} while( !filter_->get( ++next_entry_no ) );

				// read up to next entry
				while( buf.used() <= pos_in_buffer ){
					pos_in_buffer -= buf.used();
					buf.setUsed( (uint32_t)std::min( (uint64_t)buf.size(), file_size - pos_in_file ) );
					underlying_->Read( pos_in_file, buf.get(), buf.used() );
					pos_in_file += buf.used();
				}

				return buf.get() + pos_in_buffer;
		}

		void FreeMemory()
		{
				if( filter_ != nullptr ){
					memory_manager.release( filter_->memSize() );
					filter_->RemoveFile();
					delete filter_;
					filter_ = nullptr;

					underlying_->Remove();
				}
		}

private:
		const int buf_size;
		StreamBuffer buf;

		MemoryManager &memory_manager;
		// only entries whose bit is set should be read
		bitfield *filter_;
		FileDisk *underlying_;
		const int entry_size_;
		const uint64_t file_size;

		// position in file where next read from file should be done
		uint64_t pos_in_file = 0;

		// position in buffer where next entry starts
		int64_t pos_in_buffer;

		// no of entry to be read
		int64_t next_entry_no = -1;

};

#endif // FILTERED_DISK_HPP_
