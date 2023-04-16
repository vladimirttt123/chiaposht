#ifndef FILTERED_DISK_HPP_
#define FILTERED_DISK_HPP_

#include <utility>
#include "bitfield.hpp"
#include "disk.hpp"


struct FilteredDisk : Disk
{
		FilteredDisk(BufferedDisk underlying, MemoryManager &memory_manager, bitfield *filter, int entry_size)
			: memory_manager(memory_manager), filter_( filter )
			, underlying_(std::move(underlying)), entry_size_(entry_size)
		{
				assert(entry_size_ > 0);
				while (!filter_->get(last_idx_)) {
						last_physical_ += entry_size_;
						++last_idx_;
				}
				assert(filter_->get(last_idx_));
				assert(last_physical_ == last_idx_ * entry_size_);
		}

		uint8_t const* Read(uint64_t begin, uint64_t length) override
		{
				// we only support a single read-pass with no going backwards
				assert(begin >= last_logical_);
				assert((begin % entry_size_) == 0);
				assert(filter_->get(last_idx_));
				assert(last_physical_ == last_idx_ * entry_size_);

				if (begin > last_logical_) {
						// last_idx_ et.al. always points to an entry we have (i.e. the bit
						// is set). So when we advance from there, we always take at least
						// one step on all counters.
						last_logical_ += entry_size_;
						last_physical_ += entry_size_;
						++last_idx_;

						while (begin > last_logical_)
						{
								if (filter_->get(last_idx_)) {
										last_logical_ += entry_size_;
								}
								last_physical_ += entry_size_;
								++last_idx_;
						}

						while (!filter_->get(last_idx_)) {
								last_physical_ += entry_size_;
								++last_idx_;
						}
				}

				assert(filter_->get(last_idx_));
				assert(last_physical_ == last_idx_ * entry_size_);
				assert(begin == last_logical_);
				return underlying_.Read(last_physical_, length);
		}

		void Write(uint64_t begin, const uint8_t *memcache, uint64_t length) override
		{
				assert(false);
				throw std::runtime_error("Write() called on read-only disk abstraction");
		}
		void Truncate(uint64_t new_size) override
		{
				underlying_.Truncate(new_size);
				if (new_size == 0) FreeMemory();
		}
		std::string GetFileName() override { return underlying_.GetFileName(); }

		void FreeMemory() override
		{
				if( filter_ != nullptr ){
					memory_manager.release( filter_->memSize() );
					filter_->FreeMemory();
					delete filter_;
					filter_ = nullptr;
				}
				underlying_.FreeMemory();
		}


private:
		MemoryManager &memory_manager;
		// only entries whose bit is set should be read
		bitfield *filter_;
		BufferedDisk underlying_;
		int entry_size_;

		// the "physical" disk offset of the last read
		uint64_t last_physical_ = 0;
		// the "logical" disk offset of the last read. i.e. the offset as if the
		// file would have been compacted based on filter_
		uint64_t last_logical_ = 0;

		// the index of the last read. This is also the index into the bitfield. It
		// could be computed as last_physical_ / entry_size_, but we want to avoid
		// the division.
		uint64_t last_idx_ = 0;
};

#endif // FILTERED_DISK_HPP_
