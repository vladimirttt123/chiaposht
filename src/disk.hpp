// Copyright 2018 Chia Network Inc

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//    http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef SRC_CPP_DISK_HPP_
#define SRC_CPP_DISK_HPP_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <thread>
#include <mutex>


// enables disk I/O logging to disk.log
// use tools/disk.gnuplot to generate a plot
#define ENABLE_LOGGING 0

using namespace std::chrono_literals; // for operator""min;

#include "chia_filesystem.hpp"

#include "./bits.hpp"
#include "./util.hpp"
#include "threading.hpp"

constexpr uint64_t write_cache = 256 * 1024;
constexpr uint64_t read_ahead = 256 * 1024;
uint64_t BUF_SIZE = 256*1024;
bool LEAVE_FILES = false;


struct Disk {
    virtual uint8_t const* Read(uint64_t begin, uint64_t length) = 0;
    virtual void Write(uint64_t begin, const uint8_t *memcache, uint64_t length) = 0;
    virtual void Truncate(uint64_t new_size) = 0;
    virtual std::string GetFileName() = 0;
    virtual void FreeMemory() = 0;
    virtual ~Disk() = default;
};

#if ENABLE_LOGGING
// logging is currently unix / bsd only: use <fstream> or update
// calls to ::open and ::write to port to windows
#include <fcntl.h>
#include <unistd.h>
#include <mutex>
#include <unordered_map>
#include <cinttypes>

enum class op_t : int { read, write};

void disk_log(fs::path const& filename, op_t const op, uint64_t offset, uint64_t length)
{
    static std::mutex m;
    static std::unordered_map<std::string, int> file_index;
    static auto const start_time = std::chrono::steady_clock::now();
    static int next_file = 0;

    auto const timestamp = std::chrono::steady_clock::now() - start_time;

    int fd = ::open("disk.log", O_WRONLY | O_CREAT | O_APPEND, 0755);

    std::unique_lock<std::mutex> l(m);

    char buffer[512];

    int const index = [&] {
        auto it = file_index.find(filename.string());
        if (it != file_index.end()) return it->second;
        file_index[filename.string()] = next_file;

        int const len = std::snprintf(buffer, sizeof(buffer)
            , "# %d %s\n", next_file, filename.string().c_str());
        ::write(fd, buffer, len);
        return next_file++;
    }();

    // timestamp (ms), start-offset, end-offset, operation (0 = read, 1 = write), file_index
    int const len = std::snprintf(buffer, sizeof(buffer)
        , "%" PRId64 "\t%" PRIu64 "\t%" PRIu64 "\t%d\t%d\n"
        , std::chrono::duration_cast<std::chrono::milliseconds>(timestamp).count()
        , offset
        , offset + length
        , int(op)
        , index);
    ::write(fd, buffer, len);
    ::close(fd);
}
#endif

struct FileDisk {
		explicit FileDisk( const fs::path &filename, bool forWrite = true )
    {
        filename_ = filename;

				if( forWrite ) {
					Open( writeFlag );
					SetCouldBeClosed();
				}
    }

    void Open(uint8_t flags = 0)
    {
			UnsetCouldBeClosed( );
        // if the file is already open, don't do anything
				if (f_)  return;

        // Opens the file for reading and writing
        do {
#ifdef _WIN32
            f_ = ::_wfopen(filename_.c_str(), (flags & writeFlag) ? L"w+b" : L"r+b");
#else
            f_ = ::fopen(filename_.c_str(), (flags & writeFlag) ? "w+b" : "r+b");
#endif
            if (f_ == nullptr) {
                std::string error_message =
										"Could not open " + filename_.string() + ": err" + std::to_string(errno) + " " + ::strerror(errno) + ".";

								if( errno == 24 && could_be_closed.size() > 0 ){
									std::cout << "Warning: Too many open file, forced to close some. Closing " << could_be_closed.size() << " files." << std::endl;
									CloseCouldBeClosed();
								}
								else if (errno == 24 || (flags & retryOpenFlag) ) {
                    std::cout << error_message << " Retrying in five minutes." << std::endl;
										// Close as much as possible to allow remount
										CloseCouldBeClosed();
										std::this_thread::sleep_for(5min);
                } else {
                    throw InvalidValueException(error_message);
                }
            }
        } while (f_ == nullptr);
    }

    FileDisk(FileDisk &&fd)
    {
        filename_ = std::move(fd.filename_);
        f_ = fd.f_;
        fd.f_ = nullptr;
				writeMax = fd.writeMax;
    }

    FileDisk(const FileDisk &) = delete;
    FileDisk &operator=(const FileDisk &) = delete;

    void Close()
    {
        if (f_ == nullptr) return;
				UnsetCouldBeClosed();
        ::fclose(f_);
        f_ = nullptr;
        readPos = 0;
        writePos = 0;
				{
					std::lock_guard<std::mutex> lk(mutFileManager);
					total_bytes_written += bytes_written;
				}
				bytes_written = 0;
    }

		void Remove( bool noWarn = false ){
			Close();
			if( LEAVE_FILES )
				RenameFileToDeleted();
			else{
				if( !fs::remove( GetFileName() ) && !noWarn )
					std::cout << "Warning: Some problem with file removing: " << GetFileName() << std::endl;
			}
		}

		~FileDisk() {
			UnsetCouldBeClosed();
			Close();
		}

    void Read(uint64_t begin, uint8_t *memcache, uint64_t length)
    {
        Open(retryOpenFlag);
#if ENABLE_LOGGING
        disk_log(filename_, op_t::read, begin, length);
#endif
        // Seek, read, and replace into memcache
        uint64_t amtread;
        do {
						if( !bReading ) Flush();
            if ((!bReading) || (begin != readPos)) {
#ifdef _WIN32
								_fseeki64(f_, begin, SEEK_SET);
#else
                // fseek() takes a long as offset, make sure it's wide enough
                static_assert(sizeof(long) >= sizeof(begin));
								auto seek_res = ::fseek( f_, begin, SEEK_SET );
								if( seek_res != 0 ){
									std::cout << "Error of seeking to position " << begin << " where writeMax= " << writeMax
														<< " in file " << filename_ << std::endl;
									throw InvalidStateException( "Cann't seek to " + std::to_string(begin) + " in file " + GetFileName() );
								}
#endif
                bReading = true;
            }
            amtread = ::fread(reinterpret_cast<char *>(memcache), sizeof(uint8_t), length, f_);
            readPos = begin + amtread;
            if (amtread != length) {
                std::cout << "Only read " << amtread << " of " << length << " bytes at offset "
                          << begin << " from " << filename_ << " with length " << writeMax
													<< ". Error " << ferror(f_) << ": " << ::strerror(ferror(f_))
													<< ". Retrying in five minutes." << std::endl;
                // Close, Reopen, and re-seek the file to recover in case the filesystem
                // has been remounted.
                Close();
								CloseCouldBeClosed();
                bReading = false;
								std::this_thread::sleep_for(5min);
                Open(retryOpenFlag);
            }
        } while (amtread != length);
			SetCouldBeClosed();
    }

    void Write(uint64_t begin, const uint8_t *memcache, uint64_t length)
    {
        Open(writeFlag | retryOpenFlag);
#if ENABLE_LOGGING
        disk_log(filename_, op_t::write, begin, length);
#endif
        // Seek and write from memcache
        uint64_t amtwritten;
        do {
            if ((bReading) || (begin != writePos)) {
#ifdef _WIN32
                _fseeki64(f_, begin, SEEK_SET);
#else
                // fseek() takes a long as offset, make sure it's wide enough
                static_assert(sizeof(long) >= sizeof(begin));
                ::fseek(f_, begin, SEEK_SET);
#endif
                bReading = false;
            }
            amtwritten =
                ::fwrite(reinterpret_cast<const char *>(memcache), sizeof(uint8_t), length, f_);
            writePos = begin + amtwritten;
            if (writePos > writeMax)
                writeMax = writePos;
            if (amtwritten != length) {
                // If an error occurs, the resulting value of the file-position indicator for the stream is unspecified.
                // https://pubs.opengroup.org/onlinepubs/007904975/functions/fwrite.html
                //
                // And in the code above if error occurs with 0 bytes written (full disk) it will not reset the pointer
                // (writePos will still be equal to begin), however it need to be reseted.
                //
                // Otherwise this causes #234 - in phase3, when this bucket is read, it goes into endless loop.
                //
                // Thanks tinodj!
                writePos = UINT64_MAX;
                std::cout << "Only wrote " << amtwritten << " of " << length << " bytes at offset "
                          << begin << " to " << filename_ << " with length " << writeMax
                          << ". Error " << ferror(f_) << ". Retrying in five minutes." << std::endl;
                // Close, Reopen, and re-seek the file to recover in case the filesystem
                // has been remounted.
                Close();
                bReading = false;
                std::this_thread::sleep_for(5min);
                Open(writeFlag | retryOpenFlag);
            }
        } while (amtwritten != length);
			bytes_written += amtwritten;
			SetCouldBeClosed();
    }

		std::string GetFileName() const { return filename_.string(); }

    uint64_t GetWriteMax() const noexcept { return writeMax; }

    void Truncate(uint64_t new_size)
    {
			Close();
			if( LEAVE_FILES && new_size == 0 )
				RenameFileToDeleted();
			else
				fs::resize_file(filename_, new_size);
    }

		void Flush() {
			if( f_)
				if( ::fflush(f_) != 0 )
					std::cout << "Fail: Cannot flush file " << filename_ << std::endl;
		}

		static uint64_t GetTotalBytesWritten() { return total_bytes_written; }
private:

    uint64_t readPos = 0;
    uint64_t writePos = 0;
    uint64_t writeMax = 0;
		uint64_t bytes_written = 0;
    bool bReading = true;

    fs::path filename_;
    FILE *f_ = nullptr;

    static const uint8_t writeFlag = 0b01;
    static const uint8_t retryOpenFlag = 0b10;

		void RenameFileToDeleted() const {
			if( !fs::exists(GetFileName() )) return;
			std::string delname = GetFileName() + ".deleted";
			for( int i = 0; fs::exists(delname); i++ )
				delname = GetFileName() + "." + std::to_string(i) + ".deleted";
			fs::rename( GetFileName(), delname );
		}

		static std::vector<FileDisk*> could_be_closed;
		static std::mutex mutFileManager;
		static uint64_t total_bytes_written;

		inline void SetCouldBeClosed(){
			std::lock_guard<std::mutex> lk(mutFileManager);
			could_be_closed.push_back( this );
		}

		inline void UnsetCouldBeClosed(){
			std::lock_guard<std::mutex> lk(mutFileManager);
			for( uint32_t i = 0; i < could_be_closed.size(); i++ ){
				if( could_be_closed[i] == this ){
					could_be_closed[i] = could_be_closed[could_be_closed.size() - 1];
					could_be_closed.resize( could_be_closed.size() - 1 );
				}
			}
		}

		static void CloseCouldBeClosed() {
			std::lock_guard<std::mutex> lk(mutFileManager);
			for( auto f : could_be_closed )
				f->Close();
			could_be_closed.clear();
		}
};

std::vector<FileDisk*> FileDisk::could_be_closed = std::vector<FileDisk*>();
std::mutex FileDisk::mutFileManager = std::mutex();
uint64_t FileDisk::total_bytes_written = 0;

uint64_t GetTotalBytesWritten(){ return FileDisk::GetTotalBytesWritten(); }

struct BufferedDisk : Disk
{
    BufferedDisk(FileDisk* disk, uint64_t file_size) : disk_(disk), file_size_(file_size) {}

    uint8_t const* Read(uint64_t begin, uint64_t length) override
    {
        assert(length < read_ahead);
				if( begin > file_size_ )
					throw InvalidValueException( "Begin of read position is byond the end of the file: begin = "
									+ std::to_string(begin) + "; file_size = " + std::to_string(file_size_)
									+ "; filename = " + GetFileName() );
				if( (begin+length ) > file_size_ )
					throw InvalidValueException( "Try to read after end of the file: begin = "
									+ std::to_string(begin) + "; length=" + std::to_string( begin + length )
									+ "; file_size = " + std::to_string(file_size_) + "; filename = " + GetFileName() );

        NeedReadCache();
        // all allocations need 7 bytes head-room, since
        // SliceInt64FromBytes() may overrun by 7 bytes
        if (read_buffer_start_ <= begin
            && read_buffer_start_ + read_buffer_size_ >= begin + length
            && read_buffer_start_ + read_ahead >= begin + length + 7)
        {
            // if the read is entirely inside the buffer, just return it
            return read_buffer_.get() + (begin - read_buffer_start_);
        }
        else if (begin >= read_buffer_start_ || begin == 0 || read_buffer_start_ == std::uint64_t(-1)) {

            // if the read is beyond the current buffer (i.e.
            // forward-sequential) move the buffer forward and read the next
            // buffer-capacity number of bytes.
            // this is also the case we enter the first time we perform a read,
            // where we haven't read anything into the buffer yet. Note that
            // begin == 0 won't reliably detect that case, sinec we may have
            // discarded the first entry and start at some low offset but still
            // greater than 0
            read_buffer_start_ = begin;
            uint64_t const amount_to_read = std::min(file_size_ - read_buffer_start_, read_ahead);
            disk_->Read(begin, read_buffer_.get(), amount_to_read);
            read_buffer_size_ = amount_to_read;
            return read_buffer_.get();
        }
        else {
            // ideally this won't happen
            std::cout << "Disk read position regressed. It's optimized for forward scans. Performance may suffer\n"
                << "   read-offset: " << begin
                << " read-length: " << length
                << " file-size: " << file_size_
                << " read-buffer: [" << read_buffer_start_ << ", " << read_buffer_size_ << "]"
                << " file: " << disk_->GetFileName()
                << '\n';
            static uint8_t temp[128];
            // all allocations need 7 bytes head-room, since
            // SliceInt64FromBytes() may overrun by 7 bytes
            assert(length <= sizeof(temp) - 7);

            // if we're going backwards, don't wipe out the cache. We assume
            // forward sequential access
            disk_->Read(begin, temp, length);
            return temp;
        }
    }

    void Write(uint64_t const begin, const uint8_t *memcache, uint64_t const length) override
    {
        NeedWriteCache();
        if (begin == write_buffer_start_ + write_buffer_size_) {
            if (write_buffer_size_ + length <= write_cache) {
                ::memcpy(write_buffer_.get() + write_buffer_size_, memcache, length);
                write_buffer_size_ += length;
                return;
            }
            FlushCache();
        }

        if (write_buffer_size_ == 0 && write_cache >= length) {
            write_buffer_start_ = begin;
            ::memcpy(write_buffer_.get() + write_buffer_size_, memcache, length);
            write_buffer_size_ = length;
            return;
        }

        disk_->Write(begin, memcache, length);
    }

    void Truncate(uint64_t const new_size) override
    {
        FlushCache();
        disk_->Truncate(new_size);
        file_size_ = new_size;
        FreeMemory();
    }

    std::string GetFileName() override { return disk_->GetFileName(); }

		void FreeMemory( ) override {
			FreeMemory( true );
		}

		void FreeMemory( bool withFlashCache )
    {
				if( withFlashCache )
					FlushCache();

        read_buffer_.reset();
        write_buffer_.reset();
        read_buffer_size_ = 0;
        write_buffer_size_ = 0;
    }

    void FlushCache()
    {
				if (write_buffer_size_ > 0) {
					disk_->Write(write_buffer_start_, write_buffer_.get(), write_buffer_size_);
					write_buffer_size_ = 0;
				}
				disk_->Flush();
    }

private:

    void NeedReadCache()
    {
        if (read_buffer_) return;
        read_buffer_.reset(new uint8_t[read_ahead]);
        read_buffer_start_ = -1;
        read_buffer_size_ = 0;
    }

    void NeedWriteCache()
    {
        if (write_buffer_) return;
        write_buffer_.reset(new uint8_t[write_cache]);
        write_buffer_start_ = -1;
        write_buffer_size_ = 0;
    }

    FileDisk* disk_;

    uint64_t file_size_;

    // the file offset the read buffer was read from
    uint64_t read_buffer_start_ = -1;
    std::unique_ptr<uint8_t[]> read_buffer_;
    uint64_t read_buffer_size_ = 0;

    // the file offset the write buffer should be written back to
    // the write buffer is *only* for contiguous and sequential writes
    uint64_t write_buffer_start_ = -1;
    std::unique_ptr<uint8_t[]> write_buffer_;
    uint64_t write_buffer_size_ = 0;
};




struct BufferedReader{
	BufferedReader( FileDisk* disk, uint64_t const input_disk_begin, uint64_t buffer_size, uint64_t bytes_to_read )
		: disk_(disk), buffer_size_(std::min(buffer_size,bytes_to_read)), file_read_position_(input_disk_begin), bytes_to_read_(bytes_to_read)
	{
		if( bytes_to_read_ ){
			buffer = new uint8_t[buffer_size_<<1]; // double size buffer creation
			read_thread = new std::thread( [this]{this->ReadBuffer();} );
		}
	}

	inline uint64_t MoveNextBuffer(){
		if( read_thread == nullptr ){
			assert( bytes_read_from_file_ == bytes_to_read_ );
			return current_buffer_size_ = 0;
		}
		read_thread->join(); // wait for read
		delete read_thread;
		read_thread = nullptr;

		current_buffer_start_byte_ = next_buffer_start_byte_;
		current_buffer_size_ = bytes_read_from_file_ - current_buffer_start_byte_;

		assert( current_buffer_size_ > 0 );
		assert( (current_buffer_start_byte_ + current_buffer_size_) <= bytes_to_read_ );

		current_buffer_ ^= 1; // switch to next buffer
		if( bytes_read_from_file_ < bytes_to_read_ )
			read_thread = new std::thread( [this]{this->ReadBuffer();} );

		return current_buffer_size_;
	}

	inline const uint8_t * GetBuffer() const { return buffer + buffer_size_ * current_buffer_;	}

	inline uint64_t GetBufferStartPosition() const { return file_read_position_ + current_buffer_start_byte_; }
	inline uint64_t BufferSize() const { return current_buffer_size_; }
	inline bool isAtEnd() const { return GetBufferStartPosition() + buffer_size_ >= bytes_to_read_; }


	~BufferedReader(){
		if( read_thread != nullptr ){
			read_thread->join();
			delete read_thread;
		} else {
			assert( bytes_read_from_file_ == bytes_to_read_ );
		}
		if( buffer != nullptr ) delete [] buffer;
	}


private:
	FileDisk *disk_;
	const uint64_t buffer_size_;
	uint8_t *buffer = nullptr;
	const uint64_t file_read_position_;
	const uint64_t bytes_to_read_;

	uint64_t current_buffer_ = 1;
	uint64_t current_buffer_size_ = 0;
	uint64_t next_buffer_start_byte_ = 0;
	uint64_t current_buffer_start_byte_ = 0;
	uint64_t bytes_read_from_file_ = 0;

	std::thread *read_thread = nullptr;


	void ReadBuffer(){
		int64_t bytes_to_read = std::min( buffer_size_, bytes_to_read_ - bytes_read_from_file_ );
		assert( bytes_to_read > 0 );

		disk_->Read( file_read_position_ + bytes_read_from_file_, buffer + buffer_size_ * (current_buffer_^1), bytes_to_read );

		next_buffer_start_byte_ = bytes_read_from_file_;
		bytes_read_from_file_ += bytes_to_read;
	}

};

#endif  // SRC_CPP_DISK_HPP_
