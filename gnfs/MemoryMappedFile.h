#ifndef MEMORYMAPPEDFILE_H
#define MEMORYMAPPEDFILE_H
#include <string>
#include <iostream>
#include <sstream>
#ifdef WIN32
#include <windows.h>
#else
#if defined(linux) || defined(sun)
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#else
#include <w32api/windows.h>
#undef max
#undef min
#endif
#endif
#if defined(linux) || defined(sun)
class MemoryMappedFile
{
public:
    MemoryMappedFile(const char* filename)
        : file_handle_(0), base_(0), pos_(0), allocation_granularity_(0), view_offset_(0), file_size_(0),
          view_size_(0), default_view_size_(DEFAULT_VIEW_SIZE), writeable_(false)
    {
        struct stat statbuf;

        file_handle_ = open(filename, O_RDONLY);

        if (file_handle_ == -1)
        {
            // std::cerr << "Failed to open file [" << filename << "], errno = <" << errno << "> : " << strerror(errno) << std::endl;
            throw std::string("MemoryMappedFile: Failed to open file");
        }

        fstat(file_handle_, &statbuf);
        file_size_ = statbuf.st_size;
        view_offset_ = 0;
        view_size_ = default_view_size_;
        allocation_granularity_ = sysconf(_SC_PAGE_SIZE);
        map_file();
    }
    MemoryMappedFile(const char* filename, unsigned long long file_size)
        : file_handle_(0), base_(0), pos_(0), allocation_granularity_(sysconf(_SC_PAGE_SIZE)), view_offset_(0),
          file_size_(file_size), view_size_(DEFAULT_VIEW_SIZE), default_view_size_(DEFAULT_VIEW_SIZE), writeable_(true)
    {
        // create a new file of given size
        file_handle_ = open(filename, O_CREAT|O_TRUNC|O_RDWR, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
        lseek(file_handle_, file_size_ - 1, SEEK_SET);
        write(file_handle_, "", 1);
        map_file();
    }
    ~MemoryMappedFile()
    {
        unmap_file();
        if (file_handle_ != -1)
        {
            close(file_handle_);
        }
    }
    const char* pos() const
    {
        return pos_;
    }
    char* at_offset(unsigned long long int o)
    {
        return at_offset(o, default_view_size_);
    }
    char* at_offset(unsigned long long int o, size_t len)
    {
        if (o < view_offset_ || o + len > view_offset_ + view_size_)
        {
            view_offset_ = o;
            size_t excess = static_cast<size_t>(view_offset_ % allocation_granularity_);
            view_offset_ -= excess;
            unmap_file();
            view_size_ = default_view_size_;
            map_file();
            pos_ += excess;
        }
        return base_ + (o - view_offset_);
    }
    void reset()
    {
        unmap_file();
        view_offset_ = 0;
        view_size_ = default_view_size_;
        map_file();
    }
    unsigned long long int size() const
    {
        return file_size_;
    }
    void set_size(unsigned long long int new_size)
    {
        if (!writeable_)
        {
            return;
        }
        if (new_size > file_size_)
        {
            unmap_file();
            file_size_ = new_size;
            lseek(file_handle_, file_size_ - 1, SEEK_SET);
            write(file_handle_, "", 1);
            map_file();
        }
        else if (new_size < file_size_)
        {
            unmap_file();
            file_size_ = new_size;
            ftruncate(file_handle_, file_size_);
            view_offset_ = 0;
            map_file();
        }
    }
    void set_default_view_size(size_t default_view_size)
    {
        default_view_size_ = default_view_size;
    }
    size_t get_default_view_size() const
    {
        return default_view_size_;
    }
    friend bool getline(MemoryMappedFile& mmfile, std::string& str);
private:
    void advance(size_t s)
    {
        pos_ += s;
    }
    void move_view(const char* hint)
    {
        view_offset_ += hint - base_;
        size_t excess = static_cast<size_t>(view_offset_ % allocation_granularity_);
        view_offset_ -= excess;
        unmap_file();
        view_size_ = default_view_size_;
        map_file();
        pos_ += excess;
    }
    void map_file()
    {
        if (view_offset_ > size())
        {
            std::cerr << "WARNING : MemoryMappedFile::map_file() : trying to map off end of file" << std::endl;
            view_size_ = default_view_size_;
            if (size() > view_size_)
            {
                view_offset_ = size() - view_size_;
            }
            else
            {
                view_offset_ = 0;
                view_size_ = size();
            }
        }
        if (view_offset_ + view_size_ > size())
        {
            //std::cerr << "MemoryMappedFile::map_file() : view_offset_ = " << view_offset_ << ", view_size_ = " << view_size_ << ", size() = " << size() << std::endl;
            view_size_ = static_cast<size_t>(size() - view_offset_);
            //std::cerr << "MemoryMappedFile::map_file() : view_offset_ = " << view_offset_ << ", view_size_ = " << view_size_ << ", size() = " << size() << std::endl;
        }

        int prot = PROT_READ;
        if (writeable_)
        {
            prot = prot | PROT_WRITE;
        }

        base_ = static_cast<char*>(mmap(0, view_size_, prot, MAP_SHARED, file_handle_, view_offset_));
        if (base_ == MAP_FAILED)
        {
            std::cerr << "mmap failed: errno = <" << errno << "> : " << strerror(errno) << " : view_size_ = " << view_size_ << ", view_offset_ = " << view_offset_ << std::endl;
            throw std::string("mmap_failed");
        }
        else
        {
            pos_ = base_;
        }
    }
    void unmap_file()
    {
        munmap(base_, view_size_);
    }
    const char* view_end() const
    {
        return base_ + view_size_;
    }
    const char* end() const
    {
        off_t end_offset = size() - view_offset_;
        if (end_offset >> 31)
        {
            return reinterpret_cast<char*>(~0);
        }
        else
        {
            return base_ + end_offset;
        }
    }

    int file_handle_;
    char* base_;
    char* pos_;
    size_t allocation_granularity_;
    //off_t view_offset_;
    unsigned long long int view_offset_;
    //off_t file_size_;
    unsigned long long int file_size_;
    size_t view_size_;
    size_t default_view_size_;
    bool writeable_;
    enum { DEFAULT_VIEW_SIZE = 32 * 1024 * 1024};
};
#else
#ifdef WIN32
namespace
{
void report_error(const std::string& msg, const std::string& extra_msg = std::string(""))
{
    DWORD error_code = GetLastError();
    LPVOID lpMsgBuf;
    FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER |
                  FORMAT_MESSAGE_FROM_SYSTEM |
                  FORMAT_MESSAGE_IGNORE_INSERTS,
                  NULL,
                  error_code,
                  MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                  (LPTSTR) &lpMsgBuf,
                  0,
                  NULL);
    std::cerr << msg << " : [" << static_cast<char*>(lpMsgBuf) << "]" << std::endl;
    if (!extra_msg.empty())
    {
        std::cerr << extra_msg << std::endl;
    }
    throw std::string(static_cast<char*>(lpMsgBuf));
}
}
class MemoryMappedFile
{
public:
    MemoryMappedFile(const char* filename, unsigned long long file_size = 0)
        : file_handle_(INVALID_HANDLE_VALUE), mapping_handle_(INVALID_HANDLE_VALUE),
          base_(0), pos_(0), allocation_granularity_(0), view_offset_(0),
          file_size_low_(0), file_size_high_(0), view_size_(DEFAULT_VIEW_SIZE), default_view_size_(DEFAULT_VIEW_SIZE), writeable_(file_size != 0)
    {
        DWORD desired_access = GENERIC_READ;
        DWORD share_mode = 0;
        DWORD creation_disposition = OPEN_ALWAYS;
        if (writeable_)
        {
            // need to create new file in read/write mode
            desired_access = GENERIC_READ|GENERIC_WRITE;
            share_mode = FILE_SHARE_READ|FILE_SHARE_WRITE;
            creation_disposition = CREATE_ALWAYS;
        }
        file_handle_ = CreateFile(filename, desired_access, share_mode, NULL, creation_disposition, FILE_ATTRIBUTE_NORMAL, NULL);

        if (file_handle_ == INVALID_HANDLE_VALUE)
        {
            report_error("CreateFile failed");
        }
        set_file_size_and_map(file_size);
    }
    ~MemoryMappedFile()
    {
        if (base_)
        {
            unmap_file();
        }
        if (mapping_handle_ != INVALID_HANDLE_VALUE)
        {
            CloseHandle(mapping_handle_);
        }
        if (file_handle_ != INVALID_HANDLE_VALUE)
        {
            CloseHandle(file_handle_);
        }
    }
    const char* pos() const
    {
        return pos_;
    }
    char* at_offset(unsigned long long int o)
    {
        return at_offset(o, default_view_size_);
    }
    char* at_offset(unsigned long long int o, size_t len)
    {
        if (o < view_offset_ || o + len > view_offset_ + view_size_)
        {
            view_offset_ = o;
            size_t excess = static_cast<size_t>(view_offset_ % allocation_granularity_);
            view_offset_ -= excess;
            unmap_file();
            view_size_ = default_view_size_;
            while (o + len > view_offset_ + view_size_)
            {
                view_size_ += default_view_size_;
            }
            map_file();
            pos_ += excess;
        }
        return base_ + (o - view_offset_);
    }
    void reset()
    {
        unmap_file();
        view_offset_ = 0;
        view_size_ = default_view_size_;
        map_file();
    }
    unsigned long long int size() const
    {
        unsigned long long int big_size = file_size_high_;
        big_size <<= 32;
        big_size += file_size_low_;
        return big_size;
    }
    void set_size(unsigned long long int new_size)
    {
        if (!writeable_)
        {
            return;
        }
        unmap_file();
        set_file_size_and_map(new_size);
    }
    void set_default_view_size(size_t default_view_size)
    {
        default_view_size_ = default_view_size;
    }
    size_t get_default_view_size() const
    {
        return default_view_size_;
    }
    friend bool getline(MemoryMappedFile& mmfile, std::string& str);
private:
    void set_file_size_and_map(unsigned long long file_size)
    {
        DWORD protect = PAGE_READONLY;
        if (writeable_)
        {
            // need to seek to end of file, given by file_size, to create file of right size
            LARGE_INTEGER file_size_li;
            file_size_li.QuadPart = file_size;
            SetFilePointerEx(file_handle_, file_size_li, NULL, FILE_BEGIN);
            SetEndOfFile(file_handle_);
            protect = PAGE_READWRITE;
        }

        file_size_low_ = GetFileSize(file_handle_, &file_size_high_);
        if (mapping_handle_ != INVALID_HANDLE_VALUE)
        {
            CloseHandle(mapping_handle_);
        }
        mapping_handle_ = CreateFileMapping(file_handle_, NULL, protect, 0, 0, NULL);
        if (!mapping_handle_)
        {
            report_error("CreateFileMapping failed");
        }

        SYSTEM_INFO system_info;
        GetSystemInfo(&system_info);
        allocation_granularity_ = system_info.dwAllocationGranularity;
        view_offset_ = 0;
        view_size_ = default_view_size_;
        map_file();
    }
    void advance(size_t s)
    {
        pos_ += s;
    }
    void move_view(const char* hint)
    {
        unmap_file();
        view_offset_ += hint - base_;
        size_t excess = static_cast<size_t>(view_offset_ % allocation_granularity_);
        view_offset_ -= excess;
        view_size_ = default_view_size_;
        map_file();
        pos_ += excess;
    }
    void map_file()
    {
        if (view_offset_ + view_size_ > size())
        {
            view_size_ = static_cast<size_t>(size() - view_offset_);
        }

        DWORD offset_high = static_cast<DWORD>(view_offset_ >> 32);
        DWORD offset_low = static_cast<DWORD>(view_offset_ & 0xFFFFFFFF);
        DWORD prot = FILE_MAP_READ;
        if (writeable_)
        {
            prot |= FILE_MAP_WRITE;
        }
        base_ = static_cast<char*>(MapViewOfFile(mapping_handle_, prot, offset_high, offset_low, view_size_));
        if (!base_)
        {
            std::ostringstream oss;
            oss << "offset_high = " << offset_high << ", offset_low = " << offset_low << ", view_size_ = " << static_cast<unsigned int>(view_size_);
            report_error("MapViewOfFile failed", oss.str());
        }
        else
        {
            pos_ = base_;
        }
    }
    void unmap_file()
    {
        UnmapViewOfFile(base_);
    }
    const char* view_end() const
    {
        return base_ + view_size_;
    }
    const char* end() const
    {
        unsigned long long int end_offset = size() - view_offset_;
        if (end_offset >> 32)
        {
            return reinterpret_cast<char*>(~0);
        }
        else
        {
            return base_ + end_offset;
        }
    }

    HANDLE file_handle_;
    HANDLE mapping_handle_;
    char* base_;
    char* pos_;
    DWORD allocation_granularity_;
    DWORD file_size_low_;
    DWORD file_size_high_;
    unsigned long long int view_offset_;
    size_t view_size_;
    size_t default_view_size_;
    bool writeable_;
    enum { DEFAULT_VIEW_SIZE = 32 * 1024 * 1024};
};
#endif
#endif

inline bool getline(MemoryMappedFile& mmfile, std::string& str)
{
    const char* p = mmfile.pos();
    const char* q = p;
    const char* e = mmfile.end();
    const char* ve = mmfile.view_end();
    while (q != e && q < ve && *q != 0x0d && *q != 0x0a) ++q;
    if (q == e)
        return false;
    if (q >= ve) mmfile.move_view(p);
    p = mmfile.pos();
    q = p;
    e = mmfile.end();
    ve = mmfile.view_end();
    while (q != e && q != ve && *q != 0x0d && *q != 0x0a) ++q;
    if (q == e)
        return false;
    size_t s = q - p + 1;
    if (q != e && *q == 0x0d) ++s;
    str.assign(p, q);
    mmfile.advance(s);
    return true;
}
#endif
