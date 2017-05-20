#ifndef FILEBASEDSPARSEMATRIX_H
#define FILEBASEDSPARSEMATRIX_H
#include <list>
#include <algorithm>
struct SparseRowWrapper
{
    size_t one_cap_;
    size_t one_size_;
    uint32_t one_[];
    static size_t memory_size(size_t cap) 
    {
        return sizeof(size_t) * 2 + sizeof(uint32_t) * cap;
    }
    static size_t capacity(size_t memory_size)
    {
        return (memory_size - 2 * sizeof(size_t)) / sizeof(uint32_t);
    }
    friend std::ostream& operator<<(std::ostream& os, const SparseRowWrapper& srw)
    {
        os << "one_cap_ = " << srw.one_cap_ << ", one_size_ = " << srw.one_size_;
        os << ", one_ = { ";

        for (size_t i = 0; i < srw.one_size_; ++i)\
        {
            os << srw.one_[i];
            if (i < srw.one_size_ - 1)
            os << ", ";
        }
        os << " }";
        return os;
    }
    private:
        SparseRowWrapper();
        SparseRowWrapper(const SparseRowWrapper&);
        SparseRowWrapper& operator=(const SparseRowWrapper&);
};

class FileBasedSparseRowManager;
class FileBasedSparseRow : public ISparseRow
{
    public:
        FileBasedSparseRow(FileBasedSparseRowManager& fbsrm, MemoryMappedFile& mmf, size_t r, off_t offset)
                : fbsrm_(fbsrm), mmf_(mmf), row_(r), offset_(offset), srwp_(0)
        {
            deserialise();
        }
        FileBasedSparseRow(const FileBasedSparseRow& rhs) 
                : fbsrm_(rhs.fbsrm_), mmf_(rhs.mmf_), row_(rhs.row_), offset_(rhs.offset_), srwp_(rhs.srwp_)
        {
        }
        ~FileBasedSparseRow()
        {
            serialise();
        }

        ISparseRow::xor_status xor(size_t col)
        {
            // Note that we must keep columns in one_ sorted.
            // one_ = 2 3 4 9 10
            // one_size_ = 5
            // Case (i) xor(3)
            // Case (ii) xor(6)
            // Case (iii) xor(11)
            // stop search when one_[i] >= col
            const int inc = SparseRow::default_inc;
            //if (col < 0L) return ISparseRow::XOR_FAILED;

            // Optimisation, check for append first, since
            // we're likely to do this often.
            if (srwp_->one_size_ && col > srwp_->one_[srwp_->one_size_ - 1])
            {
                append(col, inc);
                return ISparseRow::XOR_ADDED;
            }

            uint32_t* pos = std::lower_bound(srwp_->one_, srwp_->one_ + srwp_->one_size_, col);
            //size_t* pos = std::lower_bound(one_, one_ + one_size_, col);
            size_t i = static_cast<long int>(pos - srwp_->one_);
            // i is now the position for col in one_
            if (i >= srwp_->one_size_)
            {
                append(col, inc);
                return ISparseRow::XOR_ADDED;
            }

            if (srwp_->one_[i] == col)
            {
                remove(col, i);
                return ISparseRow::XOR_REMOVED;
            }

            insert(col, i, inc);
            return ISparseRow::XOR_ADDED;
        }

        void set_row(const int* first_col, int col_count)
        {
            // ??? NOT YET IMPLEMENTED ???
        }

        void set_row(const std::vector<size_t>& columns)
        {
            // ??? NOT YET IMPLEMENTED ???
        }

        // NOTE: this is an optimisation - only call add(col) if you know that
        // (a) enough space is already allocated for the row
        // (b) the columns are being added in the correct order
        ISparseRow::xor_status add_next(size_t col)
        {
            srwp_->one_[srwp_->one_size_] = col;
            ++srwp_->one_size_;
            return ISparseRow::XOR_ADDED;
        }

        void sort()
        {
            std::sort(srwp_->one_, srwp_->one_ + srwp_->one_size_);
        }

        void clear() 
        {
            srwp_->one_size_ = 0;
        }

        long int highest_column() const
        {
            if (srwp_->one_size_ == 0)
                return -1;
            return static_cast<long int>(srwp_->one_[srwp_->one_size_ - 1]);
        }
        size_t size() const
        {
            return srwp_->one_size_;
        }
        virtual size_t memory_usage()
        {
            return 0;
        }
        typedef SparseRowIterator const_iterator;
        virtual const_iterator begin() const
        {
            return const_iterator(srwp_->one_);
        }
        virtual const_iterator end() const
        {
            return const_iterator(srwp_->one_ + srwp_->one_size_);
        }
        void copy(ISparseRow& copy_of_row)
        {
            for (auto it = begin();
                 it != end();
                 ++it)
            {
                copy_of_row.add_next(*it);
            }
        }
        void compress() {}
    private:
        FileBasedSparseRowManager& fbsrm_;
        MemoryMappedFile& mmf_;
        size_t row_;
        off_t offset_;
        SparseRowWrapper* srwp_;
        void serialise() {}
        void deserialise()
        {
            srwp_ = reinterpret_cast<SparseRowWrapper*>(mmf_.at_offset(offset_));
        }

        void extend(int inc = SparseRow::default_inc);

        void append(size_t col, int inc = SparseRow::default_inc)
        {
            extend(inc);
            srwp_->one_[srwp_->one_size_] = col;
            ++srwp_->one_size_;
        }

        void insert(size_t col, size_t pos, int inc = SparseRow::default_inc)
        {
            extend(inc);
            size_t j = srwp_->one_size_;
            while (j > pos)
            {
                srwp_->one_[j] = srwp_->one_[j - 1];
                --j;
            }
            srwp_->one_[pos] = col;
            ++srwp_->one_size_;
        }

        void remove(size_t col, size_t pos)
        {
            //std::cerr << "remove(" << col << "," << pos << ")" << std::endl;
            while (pos < srwp_->one_size_ - 1)
            {
                srwp_->one_[pos] = srwp_->one_[pos + 1];
                ++pos;
            }
            --srwp_->one_size_;
        }

};

class FileBasedSparseRowManager
{
    private:
        static const size_t initial_row_size = sizeof(size_t) * 2 + sizeof(uint32_t) * 20;
    public:
        FileBasedSparseRowManager(const std::string& filename, size_t row_offset, size_t row_count)
                : filename_(filename), mmf_(filename.c_str(), row_count * initial_row_size), row_offset_(row_offset), row_count_(row_count), row_index_(row_count)
        {
            mmf_.set_default_view_size(1024 * 64);
            off_t offset = 0;
            for (size_t i = 0; i < row_count_; ++i)
            {
                row_index_[i] = offset;
                SparseRowWrapper* srwp = reinterpret_cast<SparseRowWrapper*>(mmf_.at_offset(offset));
                srwp->one_size_ = 0;
                srwp->one_cap_ = 20;
                offset += initial_row_size;
            }
        }
        ~FileBasedSparseRowManager()
        {
            ::unlink(filename_.c_str());
        }

        void set_row(long int r, const int* first_col, int col_count)
        {
            // initialise row to have num_cols columns parsed from s
            if (r < (long int)row_offset_)
            {
                std::ostringstream oss;
                oss << "FileBasedSparseRowManager::row() : row index out of range : " << r;
                throw oss.str();
            }
            size_t index = r - row_offset_;
            if (index >= row_count_)
            {
                extend_file(index);
            }
            extend(r, col_count);
            if (r == 8001489)
            {
                int x = 0;
                ++x;
            }
            FileBasedSparseRow sr(*this, mmf_, r, row_index_[index]);
            for (const int* it = first_col;
                 it != first_col + col_count;
                 ++it)
            {
                sr.add_next(*it);
            }
            sr.sort();
        }

        void set_row(long int r, const std::vector<size_t>& columns)
        {
            // initialise row to have num_cols columns parsed from s
            if (r < (long int)row_offset_)
            {
                std::ostringstream oss;
                oss << "FileBasedSparseRowManager::row() : row index out of range : " << r;
                throw oss.str();
            }
            size_t index = r - row_offset_;
            if (index >= row_count_)
            {
                extend_file(index);
            }
            extend(r, columns.size());
            FileBasedSparseRow sr(*this, mmf_, r, row_index_[index]);
            for (auto it = columns.begin();
                 it != columns.end();
                 ++it)
            {
                sr.add_next(*it);
            }
            sr.sort();
        }

        FileBasedSparseRow row(long int r)
        {
            if (r < (long int)row_offset_)
            {
                std::ostringstream oss;
                oss << "FileBasedSparseRowManager::row() : row index out of range : " << r;
                throw oss.str();
            }
            size_t index = r - row_offset_;
            if (index >= row_count_)
            {
                extend_file(index);
            }
            FileBasedSparseRow sr(*this, mmf_, r, row_index_[index]);
            return sr;
        }

        FileBasedSparseRow row(long int r, long int num_cols, char* s)
        {
            // initialise row to have num_cols columns parsed from s
            if (r < (long int)row_offset_)
            {
                std::ostringstream oss;
                oss << "FileBasedSparseRowManager::row() : row index out of range : " << r;
                throw oss.str();
            }
            size_t index = r - row_offset_;
            if (index >= row_count_)
            {
                extend_file(index);
            }
            extend(r, num_cols);
            FileBasedSparseRow sr(*this, mmf_, r, row_index_[index]);
            while ((s = strtok(0, " ")))
            {
                int n = std::atol(s);
                sr.add_next(n);
            }
            return sr;
        }

        off_t extend(long int r, size_t cap)
        {
            //std::cerr << "FileBasedSparseRowManager::extend, r = " << r << ", cap = " << cap << std::endl;
            size_t index = r - row_offset_;
            if (index >= row_count_)
            {
                // need to create a new row
                // initialise any new rows to default capacity
                extend_file(index);
                off_t offset = row_index_[index];
                SparseRowWrapper* srwp = reinterpret_cast<SparseRowWrapper*>(mmf_.at_offset(offset, SparseRowWrapper::memory_size(cap)));
                srwp->one_cap_ = cap;
                srwp->one_size_ = 0;
                return offset; 
            }

            off_t offset = row_index_[index];
            SparseRowWrapper* old_srwp_ptr = reinterpret_cast<SparseRowWrapper*>(mmf_.at_offset(offset));
            old_srwp_ptr = reinterpret_cast<SparseRowWrapper*>(mmf_.at_offset(offset, SparseRowWrapper::memory_size(old_srwp_ptr->one_cap_)));

            //std::cerr << "FileBasedSparseRowManager::extend, offset = " << offset << ", *old_srwp_ptr = " << *old_srwp_ptr << std::endl;
            std::vector<size_t> row_data(old_srwp_ptr->one_size_);
            for (size_t i = 0; i < old_srwp_ptr->one_size_; ++i)
            {
                row_data[i] = old_srwp_ptr->one_[i];
            }
            add_to_free_list(offset);

            off_t new_offset;
            SparseRowWrapper* new_srwp = find_in_free_list(cap, new_offset);
            if (!new_srwp)
            {
                new_offset = mmf_.size();
                mmf_.set_size(mmf_.size() + SparseRowWrapper::memory_size(cap));
                new_srwp = reinterpret_cast<SparseRowWrapper*>(mmf_.at_offset(new_offset, SparseRowWrapper::memory_size(cap)));
                new_srwp->one_cap_ = cap;
            }
            row_index_[index] = new_offset;
            //std::cerr << "FileBasedSparseRowManager::extend, new_offset = " << new_offset << ", new_srwp = " << *new_srwp << std::endl;
            new_srwp->one_size_ = row_data.size();
            for (size_t i = 0; i < new_srwp->one_size_; ++i)
            {
                new_srwp->one_[i] = row_data[i];
            }
            //std::cerr << "FileBasedSparseRowManager::extend, new_offset = " << new_offset << ", new_srwp = " << *new_srwp << std::endl;
            return new_offset;
        }

        void remove_row(long int r)
        {
            size_t index = r - row_offset_;
            if (index < row_count_)
            {
                off_t offset = row_index_[index];
                row_index_[index] = 0; // ???
                add_to_free_list(offset);
            }
        }

        void replace_row(long int r, long int new_r)
        {
            size_t index = r - row_offset_;
            size_t new_index = new_r - row_offset_;
            if (index < row_count_ && new_index < row_count_)
            {
                off_t offset = row_index_[index];
                row_index_[index] = 0;
                row_index_[new_index] = offset;
            }
        }

    private:
        std::string filename_;
        MemoryMappedFile mmf_;
        size_t row_offset_;
        size_t row_count_;
        std::vector<off_t> row_index_;
        typedef std::list<std::pair<off_t, size_t> > free_list_type;
        free_list_type free_list_;
        void add_to_free_list(off_t offset)
        {
            SparseRowWrapper* srwp_ptr = reinterpret_cast<SparseRowWrapper*>(mmf_.at_offset(offset));
            srwp_ptr->one_size_ = 0;
            size_t s = SparseRowWrapper::memory_size(srwp_ptr->one_cap_);
            if (free_list_.empty())
            {
                free_list_.push_back(std::pair<off_t, size_t>(offset, s));
                return;
            }
            auto it = free_list_.begin();
            free_list_type::iterator prev_it;
            while (it != free_list_.end() && it->first < offset)
            {
                prev_it = it;
                ++it;
            }
            // it is at end of list, or it->first >= offset (since offsets should be unique in free list, this will be it->first > offset)
            if (it == free_list_.end())
            {
                // Can we merge with previous list item?
                if (prev_it->first + prev_it->second == (size_t)offset)
                {
                    prev_it->second += s;
                }
                else
                {
                    free_list_.push_back(std::pair<off_t, size_t>(offset, s));
                }
                return;
            }
            if (it == free_list_.begin())
            {
                // Can we merge with next list item, if any?
                //if (offset + s == (off_t)it->first)
                if ((off_t)s == it->first - offset)
                {
                    it->first = offset;
                    it->second += s;
                }
                else
                {
                    free_list_.insert(it, std::pair<off_t, size_t>(offset, s));
                }
                return;
            }
            // we're somewhere in the middle
            // Can we merge with prev item and/or next item
            //if (prev_it->first + prev_it->second == offset)
            if ((off_t)prev_it->second == offset - prev_it->first)
            {
                prev_it->second += s;
                //if (prev_it->first + prev_it->second == it->first)
                if ((off_t)prev_it->second == it->first - prev_it->first)
                {
                    prev_it->second += it->second;
                    free_list_.erase(it);
                }
            }
            //else if (offset + s == it->first)
            else if ((off_t)s == it->first - offset)
            {
                it->first = offset;
                it->second += s;
            }
            else
            {
                free_list_.insert(it, std::pair<off_t, size_t>(offset, s));
            }
        }

        SparseRowWrapper* find_in_free_list(size_t cap, off_t& offset)
        {
            if (free_list_.empty())
            {
                return 0;
            } 
            auto it = free_list_.begin();
            while (it != free_list_.end() &&
                   SparseRowWrapper::capacity(it->second) < cap)
            {
                ++it;
            }
            if (it == free_list_.end())
            {
                return 0; 
            }
            SparseRowWrapper* srwp = reinterpret_cast<SparseRowWrapper*>(mmf_.at_offset(it->first, it->second));
            srwp->one_cap_ = SparseRowWrapper::capacity(it->second);
            srwp->one_size_ = 0;
            offset = it->first;
            free_list_.erase(it);
            return srwp;
        }

        void extend_from_free_list(size_t index)
        {
            if (free_list_.empty())
            {
                return;
            }
            // 
            // e.g. row_count_ = 1000, index = 1010, need to add 11, free_list_.size() = 5 : to_get <- 5
            //      row_count_ = 1000, index = 1010, need to add 11, free_list_.size() = 10 : to_get <- 10
            //      row_count_ = 1000, index = 1010, need to add 11, free_list_.size() = 20 : to_get <- 11

            size_t to_get = std::min(free_list_.size(), index - row_count_ + 1);
            size_t old_row_count = row_count_;
            row_count_ += to_get;
            row_index_.resize(row_count_);
            auto it = free_list_.begin();
            for (size_t i = 0; i < to_get; ++i, ++it)
            {
                row_index_[old_row_count + i] = it->first;
                SparseRowWrapper* srwp = reinterpret_cast<SparseRowWrapper*>(mmf_.at_offset(it->first, it->second));
                srwp->one_cap_ = SparseRowWrapper::capacity(it->second);
                srwp->one_size_ = 0;
            }
            free_list_.erase(free_list_.begin(), it);
        }

        void extend_file(size_t index)
        {
            //std::cerr << "In FileBasedSparseRowManager::extend_file" << std::endl;
            if (index < row_count_)
                return;
            extend_from_free_list(index);
            if (index < row_count_)
                return;
            size_t old_row_count = row_count_;
            row_count_ = index + 1;
            size_t old_size = mmf_.size();
            mmf_.set_size(mmf_.size() + (index - old_row_count + 1) * initial_row_size);
            row_index_.resize(row_count_);
            //off_t offset = row_index_[old_row_count - 1];
            off_t offset = old_size;
            for (size_t j = old_row_count; j < row_count_; ++j)
            {
                //offset += initial_row_size;
                row_index_[j] = offset;
                SparseRowWrapper* srwp = reinterpret_cast<SparseRowWrapper*>(mmf_.at_offset(offset));
                srwp->one_size_ = 0;
                srwp->one_cap_ = 20;
                offset += initial_row_size;
            }
        }
};

inline void FileBasedSparseRow::extend(int inc)
{
    if (srwp_->one_cap_ <= srwp_->one_size_)
    {
        size_t new_memory_size = SparseRowWrapper::memory_size(srwp_->one_cap_ + inc);
        offset_ = fbsrm_.extend(row_, srwp_->one_cap_ + inc);
        srwp_ = reinterpret_cast<SparseRowWrapper*>(mmf_.at_offset(offset_, new_memory_size));
    }
}
#endif
