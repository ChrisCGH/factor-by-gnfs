#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H
#include "MemoryMappedFile.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <memory.h>
#include <stdint.h>
#ifndef WIN32
#include <unistd.h>
#include <uuid/uuid.h>
#endif
#include "BitOperations.h"

// class to implement a sparse matrix over F_2

class SparseRowIterator
{
public:
    SparseRowIterator(uint32_t* ptr) : ptr_(ptr)
    {}
    SparseRowIterator(const SparseRowIterator& sri) : ptr_(sri.ptr_)
    {}
    SparseRowIterator& operator++();
    size_t operator*() const;

    bool operator==(const SparseRowIterator& sri) const;
    bool operator!=(const SparseRowIterator& sri) const
    {
        return !(*this == sri);
    }
private:
    uint32_t* ptr_;
    //size_t* ptr_;
};

class ISparseRow
{
public:
    enum xor_status { XOR_FAILED = 0, XOR_ADDED, XOR_REMOVED };
    virtual xor_status xor(size_t col) = 0;
    virtual xor_status add_next(size_t col) = 0;
    virtual void set_row(const std::vector<size_t>& columns) = 0;
    virtual void set_row(const int* first_col, int col_count) = 0;
    virtual void clear() = 0;
    virtual long int highest_column() const = 0;
    virtual size_t size() const = 0;
    virtual size_t memory_usage() = 0;
    typedef SparseRowIterator const_iterator;
    virtual const_iterator begin() const = 0;
    virtual const_iterator end() const = 0;
    virtual void compress() = 0;
    virtual ~ISparseRow() {};
};

class SparseMatrix;
class SparseRow : public ISparseRow
{
public:
    friend class SparseMatrix;
    friend class SparseMatrix2;
    friend class SparseMatrix3;
    friend class SparseMatrix4;
    friend void transpose(const SparseMatrix& A, SparseMatrix& A_t, long int max_row_size, bool clear_A);
    enum { default_inc = 50 };
    SparseRow() : one_(0), one_size_(0), one_cap_(0)
    {}
    SparseRow(long int cols) : one_(0), one_size_(0), one_cap_(cols)
    {
        one_ = new uint32_t [ one_cap_ ];
        //one_ = new size_t [ one_cap_ ];
    }
private:

    SparseRow(const std::string& str) : one_(0), one_size_(0), one_cap_(0)
    {
        long int num_cols = 0;
        char* s = SparseRow::begin_parse(str, num_cols);
        if (!s) return;
        one_ = new uint32_t [ num_cols ];
        //one_ = new size_t [ num_cols ];
        one_size_ = num_cols;
        one_cap_ = num_cols;
        int i = 0;
        while ((s = strtok(0, " ")))
        {
            int n = std::atol(s);
            one_[i] = n;
            i++;
        }
    }

    SparseRow(long int num_cols, char* s) : one_(0), one_size_(0), one_cap_(0)
    {
        one_ = new uint32_t [ num_cols ];
        //one_ = new size_t [ num_cols ];
        one_size_ = num_cols;
        one_cap_ = num_cols;
        int i = 0;
        while ((s = strtok(0, " ")))
        {
            int n = std::atol(s);
            one_[i] = n;
            i++;
        }
    }
public:
    virtual ~SparseRow()
    {
        delete [] one_;
    }
private:
    SparseRow& operator=(const SparseRow& sr)
    {
        if (sr.one_size_ > one_cap_)
        {
            delete [] one_;
            one_ = new uint32_t [ sr.one_size_ ];
            //one_ = new size_t [ sr.one_size_ ];
            if (!one_)
            {
                std::cerr << "operator new failed!" << std::endl;
                std::exit(0);
            }
            one_cap_ = sr.one_size_;
        }
        one_size_ = sr.one_size_;
        memcpy(one_, sr.one_, one_size_ * sizeof(uint32_t));
        //memcpy(one_, sr.one_, one_size_ * sizeof(long int));
        return *this;
    }
    friend std::ostream& operator<<(std::ostream& os, const SparseRow& r);
    long int highest_column() const
    {
        if (one_size_ == 0) return -1;
        return static_cast<long int>(one_[one_size_ - 1]);
    }
    void clear()
    {
        one_size_ = 0;
    }

    static char* begin_parse(const std::string& str, long int& column_count)
    {
        // This function gets the column count from str using strtok(buf, " "),
        // and returns a char* buf suitable for subsequent calls to strtok(0, " ")
        static char* buf = 0;
        static std::string::size_type buflen = 0;
        column_count = 0;
        if (str.empty()) return 0;
        if (str.size() > buflen)
        {
            delete [] buf;
            buflen = str.size();
            buf = new char [ buflen + 1 ];
        }
        strcpy(buf, str.c_str());
        char* s = strtok(buf, " ");
        column_count = std::atol(s);
        return buf;
    }

    void compress();
private:
    SparseRow(const SparseRow& sr);
    bool operator==(const SparseRow& sr) const;
    bool operator!=(const SparseRow& sr) const;

    void extend(int inc = default_inc)
    {
        try
        {
            if (one_cap_ <= one_size_)
            {
                one_cap_ += inc;
                uint32_t* new_one_ = new uint32_t [ one_cap_ ];
                //size_t* new_one_ = new size_t [ one_cap_ ];
                memcpy(new_one_, one_, one_size_ * sizeof(uint32_t));
                //memcpy(new_one_, one_, one_size_ * sizeof(size_t));
                delete [] one_;
                one_ = new_one_;
            }
        }
        catch (std::bad_alloc& e)
        {
            std::ostringstream oss;
            oss << "SparseRow::extend() : Bad allocation - tried to allocate " << static_cast<unsigned int>(one_cap_) << " long ints : " << e.what();
            throw oss.str();
        }
    }

    void append(size_t col, int inc = default_inc)
    {
        extend(inc);
        one_[one_size_] = col;
        ++one_size_;
    }

    void insert(size_t col, size_t pos, int inc = default_inc)
    {
        extend(inc);
        size_t j = one_size_;
        while (j > pos)
        {
            one_[j] = one_[j - 1];
            --j;
        }
        one_[pos] = col;
        ++one_size_;
    }

    void remove(size_t col, size_t pos)
    {
        //std::cerr << "remove(" << col << "," << pos << ")" << std::endl;
        while (pos < one_size_ - 1)
        {
            one_[pos] = one_[pos + 1];
            ++pos;
        }
        --one_size_;
    }

    bool add(size_t col)
    {
        if (one_size_ > 0 && col <= one_[one_size_ - 1]) return false;
        append(col);
        return true;
    }

    //public:

    void set_row(const int* first_col, int col_count)
    {
        if (col_count > (int)one_cap_)
        {
            one_cap_ = col_count;
            delete [] one_;
            one_ = new uint32_t [ one_cap_ ];
        }
        one_size_ = 0;
        for (const int* it = first_col;
                it != first_col + col_count;
                ++it)
        {
            one_[one_size_] = *it;
            ++one_size_;
        }
        std::sort(one_, one_ + one_size_);
    }

    void set_row(const std::vector<size_t>& columns)
    {
        if (columns.size() > one_cap_)
        {
            one_cap_ = columns.size();
            delete [] one_;
            one_ = new uint32_t [ one_cap_ ];
        }
        one_size_ = 0;
        for (auto& col: columns)
        {
            one_[one_size_] = col;
            ++one_size_;
        }
        std::sort(one_, one_ + one_size_);
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
        const int inc = default_inc;
        //if (col < 0L) return ISparseRow::XOR_FAILED;

        // Optimisation, check for append first, since
        // we're likely to do this often.
        if (one_size_ && col > one_[one_size_ - 1])
        {
            append(col, inc);
            return ISparseRow::XOR_ADDED;
        }

        uint32_t* pos = std::lower_bound(one_, one_ + one_size_, col);
        //size_t* pos = std::lower_bound(one_, one_ + one_size_, col);
        size_t i = static_cast<long int>(pos - one_);
        // i is now the position for col in one_
        if (i >= one_size_)
        {
            append(col, inc);
            return ISparseRow::XOR_ADDED;
        }

        if (one_[i] == col)
        {
            remove(col, i);
            return ISparseRow::XOR_REMOVED;
        }

        insert(col, i, inc);
        return ISparseRow::XOR_ADDED;
    }
    // NOTE: this is an optimisation - only call add(col) if you know that
    // (a) enough space is already allocated for the row
    // (b) the columns are being added in the correct order
    ISparseRow::xor_status add_next(size_t col)
    {
        one_[one_size_] = col;
        ++one_size_;
        return ISparseRow::XOR_ADDED;
    }
public:
    const_iterator begin() const
    {
        return SparseRowIterator(one_);
    }
    const_iterator end() const
    {
        return SparseRowIterator(one_ + one_size_);
    }
private:
    size_t size() const
    {
        return one_size_;
    }
    size_t memory_usage()
    {
        size_t s = sizeof(*this);
        s += one_cap_ * sizeof(uint32_t);
        //s += one_cap_ * sizeof(long int);
        return s;
    }
private:
    uint32_t* one_;
    size_t one_size_;
    size_t one_cap_;
};

#include "FileBasedSparseMatrix.h"

#define FAST_BM 1
#ifdef FAST_BM
struct FastVector
{
    FastVector() : vec_(0), size_(0), capacity_(0)
    {}
    FastVector(const FastVector& fv) : vec_(0), size_(0), capacity_(0)
    {
        resize(fv.size_);
        memcpy(vec_, fv.vec_, size_ * sizeof(uint32_t));
    }
    FastVector& operator=(const FastVector& fv)
    {
        if (this == &fv) return *this;
        resize(fv.size_);
        memcpy(vec_, fv.vec_, size_ * sizeof(uint32_t));
        return *this;
    }
    ~FastVector()
    {
        delete [] vec_;
    }

    uint32_t* vec_;
    size_t size_;
    size_t capacity_;
    void clear()
    {
        size_ = 0;
    }
    size_t size() const
    {
        return size_;
    }
    void resize(size_t s, bool keep = false)
    {
        const size_t incr = 100;
        if (s > capacity_)
        {
            capacity_ = s + incr;
            uint32_t* tmp = new uint32_t [ capacity_ ];
            memset(tmp, 0, capacity_ * sizeof(uint32_t));
            if (keep)
            {
                memcpy(tmp, vec_, size_ * sizeof(uint32_t));
            }
            delete [] vec_;
            vec_ = tmp;
        }
        for (size_t i = size_; i < s; i++) vec_[i] = 0UL;
        size_ = s;
    }
    uint32_t* begin() const
    {
        return vec_;
    }
    uint32_t* end() const
    {
        return vec_ + size_;
    }
    uint32_t& operator[] (size_t i)
    {
        return vec_[i];
    }
    const uint32_t& operator[] (size_t i) const
    {
        return vec_[i];
    }
};
struct FastVector64
{
    FastVector64() : vec_(0), size_(0), capacity_(0)
    {}
    FastVector64(const FastVector64& fv) : vec_(0), size_(0), capacity_(0)
    {
        resize(fv.size_);
        memcpy(vec_, fv.vec_, size_ * sizeof(unsigned long long int));
    }
    FastVector64& operator=(const FastVector64& fv)
    {
        if (this == &fv) return *this;
        resize(fv.size_);
        memcpy(vec_, fv.vec_, size_ * sizeof(unsigned long long int));
        return *this;
    }
    ~FastVector64()
    {
        delete [] vec_;
    }

    unsigned long long int* vec_;
    size_t size_;
    size_t capacity_;
    void clear()
    {
        size_ = 0;
    }
    size_t size() const
    {
        return size_;
    }
    void resize(size_t s, bool keep = false)
    {
        const size_t incr = 100;
        if (s > capacity_)
        {
            capacity_ = s + incr;
            unsigned long long int* tmp = new unsigned long long int [ capacity_ ];
            memset(tmp, 0, capacity_ * sizeof(unsigned long long int));
            if (keep)
            {
                memcpy(tmp, vec_, size_ * sizeof(unsigned long long int));
            }
            delete [] vec_;
            vec_ = tmp;
        }
        for (size_t i = size_; i < s; i++) vec_[i] = 0UL;
        size_ = s;
    }
    unsigned long long int* begin() const
    {
        return vec_;
    }
    unsigned long long int* end() const
    {
        return vec_ + size_;
    }
    unsigned long long int& operator[] (size_t i)
    {
        return vec_[i];
    }
    const unsigned long long int& operator[] (size_t i) const
    {
        return vec_[i];
    }
};
typedef FastVector BitMatrixRow;
typedef uint32_t* BitMatrixRowIterator;
typedef FastVector64 BitMatrix64Row;
typedef unsigned long long int* BitMatrix64RowIterator;
#else
typedef std::vector<unsigned long int> BitMatrixRow;
typedef std::vector<unsigned long int>::iterator BitMatrixRowIterator;
#endif
struct BitMatrix
{
    BitMatrix(size_t n = 0L);

    BitMatrix(size_t rows, size_t cols) : cols_(cols)
    {
        if (cols > BitOperations::BITS_IN_WORD) throw "cols must be in 0 ... BITS_IN_WORD";
        row_.resize(rows);
        for (size_t i = 0; i < rows; i++) row_[i] = 0UL;
    }

    BitMatrix(const BitMatrix& bm)
    {
        row_ = bm.row_;
        cols_ = bm.cols_;
    }

    ~BitMatrix()
    {
    }

    BitMatrix& operator=(const BitMatrix& bm)
    {
        row_ = bm.row_;
        cols_ = bm.cols_;
        return *this;
    }

    BitMatrix& operator-=(const BitMatrix& bm)
    {
        if (rows() != bm.rows()) throw "BitMatrix::operator-= : incompatible BitMatrices";
        for (size_t i = 0; i < rows(); i++)
        {
            row_[i] ^= bm.row_[i];
        }
        return *this;
    }

    BitMatrix& operator+=(const BitMatrix& bm)
    {
        if (rows() != bm.rows()) throw "BitMatrix::operator+= : incompatible BitMatrices";
        for (size_t i = 0; i < rows(); i++)
        {
            row_[i] ^= bm.row_[i];
        }
        return *this;
    }

    int operator==(const BitMatrix& bm) const
    {
        if (rows() != bm.rows() || cols() != bm.cols()) return 0;
        unsigned long int colmask = BitOperations::colMask(cols());
        for (size_t i = 0; i < rows(); i++)
        {
            if ((row_[i] & colmask) != (bm.row_[i] & colmask))
            {
                return 0;
            }
        }
        return 1;
    }

    int operator!=(const BitMatrix& bm) const
    {
        return !(*this == bm);
    }

    void exchange(size_t r1, size_t r2)
    {
        if (r1 == r2 || r1 >= rows() || r2 >= rows()) return;
        unsigned long int tmp = row_[r1];
        row_[r1] = row_[r2];
        row_[r2] = tmp;
    }

    int isZero()
    {
        unsigned long int colmask = BitOperations::colMask(cols());
        for (size_t i = 0; i < row_.size(); i++)
        {
            if ((row_[i] & colmask) != 0UL) return 0;
        }
        return 1;
    }

    size_t rows() const
    {
        return row_.size();
    }
    size_t cols() const
    {
        return cols_;
    }

    void xor(size_t row1, size_t row2)
    {
        row_[row1] ^= row_[row2];
    }

    BitMatrixRow row_;
    size_t cols_;
    friend void innerProduct(const BitMatrix& X_T, const BitMatrix& Y, BitMatrix& X_tY);
    friend void multiply(const BitMatrix& bm1, const BitMatrix& bm2, BitMatrix& prod);
    friend void multiply(const BitMatrix& ZL, const BitMatrix& ZR, const BitMatrix& UL, const BitMatrix& UR, BitMatrix& ZUL, BitMatrix& ZUR);
    friend void sym_multiply(const BitMatrix& X, BitMatrix& BBt);
    friend void invert(const BitMatrix& bm, BitMatrix& bm_inv);
    friend void randomise(BitMatrix& bm);
    friend void chooseS(const BitMatrix& T, const BitMatrix& Sim1, BitMatrix& Si, BitMatrix& Winvi);
    friend void kernel(const BitMatrix& M, BitMatrix& kerM);
    friend void kernel(const BitMatrix& L, const BitMatrix& R,
                       BitMatrix& kerL, BitMatrix& kerR);
    friend void kernel(std::vector<BitMatrix>& M, BitMatrix& kerM);
    friend std::ostream& operator<<(std::ostream& os, const BitMatrix& bm);
    friend std::istream& operator>>(std::istream& is, BitMatrix& bm);
    friend MemoryMappedFile& operator>>(MemoryMappedFile& is, BitMatrix& bm);
    void printTranspose(std::ostream& os);
    void printTransposeAsSparseMatrix(std::ostream& os);
    void write(std::ostream& os) const;
    void read(std::istream& is);
    void read(MemoryMappedFile& is);
    void readTranspose(std::istream& is);
    void readTranspose(MemoryMappedFile& is);
};
void kernel(const BitMatrix& M, BitMatrix& kerM);
void kernel(const BitMatrix& L, const BitMatrix& R,
            BitMatrix& kerL, BitMatrix& kerR);
void kernel(std::vector<BitMatrix>& M, BitMatrix& kerM);

struct BitMatrix64
{
    BitMatrix64(size_t n = 0L);

    BitMatrix64(size_t rows, size_t cols) : cols_(cols)
    {
        if (cols > BitOperations64::BITS_IN_WORD) throw "cols must be in 0 ... BITS_IN_WORD";
        row_.resize(rows);
        for (size_t i = 0; i < rows; i++) row_[i] = 0ULL;
    }

    BitMatrix64(const BitMatrix64& bm)
    {
        row_ = bm.row_;
        cols_ = bm.cols_;
    }

    ~BitMatrix64()
    {}

    BitMatrix64& operator=(const BitMatrix64& bm)
    {
        row_ = bm.row_;
        cols_ = bm.cols_;
        return *this;
    }

    BitMatrix64& operator-=(const BitMatrix64& bm)
    {
        if (rows() != bm.rows()) throw "BitMatrix64::operator-= : incompatible BitMatrices";
        for (size_t i = 0; i < rows(); i++)
        {
            row_[i] ^= bm.row_[i];
        }
        return *this;
    }

    BitMatrix64& operator+=(const BitMatrix64& bm)
    {
        if (rows() != bm.rows()) throw "BitMatrix64::operator+= : incompatible BitMatrices";
        for (size_t i = 0; i < rows(); i++)
        {
            row_[i] ^= bm.row_[i];
        }
        return *this;
    }

    int operator==(const BitMatrix64& bm) const
    {
        if (rows() != bm.rows() || cols() != bm.cols()) return 0;
        unsigned long long int colmask = BitOperations64::colMask(cols());
        for (size_t i = 0; i < rows(); i++)
        {
            if ((row_[i] & colmask) != (bm.row_[i] & colmask))
            {
                return 0;
            }
        }
        return 1;
    }

    int operator!=(const BitMatrix64& bm) const
    {
        return !(*this == bm);
    }

    void exchange(size_t r1, size_t r2)
    {
        if (r1 == r2 || r1 >= rows() || r2 >= rows()) return;
        unsigned long long int tmp = row_[r1];
        row_[r1] = row_[r2];
        row_[r2] = tmp;
    }

    int isZero()
    {
        unsigned long long int colmask = BitOperations64::colMask(cols());
        for (size_t i = 0; i < row_.size(); i++)
        {
            if ((row_[i] & colmask) != 0ULL) return 0;
        }
        return 1;
    }

    size_t rows() const
    {
        return row_.size();
    }
    size_t cols() const
    {
        return cols_;
    }

    //std::vector<unsigned long int> row_;
    BitMatrix64Row row_;
    size_t cols_;
    friend void innerProduct(const BitMatrix64& X_T, const BitMatrix64& Y, BitMatrix64& X_tY);
    friend void multiply(const BitMatrix64& bm1, const BitMatrix64& bm2, BitMatrix64& prod);
    friend void multiply(const BitMatrix64& ZL, const BitMatrix64& ZR, const BitMatrix64& UL, const BitMatrix64& UR, BitMatrix64& ZUL, BitMatrix64& ZUR);
    friend void sym_multiply(const BitMatrix64& X, BitMatrix64& BBt);
    friend void invert(const BitMatrix64& bm, BitMatrix64& bm_inv);
    friend void randomise(BitMatrix64& bm);
    friend void chooseS(const BitMatrix64& T, const BitMatrix64& Sim1, BitMatrix64& Si, BitMatrix64& Winvi);
    friend void kernel(const BitMatrix64& M, BitMatrix64& kerM);
    friend void kernel(const BitMatrix64& L, const BitMatrix64& R,
                       BitMatrix64& kerL, BitMatrix64& kerR);
    friend void kernel(std::vector<BitMatrix64>& M, BitMatrix64& kerM);
    friend std::ostream& operator<<(std::ostream& os, const BitMatrix64& bm);
    friend std::istream& operator>>(std::istream& is, BitMatrix64& bm);
    friend MemoryMappedFile& operator>>(MemoryMappedFile& is, BitMatrix64& bm);
    void printTranspose(std::ostream& os);
    void printTransposeAsSparseMatrix(std::ostream& os);
    void write(std::ostream& os) const;
    void read(std::istream& is);
    void read(MemoryMappedFile& is);
    void readTranspose(std::istream& is);
    void readTranspose(MemoryMappedFile& is);
};
void kernel(const BitMatrix64& M, BitMatrix64& kerM);
void kernel(const BitMatrix64& L, const BitMatrix64& R,
            BitMatrix64& kerL, BitMatrix64& kerR);
void kernel(std::vector<BitMatrix64>& M, BitMatrix64& kerM);

class ISparseMatrix
{
public:
    virtual ~ISparseMatrix() {}
    virtual void set_size(int rows, int cols) = 0;
    virtual void clear() = 0;
    virtual size_t rows() const = 0;
    virtual size_t cols() const = 0;
    virtual ISparseRow::xor_status xor(size_t row, size_t col) = 0;
    virtual void set_row(size_t row, const std::vector<size_t>& columns) = 0;
    virtual void set_row(size_t row, const int* first_col, int col_count) = 0;
    virtual ISparseRow::const_iterator begin(long int i) const = 0;
    virtual ISparseRow::const_iterator end(long int i) const = 0;
    virtual size_t row_size(long int i) const = 0;
    virtual void removeEmptyRows() = 0;
    virtual void clear_row(long int row) = 0;
    virtual void set_cols() = 0;
    virtual void set_write_row_count(bool write) = 0;
};

//#define FILE_BASED_SPARSE_MATRIX 1

class SparseMatrix : public ISparseMatrix
{
public:
    SparseMatrix()
        : rows_(0), allocated_rows_(0), cols_(0), sparse_row_(0), write_row_count_(true)
#ifdef FILE_BASED_SPARSE_MATRIX
        ,first_row_allocated_on_disc_(-1), fbsrm_(0), max_rows_in_memory_(2000000)
#endif
    {}
    SparseMatrix(int rows, int cols = 0)
        : rows_(rows), cols_(cols), write_row_count_(true)
#ifdef FILE_BASED_SPARSE_MATRIX
        ,first_row_allocated_on_disc_(-1), fbsrm_(0), max_rows_in_memory_(2000000)
#endif
    {
        allocate();
    }
#ifdef FILE_BASED_SPARSE_MATRIX
    void set_max_rows_in_memory(size_t max_rows_in_memory)
    {
        max_rows_in_memory_ = max_rows_in_memory;
    }
#endif
    ~SparseMatrix()
    {
        if (sparse_row_)
        {
            for (size_t row = 0; row < rows_; row++)
            {
                delete sparse_row_[row];
            }
            delete [] sparse_row_;
        }
#ifdef FILE_BASED_SPARSE_MATRIX
        delete fbsrm_;
#endif
    }

    void set_size(int rows, int cols)
    {
        this->clear();
        rows_ = rows;
        cols_ = cols;
        allocate();
    }

    void clear()
    {
        for (size_t row = 0; row < rows_; row++)
        {
            delete sparse_row_[row];
        }
        delete [] sparse_row_;
        sparse_row_ = 0;
        rows_ = 0;
        allocated_rows_ = 0;
        cols_ = 0;
#ifdef FILE_BASED_SPARSE_MATRIX
        delete fbsrm_;
        fbsrm_ = 0;
        first_row_allocated_on_disc_ = -1;
#endif
    }
    size_t rows() const
    {
        return rows_;
    }
    size_t cols() const
    {
        return cols_;
    }

    ISparseRow::xor_status xor(size_t row, size_t col)
    {
        //if (row < 0) return SparseRow::XOR_FAILED;
        extend(row + 1);
        //ISparseRow::xor_status rc = sparse_row_[row]->xor(col);
        ISparseRow::xor_status rc;
#ifdef FILE_BASED_SPARSE_MATRIX
        if (first_row_allocated_on_disc_ >= 0 && row >= (size_t)first_row_allocated_on_disc_)
        {
            rc = fbsrm_->row(row).xor(col);
        }
        else
#endif
        {
            rc = sparse_row_[row]->xor(col);
        }
        if (SparseRow::XOR_ADDED == rc)
        {
            if (col + 1 > cols_) cols_ = col + 1;
        }
        return rc;
    }

    void set_row(size_t row, const int* first_col, int col_count)
    {
        extend(row + 1);
        int highest_col = -1;
#ifdef FILE_BASED_SPARSE_MATRIX
        if (first_row_allocated_on_disc_ >= 0 && row >= (size_t)first_row_allocated_on_disc_)
        {
            fbsrm_->set_row(row, first_col, col_count);
            highest_col = fbsrm_->row(row).highest_column();
        }
        else
#endif
        {
            sparse_row_[row]->set_row(first_col, col_count);
            highest_col = sparse_row_[row]->highest_column();
        }

        if (col_count)
        {
            if (highest_col + 1 > (int)cols_) cols_ = highest_col + 1;
        }
    }

    void set_row(size_t row, const std::vector<size_t>& columns)
    {
        extend(row + 1);
        int highest_col = -1;
#ifdef FILE_BASED_SPARSE_MATRIX
        if (first_row_allocated_on_disc_ >= 0 && row >= (size_t)first_row_allocated_on_disc_)
        {
            fbsrm_->set_row(row, columns);
            highest_col = sparse_row_[row]->highest_column();
        }
        else
#endif
        {
            sparse_row_[row]->set_row(columns);
            highest_col = sparse_row_[row]->highest_column();
        }

        if (!columns.empty())
        {
            if (highest_col + 1 > (int)cols_) cols_ = highest_col + 1;
        }
    }

    //friend void transpose(const SparseMatrix& A, SparseMatrix& A_t, long int max_row_size = 0, bool clear_A = false);
    friend void transpose(const SparseMatrix& A, SparseMatrix& A_t, long int max_row_size, bool clear_A);
    friend std::ostream& operator<<(std::ostream& os, const SparseMatrix& A);
    friend std::istream& operator>>(std::istream& is, SparseMatrix& A);

    ISparseRow::const_iterator begin(long int i) const
    {
#ifdef FILE_BASED_SPARSE_MATRIX
        if (first_row_allocated_on_disc_ >= 0 && i >= first_row_allocated_on_disc_)
        {
            return fbsrm_->row(i).begin();
        }
#endif

        return sparse_row_[i]->begin();
    }

    ISparseRow::const_iterator end(long int i) const
    {
#ifdef FILE_BASED_SPARSE_MATRIX
        if (first_row_allocated_on_disc_ >= 0 && i >= first_row_allocated_on_disc_)
        {
            return fbsrm_->row(i).end();
        }
#endif

        return sparse_row_[i]->end();
    }

    size_t row_size(long int i) const
    {
#ifdef FILE_BASED_SPARSE_MATRIX
        if (first_row_allocated_on_disc_ >= 0 && i >= first_row_allocated_on_disc_)
        {
            return fbsrm_->row(i).size();
        }
#endif
        if (!sparse_row_ || !sparse_row_[i]) return 0;
        return sparse_row_[i]->size();
    }

    void copy_row(long int i, ISparseRow& copy_of_row) const
    {
#ifdef FILE_BASED_SPARSE_MATRIX
        if (first_row_allocated_on_disc_ >= 0 && i >= first_row_allocated_on_disc_)
        {
            fbsrm_->row(i).copy(copy_of_row);
            return;
        }
#endif
        if (sparse_row_ && sparse_row_[i])
        {
            for (auto it = sparse_row_[i]->begin();
                    it != sparse_row_[i]->end();
                    ++it)
            {
                copy_of_row.add_next(*it);
            }
        }
    }

    void removeEmptyRows();
    void clear_row(long int row)
    {
#ifdef FILE_BASED_SPARSE_MATRIX
        if (first_row_allocated_on_disc_ >= 0 && row >= first_row_allocated_on_disc_)
        {
            fbsrm_->row(row).clear();
            return;
        }
#endif
        sparse_row_[row]->clear();
    }

    size_t memory_usage(const char* s = "");

    void compress();
    void set_cols();
    void read(MemoryMappedFile& is);

    void set_write_row_count(bool write)
    {
        write_row_count_ = write;
    }

private:
    SparseMatrix(const SparseMatrix& sm);
    SparseMatrix& operator=(const SparseMatrix& sm);
    bool operator==(const SparseMatrix& A) const;
    bool operator!=(const SparseMatrix& A) const;

    void allocate()
    {
        allocated_rows_ = 2 * rows_;
        sparse_row_ = new ISparseRow* [ allocated_rows_ ];
        for (size_t row = 0; row < rows_; row++)
        {
#ifdef FILE_BASED_SPARSE_MATRIX
            try
            {
                if (first_row_allocated_on_disc_ >= 0)
                {
                    sparse_row_[row] = 0;
                }
                else
#endif
                {
#ifdef FILE_BASED_SPARSE_MATRIX
                    if (row >= max_rows_in_memory_)
                    {
                        allocate_fbsrm(row);
                    }
                    else
#endif
                    {
                        sparse_row_[row] = new SparseRow(20L);
                    }
                }
            }
#ifdef FILE_BASED_SPARSE_MATRIX
            catch (const std::bad_alloc& e)
            {
                allocate_fbsrm(row);
            }
        }
#endif
    }

    void extend(size_t rows)
    {
        const long int extra = 1000L;
        if (rows > allocated_rows_)
        {
            allocated_rows_ = rows + extra;
            ISparseRow** new_row = new ISparseRow* [ allocated_rows_ ];
            for (size_t row = 0; row < rows_; ++row)
            {
                new_row[row] = sparse_row_[row];
            }
            delete [] sparse_row_;
            sparse_row_ = new_row;
        }
        if (rows > rows_)
        {
            for (size_t row = rows_; row < rows; ++row)
            {
#ifdef FILE_BASED_SPARSE_MATRIX
                try
                {
                    if (first_row_allocated_on_disc_ >= 0)
                    {
                        sparse_row_[row] = 0;
                    }
                    else
#endif
                    {
#ifdef FILE_BASED_SPARSE_MATRIX
                        if (row >= max_rows_in_memory_)
                        {
                            allocate_fbsrm(row);
                        }
                        else
#endif
                        {
                            sparse_row_[row] = new SparseRow;
                        }
                    }
                }
#ifdef FILE_BASED_SPARSE_MATRIX
                catch (const std::bad_alloc& e)
                {
                    allocate_fbsrm(row);
                }
            }
#endif
            rows_ = rows;
        }
    }

#ifdef FILE_BASED_SPARSE_MATRIX
    void allocate_fbsrm(size_t row)
    {
        if (first_row_allocated_on_disc_ < 0)
        {
            first_row_allocated_on_disc_ = row;
            std::string fbsrm_filename(get_tmp_filename());
            fbsrm_ = new FileBasedSparseRowManager(fbsrm_filename, first_row_allocated_on_disc_, 1000);
            sparse_row_[row] = 0;
        }
    }
#endif

#ifdef FILE_BASED_SPARSE_MATRIX
    static std::string get_tmp_filename()
    {
#ifdef WIN32
        char str[MAX_PATH];
        ::GetTempFileName(".", "", 0, str);
        std::cerr << "get_tmp_filename() : str = " << str << std::endl;
        return str;
#else
        uuid_t uuid;
        uuid_generate(uuid);
        char uuid_str[37];
        uuid_unparse(uuid, uuid_str);
        return std::string(uuid_str);
#endif
    }
#endif

    void read(std::istream& is);
    void parse(const std::string& str, size_t row);
    size_t rows_;
    size_t allocated_rows_;
    size_t cols_;
    ISparseRow** sparse_row_;
    bool write_row_count_;
#ifdef FILE_BASED_SPARSE_MATRIX
    long int first_row_allocated_on_disc_;
    FileBasedSparseRowManager* fbsrm_;
    size_t max_rows_in_memory_;
#endif
};
inline SparseRowIterator& SparseRowIterator::operator++()
{
    ++ptr_;
    return *this;
}

inline size_t SparseRowIterator::operator*() const
{
    return *ptr_;
}
inline bool SparseRowIterator::operator==(const SparseRowIterator& sri) const
{
    return (ptr_ == sri.ptr_);
}
inline void SparseRow::compress()
{
    if (one_cap_ > one_size_)
    {
        uint32_t* new_one_ = new uint32_t [ one_size_ ];
        //size_t* new_one_ = new size_t [ one_size_ ];
        memcpy(new_one_, one_, sizeof(uint32_t)*one_size_);
        //memcpy(new_one_, one_, sizeof(size_t)*one_size_);
        one_cap_ = one_size_;
        delete [] one_;
        one_ = new_one_;
    }
}

// A different implementation of SparseMatrix using a list of indices
// perhaps sorted
// This should make operations used in blockLanczos (multiply and multiplyt) more efficient
class SparseMatrix3;
class SparseMatrix2
{
    friend class SparseMatrix3;
public:
    SparseMatrix2();
    SparseMatrix2(size_t allocated_points);
    SparseMatrix2(const std::string& file, bool split = false);
    SparseMatrix2(const std::vector<long int>& points, long int rows, long int columns);
    ~SparseMatrix2();
    void add_row(size_t row, const std::string& str);
    void add_row(size_t row, size_t num_cols, char* s);
    friend void multiply(const SparseMatrix2& A, const BitMatrix& X, BitMatrix& AX);
    friend void multiply(const SparseMatrix2& A, const BitMatrix64& X, BitMatrix64& AX);
    friend void multiplyt(const SparseMatrix2& A, const BitMatrix& X, BitMatrix& AtX);
    friend void multiplyt(const SparseMatrix2& A, const BitMatrix64& X, BitMatrix64& AtX);
    friend void sym_multiply(const SparseMatrix2& B, const BitMatrix& X, BitMatrix& AX);
    friend void sym_multiply(const SparseMatrix2& B, const BitMatrix64& X, BitMatrix64& AX);
    friend std::ostream& operator<<(std::ostream& os, const SparseMatrix2& sm);
    void multiply_dense_part_by_bit_matrix(const BitMatrix& L, BitMatrix& BL) const;
    void multiply_dense_part_by_bit_matrix(const BitMatrix& L, const BitMatrix& R, BitMatrix& BL, BitMatrix& BR) const;
    void multiply_dense_part_by_bit_matrix(const BitMatrix64& L, BitMatrix64& BL) const;
    void multiply_dense_part_by_bit_matrix(const BitMatrix64& L, const BitMatrix64& R, BitMatrix64& BL, BitMatrix64& BR) const;
    void clear();
    size_t rows() const
    {
        return rows_;
    }
    size_t cols() const
    {
        return cols_;
    }

private:
    void write_dense_row(const std::string& str);
    std::fstream* dense_file_;
    bool split_;
    bool parse(const std::string& str, size_t row);
    void extend(size_t row, size_t cols);
    size_t rows_;
    size_t cols_;
    size_t allocated_points_;
    size_t dense_rows_;
    typedef long int Point;
    mutable Point* set_points_;
    Point* next_point_;
    Point* last_point_;
    friend void multiply(const SparseMatrix3& A, const BitMatrix& X, BitMatrix& AX);
    friend void multiply(const SparseMatrix3& A, const BitMatrix64& X, BitMatrix64& AX);
    friend void multiplyt(const SparseMatrix3& A, const BitMatrix& X, BitMatrix& AtX);
    friend void multiplyt(const SparseMatrix3& A, const BitMatrix64& X, BitMatrix64& AtX);
};

class SparseMatrix4
{
    friend class SparseMatrix3;
    friend std::ostream& operator<<(std::ostream& os, const SparseMatrix3& sm);
public:
    SparseMatrix4();
    SparseMatrix4(size_t rows, size_t allocated_points);
    SparseMatrix4(const std::string& file, bool split = false);
    SparseMatrix4(const std::vector<long int>& points, long int rows, long int columns);
    ~SparseMatrix4();
    void add_row(size_t row, const std::string& str);
    void add_row(size_t row, size_t num_cols, char* s);
    friend void multiply(const SparseMatrix4& A, const BitMatrix& X, BitMatrix& AX);
    friend void multiply(const SparseMatrix4& A, const BitMatrix64& X, BitMatrix64& AX);
    friend void multiplyt(const SparseMatrix4& A, const BitMatrix& X, BitMatrix& AtX);
    friend void multiplyt(const SparseMatrix4& A, const BitMatrix64& X, BitMatrix64& AtX);
    friend void sym_multiply(const SparseMatrix4& B, const BitMatrix& X, BitMatrix& AX);
    friend void sym_multiply(const SparseMatrix4& B, const BitMatrix64& X, BitMatrix64& AX);
    friend std::ostream& operator<<(std::ostream& os, const SparseMatrix4& sm);
    void multiply_dense_part_by_bit_matrix(const BitMatrix& L, BitMatrix& BL) const;
    void multiply_dense_part_by_bit_matrix(const BitMatrix& L, const BitMatrix& R, BitMatrix& BL, BitMatrix& BR) const;
    void multiply_dense_part_by_bit_matrix(const BitMatrix64& L, BitMatrix64& BL) const;
    void multiply_dense_part_by_bit_matrix(const BitMatrix64& L, const BitMatrix64& R, BitMatrix64& BL, BitMatrix64& BR) const;
    void clear();
    size_t rows() const
    {
        return rows_;
    }
    size_t cols() const
    {
        return cols_;
    }

private:
    void write_dense_row(const std::string& str);
    std::fstream* dense_file_;
    bool split_;
    bool parse(const std::string& str, size_t row);
    void extend(long int cols);
    long int rows_;
    long int rows_added_;
    long int cols_;
    long int allocated_points_;
    size_t dense_rows_;
    typedef long int Point;
    mutable Point* set_points_;
    Point* next_point_;
    Point* last_point_;
    friend void multiply(const SparseMatrix3& A, const BitMatrix& X, BitMatrix& AX);
    friend void multiply(const SparseMatrix3& A, const BitMatrix64& X, BitMatrix64& AX);
    friend void multiplyt(const SparseMatrix3& A, const BitMatrix& X, BitMatrix& AtX);
    friend void multiplyt(const SparseMatrix3& A, const BitMatrix64& X, BitMatrix64& AtX);
};

// Sparse matrix implementation where the columns are divided into stripes in order
// to be more cache-friendly.
//
class SparseMatrix3
{
public:
    SparseMatrix3(const std::string& file, bool split = false);
    ~SparseMatrix3();
    void clear();
    size_t rows() const
    {
        //return rows_ + very_dense_count_;
        return rows_;
    }
    size_t cols() const
    {
        return cols_;
    }

    friend void multiply(const SparseMatrix3& A, const BitMatrix& X, BitMatrix& AX);
    friend void multiply(const SparseMatrix3& A, const BitMatrix64& X, BitMatrix64& AX);
    friend void multiplyt(const SparseMatrix3& A, const BitMatrix& X, BitMatrix& AtX);
    friend void multiplyt(const SparseMatrix3& A, const BitMatrix64& X, BitMatrix64& AtX);
    friend void sym_multiply(const SparseMatrix3& B, const BitMatrix& X, BitMatrix& AX);
    friend void sym_multiply(const SparseMatrix3& B, const BitMatrix64& X, BitMatrix64& AX);
    void multiply_dense_part_by_bit_matrix(const BitMatrix& L, BitMatrix& BL) const;
    void multiply_dense_part_by_bit_matrix(const BitMatrix& L, const BitMatrix& R, BitMatrix& BL, BitMatrix& BR) const;
    void multiply_dense_part_by_bit_matrix(const BitMatrix64& L, BitMatrix64& BL) const;
    void multiply_dense_part_by_bit_matrix(const BitMatrix64& L, const BitMatrix64& R, BitMatrix64& BL, BitMatrix64& BR) const;

    friend std::ostream& operator<<(std::ostream& os, const SparseMatrix3& sm);

private:
    size_t rows_;
    size_t cols_;
    // The really sparse rows
    SparseMatrix2* sparse_;
    size_t sparse_count_;
    size_t sparse_allocated_points_;

    // Rows which are fairly sparse, divided into a number of stripes
    std::vector<SparseMatrix4* > medium_;
    size_t medium_count_;
    size_t number_of_stripes_;
    std::unordered_map<size_t, size_t> stripe_allocated_points_;

    // Very dense rows which are to be processed later
    std::fstream* very_dense_file_;
    size_t very_dense_count_;

private:
    bool parse(const std::string& str, size_t row);
    bool parse_for_sizing(const std::string& str, long int row);
    void write_very_dense_row(const std::string& str);
    void add_to_medium_dense_rows(long int num_cols, char* s);
    void add_to_size_of_medium_dense_rows(long int num_cols, char* s);
    void extend_dense(size_t stripe);
};
#endif
