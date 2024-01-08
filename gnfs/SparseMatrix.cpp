#include "SparseMatrix.h"
#include <unordered_map>
#include <algorithm>
#include <sstream>
#include "timings.h"
#include <memory>
#include <ctype.h>
#include "Logger.h"

//static Timing Timer("bl.tim");
namespace
{
const double medium_density_bound = 0.0001;
const long int stripe_size = 32768;

const double very_dense_density_bound = 0.4;
const size_t max_very_dense_rows = 50;
};

extern "C" unsigned long int genrand();

std::ostream& operator<<(std::ostream& os, const SparseRow& r)
{
    if (r.size())
    {
        os << static_cast<unsigned int>(r.size());
        for (size_t j = 0; j < r.size(); j++)
        {
            os << " " << static_cast<unsigned int>(r.one_[j]);
        }
        os << std::endl;
    }
    return os;
}

void SparseMatrix::parse(const std::string& str, size_t row)
{
    long int highest_column = 0;
    long int num_cols = SparseRow::begin_parse__(str);
    sparse_row_[row] = 0;
#ifdef FILE_BASED_SPARSE_MATRIX
    if (row >= max_rows_in_memory_)
    {
        allocate_fbsrm(row);
        FileBasedSparseRow r = fbsrm_->row(row, num_cols, s);
        highest_column = r.highest_column();
    }
    else
#endif
    {
        sparse_row_[row] = new SparseRow(num_cols, true);
        highest_column = sparse_row_[row]->highest_column();
    }
    if (highest_column > static_cast<long int>(cols_))
    {
        cols_ = highest_column;
    }
}

void SparseMatrix::read(std::istream& is)
{
    for (size_t row = 0; row < rows_; row++) delete sparse_row_[row];
    delete [] sparse_row_;
    cols_ = 0;
    long int rows = 0;
    std::string str;
    if (!getline(is, str)) return;
    rows_ = std::atol(str.c_str());
    allocated_rows_ = rows_;
    sparse_row_ = new ISparseRow* [ allocated_rows_ ];

    while (getline(is, str))
    {
        parse(str, rows);
        ++rows;
    }
    rows_ = rows;
    ++cols_;
}

void SparseMatrix::read(MemoryMappedFile& mmf)
{
    for (size_t row = 0; row < rows_; row++) delete sparse_row_[row];
    delete [] sparse_row_;
    cols_ = 0;
    long int rows = 0;
    std::string str;
    if (!getline(mmf, str)) return;
    rows_ = std::atol(str.c_str());
    allocated_rows_ = rows_;
    sparse_row_ = new ISparseRow* [ allocated_rows_ ];

    while (getline(mmf, str))
    {
        parse(str, rows);
        ++rows;
    }
    rows_ = rows;
    ++cols_;
}

std::istream& operator>>(std::istream& istr, SparseMatrix& A)
{
    A.read(istr);
    return istr;
}

std::ostream& operator<<(std::ostream& os, const SparseMatrix& sm)
{
    if (sm.write_row_count_)
    {
        os << static_cast<unsigned int>(sm.rows()) << std::endl;
    }
    for (size_t i = 0; i < sm.rows(); i++)
    {
        if (sm.row_size(i))
        {
            os << static_cast<unsigned int>(sm.row_size(i));
            for (auto it = sm.begin(i);
                    it != sm.end(i);
                    ++it)
            {
                os << " " << static_cast<unsigned int>(*it);
            }
        }
        else
        {
            os << '0';
        }
        os << std::endl;
    }
    return os;
}

void transpose(const SparseMatrix& A, SparseMatrix& transA, long int max_row_size, bool clear_A)
{
    LOG_DEBUG("Entering transpose ...");
    for (size_t row = 0; row < transA.rows_; row++) delete transA.sparse_row_[row];
    delete [] transA.sparse_row_;
    transA.sparse_row_ = new ISparseRow* [ A.cols() ];
    transA.rows_ = A.cols();
    transA.allocated_rows_ = A.cols();
    transA.cols_ = A.rows();

    // sizes records the number of entries that
    // will be in each row of transA, which is
    // the same as the number of entries in the
    // corresponding "column" of A
    // initialise it and each row of transA to zero
    LOG_DEBUG("Initialise sizes ... ");
    long int* sizes = new long int [ A.cols() ];
    for (size_t i = 0; i < A.cols(); ++i)
    {
        sizes[i] = 0;
        transA.sparse_row_[i] = 0;
    }

    // read each row of A and update the corresponding
    // entries in sizes
    LOG_DEBUG("Populate sizes ... ");
    for (size_t i = 0; i < A.rows(); i++)
    {
        if (A.row_size(i))
        {
            auto itend = A.end(i);
            for (auto it = A.begin(i);
                    it != itend;
                    ++it)
            {
                size_t row = *it;
                ++sizes[row];
            }
        }
    }

#ifdef FILE_BASED_SPARSE_MATRIX
    if (transA.rows() > transA.max_rows_in_memory_)
    {
        LOG_DEBUG("Allocating FileBasedSparseMatrix, transA.rows() = " << transA.rows() << ", transA.max_rows_in_memory_ = " << transA.max_rows_in_memory_);
        transA.allocate_fbsrm(transA.max_rows_in_memory_);
    }
#endif

    // sizes is now set up
    // read A and populate rows of transA
    for (size_t i = 0; i < A.rows(); i++)
    {
        if (i % 10000 == 0)
        {
            LOG_DEBUG("row = " << i);
        }
        if (A.row_size(i))
        {
            auto itend = A.end(i);
            for (auto it = A.begin(i);
                    it != itend;
                    ++it)
            {
                size_t row = *it;
                if (!transA.sparse_row_[row])
                {
                    if (!max_row_size || sizes[row] <= max_row_size)
                    {
#ifdef FILE_BASED_SPARSE_MATRIX
                        if (row >= transA.max_rows_in_memory_)
                        {
                            transA.allocate_fbsrm(row);
                        }
                        else
#endif
                        {
                            transA.sparse_row_[row] = new SparseRow(sizes[row]);
                        }
                    }
                }
                if (transA.sparse_row_[row])
                {
                    transA.sparse_row_[row]->add_next(i);
                }
                else
                {
                    if (!max_row_size || sizes[row] <= max_row_size)
                    {
#ifdef FILE_BASED_SPARSE_MATRIX
                        if (row >= transA.max_rows_in_memory_)
#endif
                        {
                            transA.do_xor(row, i);
                        }
                    }
                }
            }
        }
        if (clear_A)
        {
            delete A.sparse_row_[i];
            A.sparse_row_[i] = 0;
        }
    }
    delete [] sizes;
    transA.memory_usage();
}

BitMatrix::BitMatrix(size_t n) : cols_(n)
{
    if (n > BitOperations::BITS_IN_WORD) throw "cols must be in 0 ... BITS_IN_WORD";
    // create an n x n identity matrix
    row_.resize(n);
    for (size_t i = 0; i < n; i++)
    {
        row_[i] = 0UL;
        BitOperations::setBit(i, row_[i]);
    }
}

void BitMatrix::printTranspose(std::ostream& os)
{
    // Print transpose of self to os
    for (size_t i = 0; i < cols(); i++)
    {
        for (size_t j = 0; j < rows(); j++)
        {
            if (BitOperations::bitSet(i, row_[j])) os << "1";
            else os << "0";
        }
        os << std::endl;
    }
}

void BitMatrix::printTransposeAsSparseMatrix(std::ostream& os)
{
    // Print transpose of self to os in sparse matrix format (count followed by columns which are set)
    for (size_t i = 0; i < cols(); i++)
    {
        std::ostringstream oss;
        size_t count = 0;
        for (size_t j = 0; j < rows(); j++)
        {
            if (BitOperations::bitSet(i, row_[j]))
            {
                oss << " " << static_cast<unsigned int>(j);
                ++count;
            }
        }
        os << static_cast<unsigned int>(count) << oss.str() << std::endl;
    }
}

void sym_multiply(const BitMatrix& B, BitMatrix& BBt)
{
    if (B.rows() > BitOperations::BITS_IN_WORD) throw "sym_multiply: too many rows, must be <= BitOperations::BITS_IN_WORD";
    //
    // Calculate B B^t, without explicitly storing B^t
    //
    //            m
    //    t      ---
    // (BB )   = \   B   B
    //      ij   /    ik  jk
    //           ---
    //           k=1
    //
    BBt.row_.clear();
    BBt.row_.resize(B.rows());
    BBt.cols_ = B.rows();

    for (size_t i = 0; i < B.rows(); ++i)
    {
        for (size_t j = 0; j < B.rows(); ++j)
        {
            if (BitOperations::bitCount(B.row_[i] & B.row_[j]) & 1)
            {
                BitOperations::setBit(j, BBt.row_[i]);
            }
        }
    }
}

void multiply(const BitMatrix& bm1, const BitMatrix& bm2, BitMatrix& prod)
{
    if (bm1.cols() != bm2.rows()) throw "multiply: Incompatible BitMatrices";
    //
    //           m
    //          ---
    // prod   = \   bm1   bm2
    //     ij   /      ik    kj
    //          ---
    //          k=1
    //
    // where m = bm1.cols() = bm2.rows()
    // and all elements are in GF(2)
    //

    prod.row_.clear();
    prod.row_.resize(bm1.rows());
    prod.cols_ = bm2.cols();
    //
    // bm2 has up to 32 rows.
    // Set up a precalculated 2-d array c[][], 4 x 256
    // which will be used like this:
    // As we iterate through bm1, each row is a 32-bit "selector" for
    // the rows of bm1. Divide the selector into 4 8-bit bytes, numbered
    // 0 - 3, and use this byte number to index into c[][] to give an
    // array of 256 32-bit words. Then use the selector byte to look up
    // in this array to give the result of xor-ing together the
    // appropriate rows from bm2.
    // We only calculate c[][] once, but use it bm1.rows() times, so if
    // bm1 has considerably more than 32 rows this will be more efficient
    // than calculating the xor for each row.
    //
    uint32_t* bm2_iter = bm2.row_.begin();
    uint32_t c[4 * 256] = {0};
    for (size_t j = 0; j < 256; ++j)
    {
        size_t selector = j;
        for (size_t k = 0; k < 8; ++k)
        {
            if (selector & 1)
            {
                c[j] ^= *(bm2_iter + k);
                c[1*256 + j] ^= *(bm2_iter + k + 8);
                c[2*256 + j] ^= *(bm2_iter + k + 16);
                c[3*256 + j] ^= *(bm2_iter + k + 24);
            }
            selector >>= 1;
        }
    }

    uint32_t* bm1_iter = bm1.row_.begin();
    uint32_t* prod_iter = prod.row_.begin();
    for (size_t i = 0; i < bm1.rows(); ++i)
    {
        uint32_t resultRow = 0UL;
        uint32_t selector = *(bm1_iter + i);
        resultRow ^= c[(unsigned char)(selector)];
        resultRow ^= c[1*256 + (unsigned char)(selector >> 8)];
        resultRow ^= c[2*256 + (unsigned char)(selector >> 16)];
        resultRow ^= c[3*256 + (unsigned char)(selector >> 24)];
        *(prod_iter + i) = resultRow;
    }
}

void multiply(const BitMatrix& ZL, const BitMatrix& ZR,
              const BitMatrix& UL, const BitMatrix& UR,
              BitMatrix& ZUL, BitMatrix& ZUR)
{
    /* Multiply the n x 2N matrix Z by the 2N x 2N matrix U
    // where Z is given by ZL and ZR and U by UL and UR.
    // The result goes into ZUL and ZUR.
    // Note that U may have fewer than 2N columns
    //
    // If we define
    //
    // (C )   = (UL)   if i < N
    //   T ij       ij
    //
    // (C )   = (UL)   if i >= N
    //   B ij       ij
    //
    // (D )   = (UR)   if i < N
    //   T ij       ij
    //
    // (D )   = (UR)   if i >= N
    //   B ij       ij
    //
    // then
    //
    // ZU = (ZL ZR) / C  D  \
    //              |  T  T |
    //              | C  D  |
    //              \  B  B /
    //
    //    = (ZL.C  + ZR.C   ZL.D  + ZR.D )
    //           T       B      T       B
    */
    int N = BitOperations::BITS_IN_WORD;
    BitMatrix CT(N, UL.cols());
    BitMatrix CB(ZR.cols(), UL.cols());
    BitMatrix DT(N, UR.cols());
    BitMatrix DB(ZR.cols(), UR.cols());
    for (int i = 0; i < N; i++)
    {
        CT.row_[i] = UL.row_[i];
        DT.row_[i] = UR.row_[i];
    }
    for (size_t i = 0; i < ZR.cols(); i++)
    {
        CB.row_[i] = UL.row_[N + i];
        DB.row_[i] = UR.row_[N + i];
    }

    BitMatrix tmp;
    multiply(ZL, CT, ZUL);
    multiply(ZR, CB, tmp);
    ZUL += tmp;
    multiply(ZL, DT, ZUR);
    multiply(ZR, DB, tmp);
    ZUR += tmp;
}

void innerProduct(const BitMatrix& X, const BitMatrix& Y, BitMatrix& X_tY)
{
    if (X.rows() != Y.rows()) throw "innerProduct: Incompatible BitMatrices";
    // X and Y are n x Ni and n x Nj matrices, stored as BitMatrix objects
    // Result is an Ni x Nj matrix
    //           n
    //  t       ---
    // X Y   = \    X   Y
    //    ij   /     ki  kj
    //          ---
    //          k=1
    //
    // and all elements are in GF(2)
    //
    X_tY.row_.clear();
    X_tY.row_.resize(X.cols());
    X_tY.cols_ = Y.cols();
    //
    // (code stolen from msieve)
    //
    uint32_t c[4 * 256] = {0};
    uint32_t* Y_iter = Y.row_.begin();
    for (uint32_t* X_iter = X.row_.begin();
            X_iter != X.row_.end();
            ++X_iter, ++Y_iter)
    {
        c[(unsigned char)(*X_iter)] ^= *Y_iter;
        c[1*256 + (unsigned char)((*X_iter) >> 8)] ^= *Y_iter;
        c[2*256 + (unsigned char)((*X_iter) >> 16)] ^= *Y_iter;
        c[3*256 + (unsigned char)((*X_iter) >> 24)] ^= *Y_iter;
    }

    for (size_t i = 0; i < 8; ++i)
    {
        uint32_t a0 = 0;
        uint32_t a1 = 0;
        uint32_t a2 = 0;
        uint32_t a3 = 0;
        for (size_t j = 0; j < 256; ++j)
        {
            if ((j >> i) & 1)
            {
                a0 ^= c[j];
                a1 ^= c[1*256 + j];
                a2 ^= c[2*256 + j];
                a3 ^= c[3*256 + j];
            }
        }
        X_tY.row_[i] = a0;
        X_tY.row_[i + 8] = a1;
        X_tY.row_[i + 16] = a2;
        X_tY.row_[i + 24] = a3;
    }
}

// adapted from Algorithm 2.2.2 (Inverse of a Matrix) in Matrix.h
void invert(const BitMatrix& MM, BitMatrix& XX)
{
    if (MM.rows() != MM.cols()) throw "invert: matrix must be square";
    if (MM.rows() == 0) throw "invert: matrix must be non-trivial";
    size_t n = MM.rows();
    BitMatrix M(MM);
    size_t j = 0;
    BitMatrix B(n); // B is n x n identity matrix

    while (j < n)
    {
        size_t i = j;
        while (i < n && !BitOperations::bitSet(j, M.row_[i])) ++i;
        if (i >= n)
        {
            throw "invert: M is not invertible";
        }
        if (i > j)
        {
            for (size_t l = j; l < n; l++)
            {
                int tmp = BitOperations::bitSet(l, M.row_[j]);
                BitOperations::copyBit(BitOperations::bitSet(l, M.row_[i]), l, M.row_[j]);
                BitOperations::copyBit(tmp, l, M.row_[i]);
            }
            uint32_t tmp = B.row_[j];
            B.row_[j] = B.row_[i];
            B.row_[i] = tmp;
        }
        // step 5 [Eliminate]
        std::vector<int> C;
        C.resize(n, 0);
        // since in GF(2) the inverse of a unit is always 1
        //int d = 1;
        for (size_t k = j + 1; k < n; k++)
        {
            C[k] = BitOperations::bitSet(j, M.row_[k]);
        }
        for (size_t k = j + 1; k < n; k++)
        {
            BitOperations::clearBit(j, M.row_[k]);
            for (size_t l = j + 1; l < n; l++)
            {
                if (C[k] && BitOperations::bitSet(l, M.row_[j]))
                {
                    if (BitOperations::bitSet(l, M.row_[k])) BitOperations::clearBit(l, M.row_[k]);
                    else BitOperations::setBit(l, M.row_[k]);
                }
            }
        }
        for (size_t k = j + 1; k < n; k++)
        {
            if (C[k])
            {
                B.row_[k] ^= B.row_[j];
            }
        }
        // Step 2. [Finished?]

        ++j;
    }

    // step 6. [Solve triangular system]
    XX.row_.resize(n);
    XX.cols_ = n;
    for (size_t i = n - 1; i + 1 != 0; --i)
    {
        XX.row_[i] = B.row_[i];
        for (size_t j = i + 1; j < n; j++)
        {
            if (BitOperations::bitSet(j, M.row_[i]))
            {
                XX.row_[i] ^= XX.row_[j];
            }
        }
    }
}

void randomise(BitMatrix& bm)
{
    for (size_t i = 0; i < bm.row_.size(); i++)
    {
        bm.row_[i] = genrand();
    }
}

void chooseS(const BitMatrix& T, const BitMatrix& Sim1, BitMatrix& Si, BitMatrix& Winvi)
{
    // Taken from P. Montgomery "A Block Lanczos Algorithm for Finding Dependencies over GF(2)"
    // ML and MR represent an N x 2N matrix M
    const size_t N = BitOperations::BITS_IN_WORD;
    BitMatrix ML(T);        // T on the left of M
    BitMatrix MR(N);  // I_N on the right of M

    // Number columns of T, with columns selected by Sim1 at the end of the list
    size_t c[N];
    size_t Nim1 = Sim1.cols();
    size_t k = 0;
    size_t l = N - Nim1;
    uint32_t colmask = BitOperations::colMask(Sim1.cols());
    for (size_t i = 0; i < N; i++)
    {
        // column i selected by Sim1 <=> row i of Sim1 is non-zero
        if ((Sim1.row_[i] & colmask) != 0UL)
        {
            c[l] = i;
            ++l;
        }
        else
        {
            c[k] = i;
            ++k;
        }
    }

    std::vector<size_t> S;
    S.clear();

    for (size_t j = 0; j < N; j++)
    {
        for (size_t k = j; k < N && !BitOperations::bitSet(c[j], ML.row_[c[j]]); k++)
        {
            if (BitOperations::bitSet(c[j], ML.row_[c[k]]))
            {
                ML.exchange(c[j], c[k]);
                MR.exchange(c[j], c[k]);
            }
        }
        if (BitOperations::bitSet(c[j], ML.row_[c[j]]))
        {
            S.push_back(c[j]);
            for (size_t i = 0; i < N; i++)
            {
                if (i != c[j])
                {
                    if (BitOperations::bitSet(c[j], ML.row_[i]))
                    {
                        ML.row_[i] ^= ML.row_[c[j]];
                        MR.row_[i] ^= MR.row_[c[j]];
                    }
                }
            }
        }
        else
        {
            for (size_t k = j; k < N && !BitOperations::bitSet(c[j], MR.row_[c[j]]); k++)
            {
                if (BitOperations::bitSet(c[j], MR.row_[c[k]]))
                {
                    ML.exchange(c[j], c[k]);
                    MR.exchange(c[j], c[k]);
                }
            }
            if (!BitOperations::bitSet(c[j], MR.row_[c[j]])) throw "chooseS: failed";
            for (size_t i = 0; i < N; i++)
            {
                if (i != c[j])
                {
                    if (BitOperations::bitSet(c[j], MR.row_[i]))
                    {
                        ML.row_[i] ^= ML.row_[c[j]];
                        MR.row_[i] ^= MR.row_[c[j]];
                    }
                }
            }
            ML.row_[c[j]] = 0UL;
            MR.row_[c[j]] = 0UL;
        }
    }

    Winvi = MR;
    Si.row_.clear();
    Si.row_.resize(N);
    Si.cols_ = S.size();
    for (size_t i = 0; i < Si.cols(); i++)
    {
        BitOperations::setBit(i, Si.row_[S[i]]);
    }
}

void kernel(const BitMatrix& MM, BitMatrix& kerM)
{
    // Find the kernel of MM using Gaussian elimination
    // (based on Algorithm 2.3.1, see Matrix.h)
    size_t m = MM.rows();
    size_t n = MM.cols();
    BitMatrix M(MM);
    size_t r = 0;
    size_t k = 0;
    std::vector<int > c;
    c.resize(m, 0);
    for (size_t i = 0; i < m; i++) c[i] = -1;
    std::vector<int > d;
    d.resize(n, 0);

    while (k < n)
    {
        // 2. [Scan column]
        size_t j = 0;
        while (j < m && (!BitOperations::bitSet(k, M.row_[j]) || c[j] != -1)) ++j;
        if (j >= m)
        {
            ++r;
            d[k] = -1;
        }
        if (j < m)
        {
            // 3. [Eliminate]
            int dd = 0;
            for (size_t i = 0; i < m; i++)
            {
                if (i != j)
                {
                    dd = BitOperations::bitSet(k, M.row_[i]);
                    BitOperations::clearBit(k, M.row_[i]);
                    if (dd) M.row_[i] ^= M.row_[j];
                }
            }
            c[j] = static_cast<int>(k);
            d[k] = static_cast<int>(j);
        }
        // 4. [Finished?]
        ++k;
    }
    // 5. [Output kernel]

    kerM.row_.resize(n);
    kerM.cols_ = r;
    size_t col = 0;
    for (k = 0; k < n; k++)
    {
        if (d[k] == -1)
        {
            for (size_t i = 0; i < n; i++)
            {
                if (d[i] >= 0)
                {
                    int bit = 0;
                    bit = BitOperations::bitSet(k, M.row_[d[i]]);
                    BitOperations::copyBit(bit, col, kerM.row_[i]);
                }
                else if (i == k)
                {
                    BitOperations::setBit(col, kerM.row_[i]);
                }
                else
                {
                    BitOperations::clearBit(col, kerM.row_[i]);
                }
            }
            ++col;
        }
    }
}

void kernel(const BitMatrix& L, const BitMatrix& R,
            BitMatrix& kerL, BitMatrix& kerR)
{
    // L and R are left and right parts of an n x 2N matrix MM
    // Find the kernel of MM using Gaussian elimination
    // (based on Algorithm 2.3.1, see Matrix.h)
    if (L.rows() != R.rows()) throw "kernel: incompatible L and R";
    size_t m = L.rows();
    size_t N = L.cols();
    size_t n = L.cols() + R.cols();
    BitMatrix ML(L);
    BitMatrix MR(R);
    size_t r = 0;
    size_t k = 0;
    std::vector<int > c;
    c.resize(m, 0);
    for (size_t i = 0; i < m; i++) c[i] = -1;
    std::vector<int > d;
    d.resize(n, 0);

    while (k < n)
    {
        // 2. [Scan column]
        size_t j = 0;
        if (k < N)
        {
            while (j < m && (!BitOperations::bitSet(k, ML.row_[j]) || c[j] != -1)) ++j;
        }
        else
        {
            while (j < m && (!BitOperations::bitSet(k - N, MR.row_[j]) || c[j] != -1)) ++j;
        }
        if (j >= m)
        {
            ++r;
            d[k] = -1;
        }
        if (j < m)
        {
            // 3. [Eliminate]
            int dd = 0;
            for (size_t i = 0; i < m; i++)
            {
                if (i != j)
                {
                    if (k < N)
                    {
                        dd = BitOperations::bitSet(k, ML.row_[i]);
                        BitOperations::clearBit(k, ML.row_[i]);
                    }
                    else
                    {
                        dd = BitOperations::bitSet(k - N, MR.row_[i]);
                        BitOperations::clearBit(k - N, MR.row_[i]);
                    }
                    if (dd)
                    {
                        ML.row_[i] ^= ML.row_[j];
                        MR.row_[i] ^= MR.row_[j];
                    }
                }
            }
            c[j] = static_cast<int>(k);
            d[k] = static_cast<int>(j);
        }
        // 4. [Finished?]
        ++k;
    }
    // 5. [Output kernel]

    kerL.row_.resize(n);
    kerR.row_.resize(n);
    if (r <= N)
    {
        kerL.cols_ = r;
        kerR.cols_ = 0;
    }
    else
    {
        kerL.cols_ = N;
        kerR.cols_ = r - N;
    }
    size_t col = 0;
    for (k = 0; k < n; k++)
    {
        if (d[k] == -1)
        {
            for (size_t i = 0; i < n; i++)
            {
                if (d[i] >= 0)
                {
                    int bit = 0;
                    if (k < N) bit = BitOperations::bitSet(k, ML.row_[d[i]]);
                    else bit = BitOperations::bitSet(k - N, MR.row_[d[i]]);
                    if (col < N)
                    {
                        BitOperations::copyBit(bit, col, kerL.row_[i]);
                    }
                    else
                    {
                        BitOperations::copyBit(bit, col - N, kerR.row_[i]);
                    }
                }
                else if (i == k)
                {
                    if (col < N)
                    {
                        BitOperations::setBit(col, kerL.row_[i]);
                    }
                    else
                    {
                        BitOperations::setBit(col - N, kerR.row_[i]);
                    }
                }
                else
                {
                    if (col < N)
                    {
                        BitOperations::clearBit(col, kerL.row_[i]);
                    }
                    else
                    {
                        BitOperations::clearBit(col - N, kerR.row_[i]);
                    }
                }
            }
            ++col;
        }
    }
}

void kernel(std::vector<BitMatrix>& M, BitMatrix& kerM)
{
    if (M.empty())
        throw "kernel: empty list of BitMatrix";
    size_t m = M[0].rows();
    size_t n = 0;
    for (auto& r: M)
    {
        if (r.rows() != m) throw "kernel: incompatible BitMatrix objects";
        n += r.cols();
    }
    size_t N = M[0].cols();
    size_t r = 0;
    size_t k = 0;
    std::vector<int> c;
    c.resize(m, 0);
    for (size_t i = 0; i < m; i++) c[i] = -1;
    std::vector<int> d;
    d.resize(n, 0);

    while (k < n)
    {
        // 2. [Scan column]
        size_t j = 0;
        size_t u = k / N;
        size_t k_ = k % N;
        BitMatrix& Mu = M[u];
        while (j < m &&
                (
                    BitOperations::bitClear(k_, Mu.row_[j])
                    ||
                    c[j] != -1
                )
              )
        {
            ++j;
        }
        if (j >= m)
        {
            ++r;
            d[k] = -1;
        }
        if (j < m)
        {
            // 3. [Eliminate]
            size_t i = 0;
            for (; i < j; ++i)
            {
                if (BitOperations::clearBit(k_, Mu.row_[i]))
                {
                    for (auto& r: M)
                    {
                        r.do_xor(i, j);
                    }
                }
            }
            for (i = j + 1; i < m; ++i)
            {
                if (BitOperations::clearBit(k_, Mu.row_[i]))
                {
                    for (auto& r: M)
                    {
                        r.do_xor(i, j);
                    }
                }
            }
            c[j] = static_cast<int>(k);
            d[k] = static_cast<int>(j);
        }
        // 4. [Finished?]
        ++k;
    }
    // 5. [Output kernel]
    kerM.row_.resize(n);
    for (size_t i = 0; i < n; ++i) kerM.row_[i] = 0UL;
    const size_t max_columns = 16;
    kerM.cols_ = max_columns;
    if (r < max_columns)
        kerM.cols_ = r;

    size_t col = 0;
    for (k = 0; k < n; ++k)
    {
        if (d[k] == -1)
        {
            for (size_t i = 0; i < n; ++i)
            {
                if (d[i] >= 0)
                {
                    int bit = 0;
                    size_t u = k / N;
                    bit = BitOperations::bitSet(k % N, M[u].row_[d[i]]);
                    BitOperations::copyBit(bit, col, kerM.row_[i]);
                }
                else if (i == k)
                {
                    BitOperations::setBit(col, kerM.row_[i]);
                }
                else
                {
                    BitOperations::clearBit(col, kerM.row_[i]);
                }
            }
            ++col;
            if (col >= kerM.cols_)
                break;
        }
    }
}

void BitMatrix::read(std::istream& is)
{
    std::string str;
    if (getline(is, str))
    {
        long int rows = std::atol(str.c_str());
        row_.resize(rows);
    }
    if (getline(is, str))
    {
        cols_ = std::atoi(str.c_str());
    }
    is >> *this;
}

void BitMatrix::read(MemoryMappedFile& is)
{
    std::string str;
    if (getline(is, str))
    {
        long int rows = std::atol(str.c_str());
        row_.resize(rows);
    }
    if (getline(is, str))
    {
        cols_ = std::atoi(str.c_str());
    }
    is >> *this;
}

void BitMatrix::readTranspose(std::istream& is)
{
    // reads a BitMatrix from a file which has been written with printTranspose
    std::string str;
    cols_ = BitOperations::BITS_IN_WORD;

    if (!getline(is, str))
    {
        return;
    }
    size_t rows = str.size();
    row_.resize(rows);
    size_t col = 0;
    bool all_zeros = true;
    for (size_t i = 0; i < rows; ++i)
    {
        if (str[i] == '1')
        {
            BitOperations::setBit(col, row_[i]);
            all_zeros = false;
        }
    }
    if (!all_zeros)
    {
        ++col;
    }

    while (col < cols_ && getline(is, str))
    {
        all_zeros = true;
        for (size_t i = 0; i < rows; ++i)
        {
            if (str[i] == '1')
            {
                BitOperations::setBit(col, row_[i]);
                all_zeros = false;
            }
        }
        if (!all_zeros)
        {
            ++col;
        }
    }

    cols_ = col;
}

void BitMatrix::readTranspose(MemoryMappedFile& is)
{
    // reads a BitMatrix from a file which has been written with printTranspose
    std::string str;
    cols_ = BitOperations::BITS_IN_WORD;

    if (!getline(is, str))
    {
        return;
    }
    size_t rows = str.size();
    row_.resize(rows);
    size_t col = 0;
    bool all_zeros = true;
    for (size_t i = 0; i < rows; ++i)
    {
        if (str[i] == '1')
        {
            BitOperations::setBit(col, row_[i]);
            all_zeros = false;
        }
    }
    if (!all_zeros)
    {
        ++col;
    }

    while (col < cols_ && getline(is, str))
    {
        all_zeros = true;
        for (size_t i = 0; i < rows; ++i)
        {
            if (str[i] == '1')
            {
                BitOperations::setBit(col, row_[i]);
                all_zeros = false;
            }
        }
        if (!all_zeros)
        {
            ++col;
        }
    }

    cols_ = col;
}

void BitMatrix::write(std::ostream& os) const
{
    os << static_cast<unsigned int>(rows()) << std::endl;
    os << static_cast<unsigned int>(cols()) << std::endl;
    os << *this;
}

std::istream& operator>>(std::istream& istr, BitMatrix& bm)
{
    // assume we have already allocated space for rows
    std::string str;
    size_t i = 0;
    while (i < bm.rows() && getline(istr, str))
    {
        // each line consists of up to BITS_IN_WORDS 0s and 1s
        for (size_t j = 0; j < bm.cols(); j++)
        {
            int bit = str.c_str()[j] - '0';
            if (bit) BitOperations::setBit(j, bm.row_[i]);
            else BitOperations::clearBit(j, bm.row_[i]);
        }

        ++i;
    }

    return istr;
}

MemoryMappedFile& operator>>(MemoryMappedFile& istr, BitMatrix& bm)
{
    // assume we have already allocated space for rows
    std::string str;
    size_t i = 0;
    while (i < bm.rows() && getline(istr, str))
    {
        // each line consists of up to BITS_IN_WORDS 0s and 1s
        for (size_t j = 0; j < bm.cols(); j++)
        {
            int bit = str.c_str()[j] - '0';
            if (bit) BitOperations::setBit(j, bm.row_[i]);
            else BitOperations::clearBit(j, bm.row_[i]);
        }

        ++i;
    }

    return istr;
}

std::ostream& operator<<(std::ostream& os, const BitMatrix& bm)
{
    if (bm.cols() == 0) return os;
    for (size_t i = 0; i < bm.rows(); i++)
    {
        for (size_t j = 0; j < bm.cols(); j++)
        {
            if (BitOperations::bitSet(j, bm.row_[i])) os << 1;
            else os << 0;
        }
        os << std::endl;
    }
    return os;
}

void SparseMatrix::removeEmptyRows()
{
    size_t nonEmptyRows = 0;
    for (size_t row = 0; row < rows(); row++)
    {
        if (row_size(row) > 0L) ++nonEmptyRows;
        //if (sparse_row_[row] && sparse_row_[row]->size() > 0L) ++nonEmptyRows;
    }
    if (nonEmptyRows == rows()) return;

    size_t new_row = 0;
    for (size_t row = 0; row < rows(); row++)
    {
        if (row_size(row) > 0L)
            //if (sparse_row_[row] && sparse_row_[row]->size() > 0L)
        {
            if (row > new_row)
            {
#ifdef FILE_BASED_SPARSE_MATRIX
                if (first_row_allocated_on_disc_ >= 0 && row >= first_row_allocated_on_disc_)
                {
                    // two cases:
                    // (i) new_row < first_row_allocated_on_disc_
                    //     copy row from disc file and allocate in memory, remove row from disc file
                    // (ii) new_row >= first_row_allocated_on_disc_
                    //      keep row in disc file but re-index it so it is new_row rather than row
                    if (new_row < first_row_allocated_on_disc_)
                    {
                        delete sparse_row_[new_row];
                        sparse_row_[new_row] = new SparseRow(row_size(row));
                        copy_row(row, *sparse_row_[new_row]);
                        fbsrm_->remove_row(row);
                        // remove row ??? NOT YET IMPLEMENTED ???
                    }
                    else
                    {
                        fbsrm_->replace_row(row, new_row);
                    }
                }
                else
#endif
                {
                    sparse_row_[new_row] = sparse_row_[row];
                    sparse_row_[row] = 0;
                }
            }
            ++new_row;
        }
        else
        {
#ifdef FILE_BASED_SPARSE_MATRIX
            if (first_row_allocated_on_disc_ >= 0 && row >= first_row_allocated_on_disc_)
            {
                // remove row from disc file
                // ??? NOT YET IMPLEMENTED ???
            }
            else
#endif
            {
                delete sparse_row_[row];
                sparse_row_[row] = 0;
            }
        }
    }
    rows_ = nonEmptyRows;
}

size_t SparseMatrix::memory_usage(const char* msg)
{
    // 1 Kb = 1024
    // 1 Mb = 1048576
    size_t s = allocated_rows_ * sizeof(SparseRow*) + sizeof(*this);
    for (size_t row = 0; row < rows_; ++row)
    {
        ISparseRow* srp = sparse_row_[row];
        if (srp)
        {
            s += srp->memory_usage();
        }
    }
    //std::cerr << "SparseMatrix::memory_usage() [" << msg << "] : " << rows() << " x " << cols() << " : " << s << " bytes used" << std::endl;
    return s;
}

void SparseMatrix::set_cols()
{
    long int cols = -1;
    for (size_t row = 0; row < rows_; ++row)
    {
        ISparseRow* srp = sparse_row_[row];
        if (srp)
        {
            if (srp->highest_column() > cols)
            {
                cols = srp->highest_column();
            }
        }
        else
        {
#ifdef FILE_BASED_SPARSE_MATRIX
            if (first_row_allocated_on_disc_ >= 0 && row >= max_rows_in_memory_)
            {
                long int hc = fbsrm_->row(row).highest_column();
                if (hc > cols)
                {
                    cols = hc;
                }
            }
#endif
        }
    }
    cols_ = cols + 1;
}

void SparseMatrix::compress()
{
    // (attempt to) reduce the amount of memory used by this matrix
    const long int extra = 1000L;
    if (allocated_rows_ > rows_ + extra)
    {
        //std::cerr << "SparseMatrix::compress(): allocation = " << allocated_rows_ << ", actual = " << rows_ << std::endl;
    }
    for (size_t row = 0; row < rows_; ++row)
    {
        ISparseRow* srp = sparse_row_[row];
        if (srp)
        {
            srp->compress();
        }
    }

}

SparseMatrix2::SparseMatrix2() : dense_file_(0), split_(false), rows_(0L), cols_(0L), allocated_points_(0), dense_rows_(0L)
    , set_points_(0), next_point_(0), last_point_(0)
{
}

SparseMatrix2::SparseMatrix2(size_t allocated_points) : dense_file_(0), split_(false), rows_(0L), cols_(0L), allocated_points_(allocated_points), dense_rows_(0L)
    , set_points_(0), next_point_(0), last_point_(0)
{
    try
    {
        set_points_ = new Point [ allocated_points_ ];
        next_point_ = set_points_;
        last_point_ = set_points_ + allocated_points_;
    }
    catch (std::bad_alloc& )
    {
        std::ostringstream os;
        os << "SparseMatrix2::SparseMatrix : Caught std::bad_alloc : allocated_points_ = " << static_cast<unsigned int>(allocated_points_);
        std::cerr << os.str() << std::endl;
        throw os.str();
    }
    catch (...)
    {
        std::ostringstream os;
        os << "SparseMatrix2::SparseMatrix : Caught unknown exception : allocated_points_ = " << static_cast<unsigned int>(allocated_points_);
        std::cerr << os.str() << std::endl;
        throw os.str();
    }
}

SparseMatrix2::SparseMatrix2(const std::string& file, bool split) : dense_file_(0), split_(split), rows_(0L), cols_(0L), allocated_points_(0), dense_rows_(0L)
    , set_points_(0), next_point_(0), last_point_(0)
{
    MemoryMappedFile mmf(file.c_str());
    {
        // 1st pass, to find out how much memory we need to allocate
        cols_ = 0;
        std::string str;
        if (!getline(mmf, str))
            return;
        rows_ = std::atol(str.c_str());

        allocated_points_ = 0;
        size_t dense_rows = 0;
        while (getline(mmf, str))
        {
            long int num_cols = SparseRow::begin_parse__(str);

            if (split_ && dense_rows < max_very_dense_rows && (double)num_cols / (double)rows_ > very_dense_density_bound)
            {
                ++dense_rows;
            }
            else
            {
                allocated_points_ += num_cols;
            }
        }
        std::cerr << "SparseMatrix2::SparseMatrix2 : allocated_points_ = " << static_cast<unsigned int>(allocated_points_) << std::endl;
    }

#ifndef USING_GCC
//   allocated_points_ = 1024UL * 145UL * 1024UL;
#else
//   allocated_points_ = 1024UL * 145UL * 1024UL;
#endif

    set_points_ = new Point [ allocated_points_ ];
    next_point_ = set_points_;
    last_point_ = set_points_ + allocated_points_;

    mmf.reset();
    cols_ = 0;
    long int rows = 0;
    std::string str;
    if (!getline(mmf, str))
        return;
    rows_ = std::atol(str.c_str());

    while (getline(mmf, str))
    {
        if (parse(str, rows))
        {
            ++rows;
        }
    }
    rows_ = rows;
    ++cols_;
    last_point_ = next_point_;
    next_point_ = set_points_;
    if (dense_file_)
    {
        delete dense_file_;
        dense_file_ = 0;
    }
}

SparseMatrix2::SparseMatrix2(const std::vector<long int>& points, long int rows, long int columns)
{
    set_points_ = new Point [ points.size() ];
    for (size_t i = 0; i < points.size(); ++i)
    {
        set_points_[i] = points[i];
    }
    last_point_ = set_points_ + points.size();
    rows_ = rows;
    cols_ = columns;
}

SparseMatrix2::~SparseMatrix2()
{
    clear();
}

void multiply(const SparseMatrix2& A, const BitMatrix& X, BitMatrix& AX)
{
    //Timer.start("multiply");
    if (A.cols() != X.rows()) throw "multiply: Incompatible SparseMatrix2 and BitMatrix";
    AX.row_.clear();
    AX.row_.resize(A.rows());
    AX.cols_ = X.cols_;
    uint32_t resultRow = 0UL;
    BitMatrixRowIterator AX_row_iter = AX.row_.begin();
    for (SparseMatrix2::Point* p = A.set_points_; p != A.last_point_; ++p)
    {
        long int col = *p;
        if (col < 0)
        {
            col = -col - 1;
            resultRow ^= *(X.row_.vec_ + col);
            *AX_row_iter = resultRow;
            ++AX_row_iter;
            resultRow = 0UL;
        }
        else
        {
            resultRow ^= *(X.row_.vec_ + col);
        }
    }
    //Timer.stop();
#if 0
    {
        static bool first_time = true;
        if (first_time)
        {
            first_time = false;
            std::fstream f1("X2.dat", std::ios::out);
            X.write(f1);
            std::fstream f2("AX2.dat", std::ios::out);
            AX.write(f2);
        }
    }
#endif
}

void multiply(const SparseMatrix2& A, const BitMatrix64& X, BitMatrix64& AX)
{
    //Timer.start("multiply");
    if (A.cols() != X.rows()) throw "multiply: Incompatible SparseMatrix2 and BitMatrix";
    AX.row_.clear();
    AX.row_.resize(A.rows());
    AX.cols_ = X.cols_;
    unsigned long long int resultRow = 0UL;
    BitMatrix64RowIterator AX_row_iter = AX.row_.begin();
    for (SparseMatrix2::Point* p = A.set_points_; p != A.last_point_; ++p)
    {
        long int col = *p;
        if (col < 0)
        {
            col = -col - 1;
            resultRow ^= *(X.row_.vec_ + col);
            *AX_row_iter = resultRow;
            ++AX_row_iter;
            resultRow = 0UL;
        }
        else
        {
            resultRow ^= *(X.row_.vec_ + col);
        }
    }
    //Timer.stop();
}

void multiplyt(const SparseMatrix2& A, const BitMatrix& X, BitMatrix& AtX)
{
    //
    // Calculate A^tX = ((X^t) A)^t which is a BitMatrix having A.cols() rows and X.cols() columns
    // but without explicitly calculating A^t
    //
    //            m
    //   t       ---               t  t         t
    // (X A)   = \   X   A     = (X A)      = (A X)
    //      ji   /    kj  ki           ij          ij
    //           ---
    //           k=1
    //
    // where m = X.rows() = A.rows()
    // and all elements are in GF(2)
    //
    if (A.rows() != X.rows()) throw "multiplyt: Incompatible SparseMatrix2 and BitMatrix";
    AtX.row_.clear();
    AtX.row_.resize(A.cols());
    AtX.cols_ = X.cols_;
    BitMatrixRowIterator X_row_iter = X.row_.begin();
    for (SparseMatrix2::Point* p = A.set_points_; p != A.last_point_; ++p)
    {
        long int col = *p;
        if (col < 0)
        {
            col = -col - 1;
            AtX.row_[col] ^= *X_row_iter;
            ++X_row_iter;
        }
        else
        {
            AtX.row_[col] ^= *X_row_iter;
        }
    }
#if 0
    {
        static bool first_time = true;
        if (first_time)
        {
            first_time = false;
            std::fstream f1("X.dat", std::ios::out);
            X.write(f1);
            std::fstream f2("AtX.dat", std::ios::out);
            AtX.write(f2);
        }
    }
#endif
}

void multiplyt(const SparseMatrix2& A, const BitMatrix64& X, BitMatrix64& AtX)
{
    //
    // Calculate A^tX = ((X^t) A)^t which is a BitMatrix having A.cols() rows and X.cols() columns
    // but without explicitly calculating A^t
    //
    //            m
    //   t       ---               t  t         t
    // (X A)   = \   X   A     = (X A)      = (A X)
    //      ji   /    kj  ki           ij          ij
    //           ---
    //           k=1
    //
    // where m = X.rows() = A.rows()
    // and all elements are in GF(2)
    //
    if (A.rows() != X.rows()) throw "multiplyt: Incompatible SparseMatrix2 and BitMatrix64";
    AtX.row_.clear();
    AtX.row_.resize(A.cols());
    AtX.cols_ = X.cols_;
    BitMatrix64RowIterator X_row_iter = X.row_.begin();
    for (SparseMatrix2::Point* p = A.set_points_; p != A.last_point_; ++p)
    {
        long int col = *p;
        if (col < 0)
        {
            col = -col - 1;
            AtX.row_[col] ^= *X_row_iter;
            ++X_row_iter;
        }
        else
        {
            AtX.row_[col] ^= *X_row_iter;
        }
    }
}

void sym_multiply(const SparseMatrix2& B, const BitMatrix& X, BitMatrix& AX)
{
    BitMatrix BX;
    multiply(B, X, BX);
    multiplyt(B, BX, AX);
}

void sym_multiply(const SparseMatrix2& B, const BitMatrix64& X, BitMatrix64& AX)
{
    BitMatrix64 BX;
    multiply(B, X, BX);
    multiplyt(B, BX, AX);
}

void SparseMatrix2::write_dense_row(const std::string& str)
{
    if (!dense_file_)
    {
        dense_file_ = new std::fstream("dense_rows.txt", std::ios::out);
    }
    *dense_file_ << str << std::endl;
}

void SparseMatrix2::add_row(size_t row, size_t num_cols)
{
    if (next_point_ + num_cols > last_point_)
    {
        //std::cerr << "Extending: num_cols = " << num_cols << std::endl;
        extend(row, num_cols);
    }
    while (char* s = strtok(0, " "))
    {
        size_t n = std::atol(s);
        *next_point_ = static_cast<long int>(n);
        ++next_point_;
        if (n > cols_)
        {
            cols_ = n;
        }
    }
    if (num_cols > 0)
    {
        *(next_point_ - 1) = - *(next_point_ - 1) - 1;
    }
    ++rows_;
}

void SparseMatrix2::add_row(size_t row, const std::string& str)
{
    long int num_cols = SparseRow::begin_parse__(str);

    add_row(row, num_cols);
}

bool SparseMatrix2::parse(const std::string& str, size_t row)
{
    long int num_cols = SparseRow::begin_parse__(str);

    if (split_ && dense_rows_ < max_very_dense_rows && (double)num_cols / (double)rows_ > very_dense_density_bound)
    {
        write_dense_row(str);
        ++dense_rows_;
        return false;
    }
    if (next_point_ + num_cols > last_point_)
    {
        //std::cerr << "Extending: num_cols = " << num_cols << std::endl;
        extend(row, num_cols);
    }

    while (char* s = strtok(0, " "))
    {
        size_t n = std::atol(s);
        *next_point_ = static_cast<long int>(n);
        ++next_point_;
        if (n > cols_)
        {
            cols_ = n;
        }
    }
    if (num_cols > 0)
    {
        *(next_point_ - 1) = - *(next_point_ - 1) - 1;
    }
    return true;
}

void SparseMatrix2::extend(size_t row, size_t cols)
{
    size_t rows_left = rows_ - row;
    if (rows_ < row)
    {
        rows_left = 100UL;
    }
    std::cerr << "SparseMatrix2::extend() : row = " << row << ", cols = " << cols << ", rows_left = " << rows_left << std::endl;
    std::cerr << "allocated_points_ = " << allocated_points_ << std::endl;
    allocated_points_ += std::max<long int>(rows_left * 20, cols + rows_left);
    std::cerr << "allocated_points_ = " << allocated_points_ << std::endl;
    size_t offset = next_point_ - set_points_;
    Point* new_set_points = new Point [ allocated_points_ ];
    next_point_ = new_set_points + offset;
    last_point_ = new_set_points + allocated_points_;
    memcpy(new_set_points, set_points_, sizeof(Point)*offset);
    delete [] set_points_;
    set_points_ = new_set_points;
}

void SparseMatrix2::clear()
{
    delete [] set_points_;
}

void SparseMatrix2::multiply_dense_part_by_bit_matrix(const BitMatrix& L, const BitMatrix& R, BitMatrix& BL, BitMatrix& BR) const
{
    if (!split_)
        return;

    delete [] set_points_;
    set_points_ = 0;

    BL.row_.clear();
    BL.row_.resize(dense_rows_);
    BL.cols_ = L.cols();
    BR.row_.clear();
    BR.row_.resize(dense_rows_);
    BR.cols_ = R.cols();

    MemoryMappedFile mmf("dense_rows.txt");
    std::string str;

    BitMatrixRowIterator BL_row_iter = BL.row_.begin();
    BitMatrixRowIterator BR_row_iter = BR.row_.begin();

    while (getline(mmf, str))
    {
        uint32_t resultRowL = 0UL;
        uint32_t resultRowR = 0UL;
        SparseRow::begin_parse__(str);
        while (char* s = strtok(0, " "))
        {
            long int n = std::atol(s);
            resultRowL ^= *(L.row_.vec_ + n);
            resultRowR ^= *(R.row_.vec_ + n);
        }
        *BL_row_iter = resultRowL;
        *BR_row_iter = resultRowR;
        ++BL_row_iter;
        ++BR_row_iter;
    }
}

std::ostream& operator<<(std::ostream& os, const SparseMatrix2& sm)
{
    long int cols = 0;
    std::ostringstream oss;
    for (SparseMatrix2::Point* p = sm.set_points_; p != sm.last_point_; ++p)
    {
        long int col = *p;
        if (col < 0)
        {
            col = -col - 1;
            oss << " " << col;
            ++cols;
            os << cols << oss.str() << std::endl;
            oss.str("");
            cols = 0;
        }
        else
        {
            oss << " " << col;
            ++cols;
        }
    }
    return os;
}

void SparseMatrix2::multiply_dense_part_by_bit_matrix(const BitMatrix& L, BitMatrix& BL) const
{
    if (!split_)
        return;

    delete [] set_points_;
    set_points_ = 0;

    BL.row_.clear();
    BL.row_.resize(dense_rows_);
    BL.cols_ = L.cols();

    MemoryMappedFile mmf("dense_rows.txt");
    std::string str;

    BitMatrixRowIterator BL_row_iter = BL.row_.begin();

    while (getline(mmf, str))
    {
        SparseRow::begin_parse__(str);
        uint32_t resultRow = 0UL;

        while (char* s = strtok(0, " "))
        {
            long int n = std::atol(s);
            resultRow ^= *(L.row_.vec_ + n);
        }
        *BL_row_iter = resultRow;
        ++BL_row_iter;
    }
}

BitMatrix64::BitMatrix64(size_t n) : cols_(n)
{
    if (n > BitOperations64::BITS_IN_WORD) throw "cols must be in 0 ... BITS_IN_WORD";
    // create an n x n identity matrix
    row_.resize(n);
    for (size_t i = 0; i < n; i++)
    {
        row_[i] = 0UL;
        BitOperations64::setBit(i, row_[i]);
    }
}

void BitMatrix64::printTranspose(std::ostream& os)
{
    // Print transpose of self to os
    for (size_t i = 0; i < cols(); i++)
    {
        for (size_t j = 0; j < rows(); j++)
        {
            if (BitOperations64::bitSet(i, row_[j])) os << "1";
            else os << "0";
        }
        os << std::endl;
    }
}

void BitMatrix64::printTransposeAsSparseMatrix(std::ostream& os)
{
    // Print transpose of self to os in sparse matrix format (count followed by columns which are set)
    for (size_t i = 0; i < cols(); i++)
    {
        std::ostringstream oss;
        size_t count = 0;
        for (size_t j = 0; j < rows(); j++)
        {
            if (BitOperations64::bitSet(i, row_[j]))
            {
                oss << " " << static_cast<unsigned int>(j);
                ++count;
            }
        }
        os << static_cast<unsigned int>(count) << oss.str() << std::endl;
    }
}

void sym_multiply(const BitMatrix64& B, BitMatrix64& BBt)
{
    if (B.rows() > BitOperations64::BITS_IN_WORD) throw "sym_multiply: too many rows, must be <= BitOperations64::BITS_IN_WORD";
    //
    // Calculate B B^t, without explicitly storing B^t
    //
    //            m
    //    t      ---
    // (BB )   = \   B   B
    //      ij   /    ik  jk
    //           ---
    //           k=1
    //
    BBt.row_.clear();
    BBt.row_.resize(B.rows());
    BBt.cols_ = B.rows();

    for (size_t i = 0; i < B.rows(); ++i)
    {
        for (size_t j = 0; j < B.rows(); ++j)
        {
            if (BitOperations64::bitCount(B.row_[i] & B.row_[j]) & 1)
            {
                BitOperations64::setBit(j, BBt.row_[i]);
            }
        }
    }
}

void multiply(const BitMatrix64& bm1, const BitMatrix64& bm2, BitMatrix64& prod)
{
    if (bm1.cols() != bm2.rows()) throw "multiply: Incompatible BitMatrices";
    //
    //           m
    //          ---
    // prod   = \   bm1   bm2
    //     ij   /      ik    kj
    //          ---
    //          k=1
    //
    // where m = bm1.cols() = bm2.rows()
    // and all elements are in GF(2)
    //

    prod.row_.clear();
    prod.row_.resize(bm1.rows());
    prod.cols_ = bm2.cols();
    //
    // bm2 has up to 64 rows.
    // Set up a precalculated 2-d array c[][], 8 x 256
    // which will be used like this:
    // As we iterate through bm1, each row is a 64-bit "selector" for
    // the rows of bm1. Divide the selector into 8 8-bit bytes, numbered
    // 0 - 7, and use this byte number to index into c[][] to give an
    // array of 256 64-bit words. Then use the selector byte to look up
    // in this array to give the result of xor-ing together the
    // appropriate rows from bm2.
    // We only calculate c[][] once, but use it bm1.rows() times, so if
    // bm1 has considerably more than 64 rows this will be more efficient
    // than calculating the xor for each row.
    //
    unsigned long long int* bm2_iter = bm2.row_.begin();
    unsigned long long int c[8 * 256] = {0};
    for (size_t j = 0; j < 256; ++j)
    {
        size_t selector = j;
        for (size_t k = 0; k < 8; ++k)
        {
            if (selector & 1)
            {
                c[j] ^= *(bm2_iter + k);
                c[1*256 + j] ^= *(bm2_iter + k + 8);
                c[2*256 + j] ^= *(bm2_iter + k + 16);
                c[3*256 + j] ^= *(bm2_iter + k + 24);
                c[4*256 + j] ^= *(bm2_iter + k + 32);
                c[5*256 + j] ^= *(bm2_iter + k + 40);
                c[6*256 + j] ^= *(bm2_iter + k + 48);
                c[7*256 + j] ^= *(bm2_iter + k + 56);
            }
            selector >>= 1;
        }
    }

    unsigned long long int* bm1_iter = bm1.row_.begin();
    unsigned long long int* prod_iter = prod.row_.begin();
    for (size_t i = 0; i < bm1.rows(); ++i)
    {
        unsigned long long int resultRow = 0UL;
        unsigned long long int selector = *(bm1_iter + i);
        resultRow ^= c[(unsigned char)selector];
        resultRow ^= c[1*256 + (unsigned char)(selector >> 8)];
        resultRow ^= c[2*256 + (unsigned char)(selector >> 16)];
        resultRow ^= c[3*256 + (unsigned char)(selector >> 24)];
        resultRow ^= c[4*256 + (unsigned char)(selector >> 32)];
        resultRow ^= c[5*256 + (unsigned char)(selector >> 40)];
        resultRow ^= c[6*256 + (unsigned char)(selector >> 48)];
        resultRow ^= c[7*256 + (unsigned char)(selector >> 56)];
        *(prod_iter + i) = resultRow;
    }
}

void multiply(const BitMatrix64& ZL, const BitMatrix64& ZR,
              const BitMatrix64& UL, const BitMatrix64& UR,
              BitMatrix64& ZUL, BitMatrix64& ZUR)
{
    /* Multiply the n x 2N matrix Z by the 2N x 2N matrix U
    // where Z is given by ZL and ZR and U by UL and UR.
    // The result goes into ZUL and ZUR.
    // Note that U may have fewer than 2N columns
    //
    // If we define
    //
    // (C )   = (UL)   if i < N
    //   T ij       ij
    //
    // (C )   = (UL)   if i >= N
    //   B ij       ij
    //
    // (D )   = (UR)   if i < N
    //   T ij       ij
    //
    // (D )   = (UR)   if i >= N
    //   B ij       ij
    //
    // then
    //
    // ZU = (ZL ZR) / C  D  \
    //              |  T  T |
    //              | C  D  |
    //              \  B  B /
    //
    //    = (ZL.C  + ZR.C   ZL.D  + ZR.D )
    //           T       B      T       B
    */
    int N = BitOperations64::BITS_IN_WORD;
    BitMatrix64 CT(N, UL.cols());
    BitMatrix64 CB(ZR.cols(), UL.cols());
    BitMatrix64 DT(N, UR.cols());
    BitMatrix64 DB(ZR.cols(), UR.cols());
    for (int i = 0; i < N; i++)
    {
        CT.row_[i] = UL.row_[i];
        DT.row_[i] = UR.row_[i];
    }
    for (size_t i = 0; i < ZR.cols(); i++)
    {
        CB.row_[i] = UL.row_[N + i];
        DB.row_[i] = UR.row_[N + i];
    }

    BitMatrix64 tmp;
    multiply(ZL, CT, ZUL);
    multiply(ZR, CB, tmp);
    ZUL += tmp;
    multiply(ZL, DT, ZUR);
    multiply(ZR, DB, tmp);
    ZUR += tmp;
}

void innerProduct(const BitMatrix64& X, const BitMatrix64& Y, BitMatrix64& X_tY)
{
    if (X.rows() != Y.rows()) throw "innerProduct: Incompatible BitMatrices";
    // X and Y are n x Ni and n x Nj matrices, stored as BitMatrix64 objects
    // Result is an Ni x Nj matrix
    //           n
    //  t       ---
    // X Y   = \    X   Y
    //    ij   /     ki  kj
    //          ---
    //          k=1
    //
    // and all elements are in GF(2)
    //
    X_tY.row_.clear();
    X_tY.row_.resize(X.cols());
    X_tY.cols_ = Y.cols();
    //
    // (code stolen from msieve)
    //
    unsigned long long int c[8 * 256] = {0};
    unsigned long long int* Y_iter = Y.row_.begin();
    for (unsigned long long int* X_iter = X.row_.begin();
            X_iter != X.row_.end();
            ++X_iter, ++Y_iter)
    {
        c[(unsigned char)(*X_iter)] ^= *Y_iter;
        c[1*256 + (unsigned char)((*X_iter) >> 8)] ^= *Y_iter;
        c[2*256 + (unsigned char)((*X_iter) >> 16)] ^= *Y_iter;
        c[3*256 + (unsigned char)((*X_iter) >> 24)] ^= *Y_iter;
        c[4*256 + (unsigned char)((*X_iter) >> 32)] ^= *Y_iter;
        c[5*256 + (unsigned char)((*X_iter) >> 40)] ^= *Y_iter;
        c[6*256 + (unsigned char)((*X_iter) >> 48)] ^= *Y_iter;
        c[7*256 + (unsigned char)((*X_iter) >> 56)] ^= *Y_iter;
    }

    for (size_t i = 0; i < 8; ++i)
    {
        unsigned long long int a0 = 0;
        unsigned long long int a1 = 0;
        unsigned long long int a2 = 0;
        unsigned long long int a3 = 0;
        unsigned long long int a4 = 0;
        unsigned long long int a5 = 0;
        unsigned long long int a6 = 0;
        unsigned long long int a7 = 0;
        for (size_t j = 0; j < 256; ++j)
        {
            if ((j >> i) & 1)
            {
                a0 ^= c[j];
                a1 ^= c[1*256 + j];
                a2 ^= c[2*256 + j];
                a3 ^= c[3*256 + j];
                a4 ^= c[4*256 + j];
                a5 ^= c[5*256 + j];
                a6 ^= c[6*256 + j];
                a7 ^= c[7*256 + j];
            }
        }
        X_tY.row_[i] = a0;
        X_tY.row_[i + 8] = a1;
        X_tY.row_[i + 16] = a2;
        X_tY.row_[i + 24] = a3;
        X_tY.row_[i + 32] = a4;
        X_tY.row_[i + 40] = a5;
        X_tY.row_[i + 48] = a6;
        X_tY.row_[i + 56] = a7;
    }
}

// adapted from Algorithm 2.2.2 (Inverse of a Matrix) in Matrix.h
void invert(const BitMatrix64& MM, BitMatrix64& XX)
{
    if (MM.rows() != MM.cols()) throw "invert: matrix must be square";
    if (MM.rows() == 0) throw "invert: matrix must be non-trivial";
    size_t n = MM.rows();
    BitMatrix64 M(MM);
    size_t j = 0;
    BitMatrix64 B(n); // B is n x n identity matrix

    while (j < n)
    {
        size_t i = j;
        while (i < n && !BitOperations64::bitSet(j, M.row_[i])) ++i;
        if (i >= n)
        {
            throw "invert: M is not invertible";
        }
        if (i > j)
        {
            for (size_t l = j; l < n; l++)
            {
                int tmp = BitOperations64::bitSet(l, M.row_[j]);
                BitOperations64::copyBit(BitOperations64::bitSet(l, M.row_[i]), l, M.row_[j]);
                BitOperations64::copyBit(tmp, l, M.row_[i]);
            }
            unsigned long long int tmp = B.row_[j];
            B.row_[j] = B.row_[i];
            B.row_[i] = tmp;
        }
        // step 5 [Eliminate]
        std::vector<int> C;
        C.resize(n, 0);
        // since in GF(2) the inverse of a unit is always 1
        //int d = 1;
        for (size_t k = j + 1; k < n; k++)
        {
            C[k] = BitOperations64::bitSet(j, M.row_[k]);
        }
        for (size_t k = j + 1; k < n; k++)
        {
            BitOperations64::clearBit(j, M.row_[k]);
            for (size_t l = j + 1; l < n; l++)
            {
                if (C[k] && BitOperations64::bitSet(l, M.row_[j]))
                {
                    if (BitOperations64::bitSet(l, M.row_[k])) BitOperations64::clearBit(l, M.row_[k]);
                    else BitOperations64::setBit(l, M.row_[k]);
                }
            }
        }
        for (size_t k = j + 1; k < n; k++)
        {
            if (C[k])
            {
                B.row_[k] ^= B.row_[j];
            }
        }
        // Step 2. [Finished?]

        ++j;
    }

    // step 6. [Solve triangular system]
    XX.row_.resize(n);
    XX.cols_ = n;
    for (size_t i = n - 1; i + 1 != 0; --i)
    {
        XX.row_[i] = B.row_[i];
        for (size_t j = i + 1; j < n; j++)
        {
            if (BitOperations64::bitSet(j, M.row_[i]))
            {
                XX.row_[i] ^= XX.row_[j];
            }
        }
    }
}

void randomise(BitMatrix64& bm)
{
    for (size_t i = 0; i < bm.row_.size(); i++)
    {
        bm.row_[i] = genrand() + ((unsigned long long int)(genrand()) << 32);
    }
}

void chooseS(const BitMatrix64& T, const BitMatrix64& Sim1, BitMatrix64& Si, BitMatrix64& Winvi)
{
    // Taken from P. Montgomery "A Block Lanczos Algorithm for Finding Dependencies over GF(2)"
    // ML and MR represent an N x 2N matrix M
    const size_t N = BitOperations64::BITS_IN_WORD;
    BitMatrix64 ML(T);        // T on the left of M
    BitMatrix64 MR(N);  // I_N on the right of M

    // Number columns of T, with columns selected by Sim1 at the end of the list
    size_t c[N];
    size_t Nim1 = Sim1.cols();
    size_t k = 0;
    size_t l = N - Nim1;
    unsigned long long int colmask = BitOperations64::colMask(Sim1.cols());
    for (size_t i = 0; i < N; i++)
    {
        // column i selected by Sim1 <=> row i of Sim1 is non-zero
        if ((Sim1.row_[i] & colmask) != 0UL)
        {
            c[l] = i;
            ++l;
        }
        else
        {
            c[k] = i;
            ++k;
        }
    }

    std::vector<size_t> S;
    S.clear();

    for (size_t j = 0; j < N; j++)
    {
        for (size_t k = j; k < N && !BitOperations64::bitSet(c[j], ML.row_[c[j]]); k++)
        {
            if (BitOperations64::bitSet(c[j], ML.row_[c[k]]))
            {
                ML.exchange(c[j], c[k]);
                MR.exchange(c[j], c[k]);
            }
        }
        if (BitOperations64::bitSet(c[j], ML.row_[c[j]]))
        {
            S.push_back(c[j]);
            for (size_t i = 0; i < N; i++)
            {
                if (i != c[j])
                {
                    if (BitOperations64::bitSet(c[j], ML.row_[i]))
                    {
                        ML.row_[i] ^= ML.row_[c[j]];
                        MR.row_[i] ^= MR.row_[c[j]];
                    }
                }
            }
        }
        else
        {
            for (size_t k = j; k < N && !BitOperations64::bitSet(c[j], MR.row_[c[j]]); k++)
            {
                if (BitOperations64::bitSet(c[j], MR.row_[c[k]]))
                {
                    ML.exchange(c[j], c[k]);
                    MR.exchange(c[j], c[k]);
                }
            }
            if (!BitOperations64::bitSet(c[j], MR.row_[c[j]])) throw "chooseS: failed";
            for (size_t i = 0; i < N; i++)
            {
                if (i != c[j])
                {
                    if (BitOperations64::bitSet(c[j], MR.row_[i]))
                    {
                        ML.row_[i] ^= ML.row_[c[j]];
                        MR.row_[i] ^= MR.row_[c[j]];
                    }
                }
            }
            ML.row_[c[j]] = 0UL;
            MR.row_[c[j]] = 0UL;
        }
    }

    Winvi = MR;
    Si.row_.clear();
    Si.row_.resize(N);
    Si.cols_ = S.size();
    for (size_t i = 0; i < Si.cols(); i++)
    {
        BitOperations64::setBit(i, Si.row_[S[i]]);
    }
}

void kernel(const BitMatrix64& MM, BitMatrix64& kerM)
{
    // Find the kernel of MM using Gaussian elimination
    // (based on Algorithm 2.3.1, see Matrix.h)
    size_t m = MM.rows();
    size_t n = MM.cols();
    BitMatrix64 M(MM);
    size_t r = 0;
    size_t k = 0;
    std::vector<int > c;
    c.resize(m, 0);
    for (size_t i = 0; i < m; i++) c[i] = -1;
    std::vector<int > d;
    d.resize(n, 0);

    while (k < n)
    {
        // 2. [Scan column]
        size_t j = 0;
        while (j < m && (!BitOperations64::bitSet(k, M.row_[j]) || c[j] != -1)) ++j;
        if (j >= m)
        {
            ++r;
            d[k] = -1;
        }
        if (j < m)
        {
            // 3. [Eliminate]
            int dd = 0;
            for (size_t i = 0; i < m; i++)
            {
                if (i != j)
                {
                    dd = BitOperations64::bitSet(k, M.row_[i]);
                    BitOperations64::clearBit(k, M.row_[i]);
                    if (dd) M.row_[i] ^= M.row_[j];
                }
            }
            c[j] = static_cast<int>(k);
            d[k] = static_cast<int>(j);
        }
        // 4. [Finished?]
        ++k;
    }
    // 5. [Output kernel]

    kerM.row_.resize(n);
    kerM.cols_ = r;
    size_t col = 0;
    for (k = 0; k < n; k++)
    {
        if (d[k] == -1)
        {
            for (size_t i = 0; i < n; i++)
            {
                if (d[i] >= 0)
                {
                    int bit = 0;
                    bit = BitOperations64::bitSet(k, M.row_[d[i]]);
                    BitOperations64::copyBit(bit, col, kerM.row_[i]);
                }
                else if (i == k)
                {
                    BitOperations64::setBit(col, kerM.row_[i]);
                }
                else
                {
                    BitOperations64::clearBit(col, kerM.row_[i]);
                }
            }
            ++col;
        }
    }
}

void kernel(const BitMatrix64& L, const BitMatrix64& R,
            BitMatrix64& kerL, BitMatrix64& kerR)
{
    // L and R are left and right parts of an n x 2N matrix MM
    // Find the kernel of MM using Gaussian elimination
    // (based on Algorithm 2.3.1, see Matrix.h)
    if (L.rows() != R.rows()) throw "kernel: incompatible L and R";
    size_t m = L.rows();
    size_t N = L.cols();
    size_t n = L.cols() + R.cols();
    BitMatrix64 ML(L);
    BitMatrix64 MR(R);
    size_t r = 0;
    size_t k = 0;
    std::vector<int > c;
    c.resize(m, 0);
    for (size_t i = 0; i < m; i++) c[i] = -1;
    std::vector<int > d;
    d.resize(n, 0);

    while (k < n)
    {
        // 2. [Scan column]
        size_t j = 0;
        if (k < N)
        {
            while (j < m && (!BitOperations64::bitSet(k, ML.row_[j]) || c[j] != -1)) ++j;
        }
        else
        {
            while (j < m && (!BitOperations64::bitSet(k - N, MR.row_[j]) || c[j] != -1)) ++j;
        }
        if (j >= m)
        {
            ++r;
            d[k] = -1;
        }
        if (j < m)
        {
            // 3. [Eliminate]
            int dd = 0;
            for (size_t i = 0; i < m; i++)
            {
                if (i != j)
                {
                    if (k < N)
                    {
                        dd = BitOperations64::bitSet(k, ML.row_[i]);
                        BitOperations64::clearBit(k, ML.row_[i]);
                    }
                    else
                    {
                        dd = BitOperations64::bitSet(k - N, MR.row_[i]);
                        BitOperations64::clearBit(k - N, MR.row_[i]);
                    }
                    if (dd)
                    {
                        ML.row_[i] ^= ML.row_[j];
                        MR.row_[i] ^= MR.row_[j];
                    }
                }
            }
            c[j] = static_cast<int>(k);
            d[k] = static_cast<int>(j);
        }
        // 4. [Finished?]
        ++k;
    }
    // 5. [Output kernel]

    kerL.row_.resize(n);
    kerR.row_.resize(n);
    if (r <= N)
    {
        kerL.cols_ = r;
        kerR.cols_ = 0;
    }
    else
    {
        kerL.cols_ = N;
        kerR.cols_ = r - N;
    }
    size_t col = 0;
    for (k = 0; k < n; k++)
    {
        if (d[k] == -1)
        {
            for (size_t i = 0; i < n; i++)
            {
                if (d[i] >= 0)
                {
                    int bit = 0;
                    if (k < N) bit = BitOperations64::bitSet(k, ML.row_[d[i]]);
                    else bit = BitOperations64::bitSet(k - N, MR.row_[d[i]]);
                    if (col < N)
                    {
                        BitOperations64::copyBit(bit, col, kerL.row_[i]);
                    }
                    else
                    {
                        BitOperations64::copyBit(bit, col - N, kerR.row_[i]);
                    }
                }
                else if (i == k)
                {
                    if (col < N)
                    {
                        BitOperations64::setBit(col, kerL.row_[i]);
                    }
                    else
                    {
                        BitOperations64::setBit(col - N, kerR.row_[i]);
                    }
                }
                else
                {
                    if (col < N)
                    {
                        BitOperations64::clearBit(col, kerL.row_[i]);
                    }
                    else
                    {
                        BitOperations64::clearBit(col - N, kerR.row_[i]);
                    }
                }
            }
            ++col;
        }
    }
}

void kernel(std::vector<BitMatrix64>& M, BitMatrix64& kerM)
{
    if (M.empty())
        throw "kernel: empty list of BitMatrix64";
    size_t m = M[0].rows();
    size_t n = 0;
    for (auto& r: M)
    {
        if (r.rows() != m) throw "kernel: incompatible BitMatrix64 objects";
        n += r.cols();
    }
    size_t N = M[0].cols();
    size_t r = 0;
    size_t k = 0;
    std::vector<int> c;
    c.resize(m, 0);
    for (size_t i = 0; i < m; i++) c[i] = -1;
    std::vector<int> d;
    d.resize(n, 0);

    while (k < n)
    {
        // 2. [Scan column]
        size_t j = 0;
        size_t u = k / N;
        //while (j < m && (!BitOperations64::bitSet(k % N, M[u].row_[j]) || c[j] != -1)) ++j;
        while (j < m && (c[j] != -1 || !BitOperations64::bitSet(k % N, M[u].row_[j]))) ++j;
        if (j >= m)
        {
            ++r;
            d[k] = -1;
        }
        if (j < m)
        {
            // 3. [Eliminate]
            int dd = 0;
            for (size_t i = 0; i < m; ++i)
            {
                if (i != j)
                {
                    size_t u = k / N;
                    dd = BitOperations64::bitSet(k % N, M[u].row_[i]);
                    BitOperations64::clearBit(k % N, M[u].row_[i]);
                    if (dd)
                    {
                        for (auto& r: M)
                        {
                            r.row_[i] ^= r.row_[j];
                        }
                    }
                }
            }
            c[j] = static_cast<int>(k);
            d[k] = static_cast<int>(j);
        }
        // 4. [Finished?]
        ++k;
    }
    // 5. [Output kernel]
    kerM.row_.resize(n);
    for (size_t i = 0; i < n; ++i) kerM.row_[i] = 0UL;
    const size_t max_columns = 16;
    kerM.cols_ = max_columns;
    if (r < max_columns)
        kerM.cols_ = r;

    size_t col = 0;
    for (k = 0; k < n; ++k)
    {
        if (d[k] == -1)
        {
            for (size_t i = 0; i < n; ++i)
            {
                if (d[i] >= 0)
                {
                    int bit = 0;
                    size_t u = k / N;
                    bit = BitOperations64::bitSet(k % N, M[u].row_[d[i]]);
                    BitOperations64::copyBit(bit, col, kerM.row_[i]);
                }
                else if (i == k)
                {
                    BitOperations64::setBit(col, kerM.row_[i]);
                }
                else
                {
                    BitOperations64::clearBit(col, kerM.row_[i]);
                }
            }
            ++col;
            if (col >= kerM.cols_)
                break;
        }
    }
}

void BitMatrix64::read(std::istream& is)
{
    std::string str;
    if (getline(is, str))
    {
        long int rows = std::atol(str.c_str());
        row_.resize(rows);
    }
    if (getline(is, str))
    {
        cols_ = std::atoi(str.c_str());
    }
    is >> *this;
}

void BitMatrix64::read(MemoryMappedFile& is)
{
    std::string str;
    if (getline(is, str))
    {
        long int rows = std::atol(str.c_str());
        row_.resize(rows);
    }
    if (getline(is, str))
    {
        cols_ = std::atoi(str.c_str());
    }
    is >> *this;
}

void BitMatrix64::readTranspose(std::istream& is)
{
    // reads a BitMatrix64 from a file which has been written with printTranspose
    std::string str;
    cols_ = BitOperations64::BITS_IN_WORD;

    if (!getline(is, str))
    {
        return;
    }
    size_t rows = str.size();
    row_.resize(rows);
    size_t col = 0;
    bool all_zeros = true;
    for (size_t i = 0; i < rows; ++i)
    {
        if (str[i] == '1')
        {
            BitOperations64::setBit(col, row_[i]);
            all_zeros = false;
        }
    }
    if (!all_zeros)
    {
        ++col;
    }

    while (col < cols_ && getline(is, str))
    {
        all_zeros = true;
        for (size_t i = 0; i < rows; ++i)
        {
            if (str[i] == '1')
            {
                BitOperations64::setBit(col, row_[i]);
                all_zeros = false;
            }
        }
        if (!all_zeros)
        {
            ++col;
        }
    }

    cols_ = col;
}

void BitMatrix64::readTranspose(MemoryMappedFile& is)
{
    // reads a BitMatrix64 from a file which has been written with printTranspose
    std::string str;
    cols_ = BitOperations64::BITS_IN_WORD;

    if (!getline(is, str))
    {
        return;
    }
    size_t rows = str.size();
    row_.resize(rows);
    size_t col = 0;
    bool all_zeros = true;
    for (size_t i = 0; i < rows; ++i)
    {
        if (str[i] == '1')
        {
            BitOperations64::setBit(col, row_[i]);
            all_zeros = false;
        }
    }
    if (!all_zeros)
    {
        ++col;
    }

    while (col < cols_ && getline(is, str))
    {
        all_zeros = true;
        for (size_t i = 0; i < rows; ++i)
        {
            if (str[i] == '1')
            {
                BitOperations64::setBit(col, row_[i]);
                all_zeros = false;
            }
        }
        if (!all_zeros)
        {
            ++col;
        }
    }

    cols_ = col;
}

void BitMatrix64::write(std::ostream& os) const
{
    os << static_cast<unsigned int>(rows()) << std::endl;
    os << static_cast<unsigned int>(cols()) << std::endl;
    os << *this;
}

std::istream& operator>>(std::istream& istr, BitMatrix64& bm)
{
    // assume we have already allocated space for rows
    std::string str;
    size_t i = 0;
    while (i < bm.rows() && getline(istr, str))
    {
        // each line consists of up to BITS_IN_WORDS 0s and 1s
        for (size_t j = 0; j < bm.cols(); j++)
        {
            int bit = str.c_str()[j] - '0';
            if (bit) BitOperations64::setBit(j, bm.row_[i]);
            else BitOperations64::clearBit(j, bm.row_[i]);
        }

        ++i;
    }

    return istr;
}

MemoryMappedFile& operator>>(MemoryMappedFile& istr, BitMatrix64& bm)
{
    // assume we have already allocated space for rows
    std::string str;
    size_t i = 0;
    while (i < bm.rows() && getline(istr, str))
    {
        // each line consists of up to BITS_IN_WORDS 0s and 1s
        for (size_t j = 0; j < bm.cols(); j++)
        {
            int bit = str.c_str()[j] - '0';
            if (bit) BitOperations64::setBit(j, bm.row_[i]);
            else BitOperations64::clearBit(j, bm.row_[i]);
        }

        ++i;
    }

    return istr;
}

std::ostream& operator<<(std::ostream& os, const BitMatrix64& bm)
{
    if (bm.cols() == 0) return os;
    for (size_t i = 0; i < bm.rows(); i++)
    {
        for (size_t j = 0; j < bm.cols(); j++)
        {
            if (BitOperations64::bitSet(j, bm.row_[i])) os << 1;
            else os << 0;
        }
        os << std::endl;
    }
    return os;
}

bool SparseMatrix3::parse(const std::string& str, size_t row)
{
    long int num_cols = SparseRow::begin_parse__(str);

    double density = (double)num_cols / (double)rows_;

    // First, write out very dense rows to be processed later
    if (very_dense_count_ < max_very_dense_rows && density > very_dense_density_bound)
    {
        write_very_dense_row(str);
        return false;
    }

    if (density > medium_density_bound)
    {
        add_to_medium_dense_rows(num_cols);
        ++medium_count_;
        return true;
    }

    sparse_->add_row(sparse_count_, num_cols);
    if (sparse_->cols() > cols_)
    {
        cols_ = sparse_->cols();
    }
    ++sparse_count_;

    return true;
}

bool SparseMatrix3::parse_for_sizing(const std::string& str, long int row)
{
    long int num_cols = SparseRow::begin_parse__(str);

    double density = (double)num_cols / (double)rows_;

    if (very_dense_count_ < max_very_dense_rows && density > very_dense_density_bound)
    {
        ++very_dense_count_;
        return false;
    }

    if (density > medium_density_bound)
    {
        add_to_size_of_medium_dense_rows(num_cols);
        ++medium_count_;
        return true;
    }

    sparse_allocated_points_ += num_cols;

    return true;
}

void SparseMatrix3::extend_dense(size_t stripe)
{
    if (stripe >= medium_.size())
    {
        medium_.reserve(stripe + 1);
        for (size_t i = medium_.size(); i < stripe + 1; ++i)
        {
            SparseMatrix4* smp = new SparseMatrix4;
            medium_.push_back(smp);
        }
    }
}

void SparseMatrix3::add_to_medium_dense_rows(long int num_cols)
{
    std::ostringstream oss;
    long int prev_stripe = -1L;
    long int cols_in_stripe = 0L;
    std::unordered_map<size_t, std::string> striped_row;
    while (char* s = strtok(0, " "))
    {
        size_t col = std::atol(s);
        size_t stripe = col / stripe_size;
        extend_dense(stripe);

        if (static_cast<long int>(stripe) != prev_stripe)
        {
            if (prev_stripe != -1)
            {
                std::string str = oss.str();
                oss.str("");
                //std::cerr << "cols_in_stripe = " << cols_in_stripe << std::endl;
                oss << cols_in_stripe << str;
                striped_row[prev_stripe] = oss.str();
                oss.str("");
            }
            prev_stripe = static_cast<long int>(stripe);
            cols_in_stripe = 0L;
        }

        oss << " " << static_cast<unsigned int>(col);
        ++cols_in_stripe;

        //std::cerr << "(stripe, row, col) = (" << stripe << ", " << row << ", " << col << ")" << std::endl;
        if (col > cols_)
        {
            cols_ = col;
        }
    }
    std::string str = oss.str();
    oss.str("");
    //std::cerr << "cols_in_stripe = " << cols_in_stripe << std::endl;
    oss << cols_in_stripe << str;
    striped_row[prev_stripe] = oss.str();
    extend_dense(prev_stripe);
    for (size_t stripe = 0; stripe < medium_.size(); ++stripe)
    {
        auto found = striped_row.find(stripe);
        if (found != striped_row.end())
        {
            medium_[stripe]->add_row(0, found->second);
        }
        else
        {
            medium_[stripe]->add_row(0, "0");
        }
    }
}

void SparseMatrix3::add_to_size_of_medium_dense_rows(long int num_cols)
{
    long int prev_stripe = -1L;
    long int cols_in_stripe = 0L;
    while (char* s = strtok(0, " "))
    {
        size_t col = std::atol(s);
        size_t stripe = col / stripe_size;

        if (static_cast<long int>(stripe) != prev_stripe)
        {
            if (prev_stripe != -1)
            {
                // reserve 1 point to hold the column count
                stripe_allocated_points_[prev_stripe] += cols_in_stripe + 1;
                if (prev_stripe > static_cast<long int>(number_of_stripes_))
                    number_of_stripes_ = prev_stripe;
            }
            prev_stripe = static_cast<long int>(stripe);
            cols_in_stripe = 0L;
        }

        ++cols_in_stripe;

    }
    stripe_allocated_points_[prev_stripe] += cols_in_stripe + 1;
    if (prev_stripe > static_cast<long int>(number_of_stripes_))
        number_of_stripes_ = prev_stripe;
}

SparseMatrix3::SparseMatrix3(const std::string& file, bool split)
    : rows_(0), cols_(0), sparse_(0), sparse_count_(0), sparse_allocated_points_(0), medium_(0), medium_count_(0), number_of_stripes_(0), very_dense_file_(0), very_dense_count_(0)
{
    MemoryMappedFile mmf(file.c_str());
    {
        // Pass 1, to work out how much memory we need to allocate
        cols_ = 0;
        long int rows = 0;
        std::string str;
        if (!getline(mmf, str))
            return;
        rows_ = std::atol(str.c_str());

        while (getline(mmf, str))
        {
            if (parse_for_sizing(str, rows))
            {
                ++rows;
            }
        }
        rows_ = rows;
        std::cerr << "About to create SparseMatrix2 : sparse_allocated_points_ = " << static_cast<unsigned int>(sparse_allocated_points_) << std::endl;
        sparse_ = new SparseMatrix2(sparse_allocated_points_);

        ++number_of_stripes_;
        medium_.reserve(number_of_stripes_);

        for (size_t i = 0; i < number_of_stripes_; ++i)
        {
            SparseMatrix4* smp = new SparseMatrix4(medium_count_, stripe_allocated_points_[i]);
            medium_.push_back(smp);
        }
        medium_count_ = 0;
        very_dense_count_ = 0;
    }

    mmf.reset();

    cols_ = 0;
    long int rows = 0;
    std::string str;
    if (!getline(mmf, str))
        return;
    rows_ = std::atol(str.c_str());

    while (getline(mmf, str))
    {
        if (parse(str, rows))
        {
            ++rows;
        }
    }
    rows_ = rows;
    ++cols_;

    for (size_t i = 0; i < medium_.size(); ++i)
    {
        SparseMatrix4* smp = medium_[i];
        smp->last_point_ = smp->next_point_;
        smp->next_point_ = smp->set_points_;
    }

    if (very_dense_file_)
    {
        delete very_dense_file_;
    }

    std::cerr << "SparseMatrix3::SparseMatrix3() : rows_ = " << static_cast<unsigned int>(rows_) << ", cols_ = " << static_cast<unsigned int>(cols_) << ", sparse_count_ = " << static_cast<unsigned int>(sparse_count_) << ", medium_count_ = " << static_cast<unsigned int>(medium_count_) << ", very_dense_count_ = " << static_cast<unsigned int>(very_dense_count_) << std::endl;
}

SparseMatrix3::~SparseMatrix3()
{
    clear();
}

void SparseMatrix3::clear()
{
    rows_ = 0;
    cols_ = 0;
    if (sparse_)
    {
        delete sparse_;
        sparse_ = 0;
    }

    for (size_t i = 0; i < medium_.size(); ++i)
    {
        delete medium_[i];
    }
    medium_.clear();
}

void multiply(const SparseMatrix3& A, const BitMatrix& X, BitMatrix& AX)
{
    if (A.cols() != X.rows()) throw "multiply: Incompatible SparseMatrix3 and BitMatrix";
    AX.row_.clear();
    AX.row_.resize(A.rows());
    AX.cols_ = X.cols_;

    // 1. Multiply the sparse rows
    BitMatrixRowIterator AX_row_iter = AX.row_.begin();
    uint32_t resultRow = 0UL;
    for (SparseMatrix2::Point* p = A.sparse_->set_points_; p != A.sparse_->last_point_; ++p)
    {
        long int col = *p;
        if (col < 0)
        {
            col = -col - 1;
            resultRow ^= *(X.row_.vec_ + col);
            *AX_row_iter = resultRow;
            ++AX_row_iter;
            resultRow = 0UL;
        }
        else
        {
            resultRow ^= *(X.row_.vec_ + col);
        }
    }

    BitMatrixRowIterator AX_row_iter_1 = AX_row_iter;
    // 2. Multiply the medium rows which are striped
    for (size_t i = 0; i < A.medium_.size(); ++i)
    {
        AX_row_iter = AX_row_iter_1;
        SparseMatrix4* sm4 = A.medium_[i];
        SparseMatrix4::Point* p = sm4->set_points_;
        while (p != sm4->last_point_)
        {
            long int col = *p;
            size_t col_count = -(col + 1);
            ++p;
            uint32_t resultRow = 0UL;
            for (size_t k = 0; k < col_count; ++k, ++p)
            {
                long int col = *p;
                resultRow ^= *(X.row_.vec_ + col);
            }
            *AX_row_iter ^= resultRow;
            ++AX_row_iter;
        }
    }
}

void multiply(const SparseMatrix3& A, const BitMatrix64& X, BitMatrix64& AX)
{
    if (A.cols() != X.rows()) throw "multiply: Incompatible SparseMatrix3 and BitMatrix64";
    AX.row_.clear();
    AX.row_.resize(A.rows());
    AX.cols_ = X.cols_;

    // 1. Multiply the sparse rows
    BitMatrix64RowIterator AX_row_iter = AX.row_.begin();
    unsigned long long int resultRow = 0UL;
    for (SparseMatrix2::Point* p = A.sparse_->set_points_; p != A.sparse_->last_point_; ++p)
    {
        long int col = *p;
        if (col < 0)
        {
            col = -col - 1;
            resultRow ^= *(X.row_.vec_ + col);
            *AX_row_iter = resultRow;
            ++AX_row_iter;
            resultRow = 0UL;
        }
        else
        {
            resultRow ^= *(X.row_.vec_ + col);
        }
    }

    BitMatrix64RowIterator AX_row_iter_1 = AX_row_iter;
    // 2. Multiply the medium rows which are striped
    for (size_t i = 0; i < A.medium_.size(); ++i)
    {
        AX_row_iter = AX_row_iter_1;
        SparseMatrix4* sm4 = A.medium_[i];
        SparseMatrix4::Point* p = sm4->set_points_;
        while (p != sm4->last_point_)
        {
            long int col = *p;
            size_t col_count = -(col + 1);
            ++p;
            unsigned long long int resultRow = 0UL;
            for (size_t k = 0; k < col_count; ++k, ++p)
            {
                long int col = *p;
                resultRow ^= *(X.row_.vec_ + col);
            }
            *AX_row_iter ^= resultRow;
            ++AX_row_iter;
        }
    }
}

void multiplyt(const SparseMatrix3& A, const BitMatrix& X, BitMatrix& AtX)
{
    if (A.rows() != X.rows()) throw "multiplyt: Incompatible SparseMatrix3 and BitMatrix";
    AtX.row_.clear();
    AtX.row_.resize(A.cols());
    AtX.cols_ = X.cols_;

    // 1. Multiply the sparse rows
    BitMatrixRowIterator X_row_iter = X.row_.begin();
    for (SparseMatrix2::Point* p = A.sparse_->set_points_;
            p != A.sparse_->last_point_;
            ++p)
    {
        long int col = *p;
        if (col < 0)
        {
            col = -col - 1;
            AtX.row_[col] ^= *X_row_iter;
            ++X_row_iter;
        }
        else
        {
            AtX.row_[col] ^= *X_row_iter;
        }
    }

    // 2. Multiply the medium rows which are striped
    BitMatrixRowIterator X_row_iter_1 = X_row_iter;
    for (size_t i = 0; i < A.medium_.size(); ++i)
    {
        X_row_iter = X_row_iter_1;
        SparseMatrix4* sm4 = A.medium_[i];
        SparseMatrix4::Point* p = sm4->set_points_;
        while (p != sm4->last_point_)
        {
            long int num_cols = - *p - 1;
            ++p;
#if 1
            uint32_t X_row = *X_row_iter;
            SparseMatrix4::Point* last_point = p + num_cols;
            while (p != last_point)
            {
                long int col = *p;
                AtX.row_[col] ^= X_row;
                ++p;
            }
#else
            for (size_t j = 0; j < num_cols; ++j, ++p)
            {
                long int col = *p;
                AtX.row_[col] ^= *X_row_iter;
            }
#endif
            ++X_row_iter;
        }
    }
}

void multiplyt(const SparseMatrix3& A, const BitMatrix64& X, BitMatrix64& AtX)
{
    if (A.rows() != X.rows()) throw "multiplyt: Incompatible SparseMatrix3 and BitMatrix64";
    AtX.row_.clear();
    AtX.row_.resize(A.cols());
    AtX.cols_ = X.cols_;

    // 1. Multiply the sparse rows
    BitMatrix64RowIterator X_row_iter = X.row_.begin();
    for (SparseMatrix2::Point* p = A.sparse_->set_points_;
            p != A.sparse_->last_point_;
            ++p)
    {
        long int col = *p;
        if (col < 0)
        {
            col = -col - 1;
            AtX.row_[col] ^= *X_row_iter;
            ++X_row_iter;
        }
        else
        {
            AtX.row_[col] ^= *X_row_iter;
        }
    }

    // 2. Multiply the medium rows which are striped
    BitMatrix64RowIterator X_row_iter_1 = X_row_iter;
    for (size_t i = 0; i < A.medium_.size(); ++i)
    {
        X_row_iter = X_row_iter_1;
        SparseMatrix4* sm4 = A.medium_[i];
        SparseMatrix4::Point* p = sm4->set_points_;
        while (p != sm4->last_point_)
        {
            long int num_cols = - *p - 1;
            ++p;
            for (long int j = 0; j < num_cols; ++j, ++p)
            {
                long int col = *p;
                AtX.row_[col] ^= *X_row_iter;
            }
            ++X_row_iter;
        }
    }
}

void sym_multiply(const SparseMatrix3& B, const BitMatrix& X, BitMatrix& AX)
{
    BitMatrix BX;
    multiply(B, X, BX);
    multiplyt(B, BX, AX);
}

void sym_multiply(const SparseMatrix3& B, const BitMatrix64& X, BitMatrix64& AX)
{
    BitMatrix64 BX;
    multiply(B, X, BX);
    multiplyt(B, BX, AX);
}

std::ostream& operator<<(std::ostream& os, const SparseMatrix3& sm)
{
    os << static_cast<unsigned int>(sm.rows()) << std::endl;
    os << *sm.sparse_;

    if (!sm.medium_.empty())
    {
        // For each row
        //std::cerr << "sm.medium_.size() = " << sm.medium_.size() << std::endl;

        std::vector<SparseMatrix4::Point*> p(sm.medium_.size());
        for (size_t j = 0; j < sm.medium_.size(); ++j)
        {
            p[j] = sm.medium_[j]->set_points_;
        }
        for (size_t i = 0; i < sm.medium_[0]->rows(); ++i)
        {
            long int cols = 0;
            std::stringstream ss;
            for (size_t j = 0; j < sm.medium_.size(); ++j)
            {
                long int col_count = *(p[j]);
                if (col_count >= 0)
                {
                    std::cerr << "Problem: col_count = " << col_count << std::endl;
                }
                ++(p[j]);
                for (size_t k = 0; k < static_cast<size_t>(-(col_count + 1)); ++k, ++(p[j]))
                {
                    long int col = *(p[j]);
                    ss << " " << col;
                    ++cols;
                }
            }
            os << cols << ss.str() << std::endl;
        }
    }

    if (sm.very_dense_count_)
    {
        MemoryMappedFile mmf("dense_rows.txt");
        std::string str;
        while (getline(mmf, str))
        {
            os << str << std::endl;
        }
    }

    return os;
}

void SparseMatrix3::write_very_dense_row(const std::string& str)
{
    if (!very_dense_file_)
    {
        very_dense_file_ = new std::fstream("dense_rows.txt", std::ios::out);
    }
    *very_dense_file_ << str << std::endl;
    ++very_dense_count_;
}

void SparseMatrix3::multiply_dense_part_by_bit_matrix(const BitMatrix& L, const BitMatrix& R, BitMatrix& BL, BitMatrix& BR) const
{
    if (!very_dense_count_)
        return;

    BL.row_.clear();
    BL.row_.resize(very_dense_count_);
    BL.cols_ = L.cols();
    BR.row_.clear();
    BR.row_.resize(very_dense_count_);
    BR.cols_ = R.cols();

    MemoryMappedFile mmf("dense_rows.txt");
    std::string str;

    BitMatrixRowIterator BL_row_iter = BL.row_.begin();
    BitMatrixRowIterator BR_row_iter = BR.row_.begin();

    while (getline(mmf, str))
    {
        SparseRow::begin_parse__(str);
        uint32_t resultRowL = 0UL;
        uint32_t resultRowR = 0UL;

        while (char* s = strtok(0, " "))
        {
            long int n = std::atol(s);
            resultRowL ^= *(L.row_.vec_ + n);
            resultRowR ^= *(R.row_.vec_ + n);
        }
        *BL_row_iter = resultRowL;
        *BR_row_iter = resultRowR;
        ++BL_row_iter;
        ++BR_row_iter;
    }
}

void SparseMatrix3::multiply_dense_part_by_bit_matrix(const BitMatrix64& L, const BitMatrix64& R, BitMatrix64& BL, BitMatrix64& BR) const
{
    if (!very_dense_count_)
        return;

    BL.row_.clear();
    BL.row_.resize(very_dense_count_);
    BL.cols_ = L.cols();
    BR.row_.clear();
    BR.row_.resize(very_dense_count_);
    BR.cols_ = R.cols();

    MemoryMappedFile mmf("dense_rows.txt");
    std::string str;

    BitMatrix64RowIterator BL_row_iter = BL.row_.begin();
    BitMatrix64RowIterator BR_row_iter = BR.row_.begin();

    while (getline(mmf, str))
    {
        SparseRow::begin_parse__(str);
        unsigned long long int resultRowL = 0UL;
        unsigned long long int resultRowR = 0UL;

        while (char* s = strtok(0, " "))
        {
            long int n = std::atol(s);
            resultRowL ^= *(L.row_.vec_ + n);
            resultRowR ^= *(R.row_.vec_ + n);
        }
        *BL_row_iter = resultRowL;
        *BR_row_iter = resultRowR;
        ++BL_row_iter;
        ++BR_row_iter;
    }
}

void SparseMatrix3::multiply_dense_part_by_bit_matrix(const BitMatrix& L, BitMatrix& BL) const
{
    if (!very_dense_count_)
        return;

    BL.row_.clear();
    BL.row_.resize(very_dense_count_);
    BL.cols_ = L.cols();

    MemoryMappedFile mmf("dense_rows.txt");
    std::string str;

    BitMatrixRowIterator BL_row_iter = BL.row_.begin();

    while (getline(mmf, str))
    {
        SparseRow::begin_parse__(str);
        uint32_t resultRow = 0UL;
        while (char* s = strtok(0, " "))
        {
            long int n = std::atol(s);
            resultRow ^= *(L.row_.vec_ + n);
        }
        *BL_row_iter = resultRow;
        ++BL_row_iter;
    }
}

void SparseMatrix3::multiply_dense_part_by_bit_matrix(const BitMatrix64& L, BitMatrix64& BL) const
{
    if (!very_dense_count_)
        return;

    BL.row_.clear();
    BL.row_.resize(very_dense_count_);
    BL.cols_ = L.cols();

    MemoryMappedFile mmf("dense_rows.txt");
    std::string str;

    BitMatrix64RowIterator BL_row_iter = BL.row_.begin();

    while (getline(mmf, str))
    {
        SparseRow::begin_parse__(str);
        unsigned long long int resultRow = 0UL;

        while (char* s = strtok(0, " "))
        {
            long int n = std::atol(s);
            resultRow ^= *(L.row_.vec_ + n);
        }
        *BL_row_iter = resultRow;
        ++BL_row_iter;
    }
}

// --- SparseMatrix4 ----
SparseMatrix4::SparseMatrix4() : dense_file_(0), split_(false), rows_(0L), rows_added_(0L), cols_(0L), allocated_points_(0), dense_rows_(0L)
    , set_points_(0), next_point_(0), last_point_(0)
{
}

SparseMatrix4::SparseMatrix4(size_t rows, size_t allocated_points) : dense_file_(0), split_(false), rows_(rows), rows_added_(0L), cols_(0L), allocated_points_(allocated_points), dense_rows_(0L)
    , set_points_(0), next_point_(0), last_point_(0)
{
    std::cerr << "SparseMatrix4::SparseMatrix4 : allocated_points_ = " << static_cast<unsigned int>(allocated_points_) << std::endl;
    try
    {
        set_points_ = new Point [ allocated_points_ ];
        next_point_ = set_points_;
        last_point_ = set_points_ + allocated_points_;
    }
    catch (std::bad_alloc& )
    {
        std::ostringstream os;
        os << "SparseMatrix4::SparseMatrix4 : Caught std::bad_alloc : allocated_points_ = " << static_cast<unsigned int>(allocated_points_);
        std::cerr << os.str() << std::endl;
        throw os.str();
    }
    catch (...)
    {
        std::ostringstream os;
        os << "SparseMatrix4::SparseMatrix4 : Caught unknown exception : allocated_points_ = " << static_cast<unsigned int>(allocated_points_);
        std::cerr << os.str() << std::endl;
        throw os.str();
    }
}

SparseMatrix4::SparseMatrix4(const std::string& file, bool split) : dense_file_(0), split_(split), rows_(0L), rows_added_(0L), cols_(0L), allocated_points_(0), dense_rows_(0L)
    , set_points_(0), next_point_(0), last_point_(0)
{
#ifndef USING_GCC
    allocated_points_ = 1024UL * 14UL * 1024UL;
#else
    allocated_points_ = 1024UL * 145UL * 1024UL;
#endif

    set_points_ = new Point [ allocated_points_ ];
    next_point_ = set_points_;
    last_point_ = set_points_ + allocated_points_;

    MemoryMappedFile mmf(file.c_str());
    cols_ = 0;
    long int rows = 0;
    std::string str;
    if (!getline(mmf, str))
        return;
    rows_ = std::atol(str.c_str());

    while (getline(mmf, str))
    {
        if (parse(str, rows))
        {
            ++rows;
        }
    }
    rows_ = rows;
    ++cols_;
    last_point_ = next_point_;
    next_point_ = set_points_;
    if (dense_file_)
    {
        delete dense_file_;
        dense_file_ = 0;
    }
}

SparseMatrix4::SparseMatrix4(const std::vector<long int>& points, long int rows, long int columns)
{
    set_points_ = new Point [ points.size() ];
    for (size_t i = 0; i < points.size(); ++i)
    {
        set_points_[i] = points[i];
    }
    last_point_ = set_points_ + points.size();
    rows_ = rows;
    rows_added_ = rows;
    cols_ = columns;
}

SparseMatrix4::~SparseMatrix4()
{
    clear();
}

void multiply(const SparseMatrix4& A, const BitMatrix& X, BitMatrix& AX)
{
    //Timer.start("multiply");
    if (A.cols() != X.rows()) throw "multiply: Incompatible SparseMatrix4 and BitMatrix";
    AX.row_.clear();
    AX.row_.resize(A.rows());
    AX.cols_ = X.cols_;
    uint32_t resultRow = 0UL;
    BitMatrixRowIterator AX_row_iter = AX.row_.begin();
    for (SparseMatrix4::Point* p = A.set_points_; p != A.last_point_; ++p)
    {
        long int col = *p;
        if (col < 0)
        {
            if (p != A.set_points_)
            {
                *AX_row_iter = resultRow;
                ++AX_row_iter;
                resultRow = 0UL;
            }
        }
        else
        {
            resultRow ^= *(X.row_.vec_ + col);
        }
    }
    *AX_row_iter = resultRow;
    //Timer.stop();
}

void multiply(const SparseMatrix4& A, const BitMatrix64& X, BitMatrix64& AX)
{
    //Timer.start("multiply");
    if (A.cols() != X.rows()) throw "multiply: Incompatible SparseMatrix4 and BitMatrix64";
    AX.row_.clear();
    AX.row_.resize(A.rows());
    AX.cols_ = X.cols_;
    unsigned long long int resultRow = 0UL;
    BitMatrix64RowIterator AX_row_iter = AX.row_.begin();
    for (SparseMatrix4::Point* p = A.set_points_; p != A.last_point_; ++p)
    {
        long int col = *p;
        if (col < 0)
        {
            if (p != A.set_points_)
            {
                *AX_row_iter = resultRow;
                ++AX_row_iter;
                resultRow = 0UL;
            }
        }
        else
        {
            resultRow ^= *(X.row_.vec_ + col);
        }
    }
    *AX_row_iter = resultRow;
    //Timer.stop();
}

void multiplyt(const SparseMatrix4& A, const BitMatrix& X, BitMatrix& AtX)
{
    //
    // Calculate A^tX = ((X^t) A)^t which is a BitMatrix having A.cols() rows and X.cols() columns
    // but without explicitly calculating A^t
    //
    //            m
    //   t       ---               t  t         t
    // (X A)   = \   X   A     = (X A)      = (A X)
    //      ji   /    kj  ki           ij          ij
    //           ---
    //           k=1
    //
    // where m = X.rows() = A.rows()
    // and all elements are in GF(2)
    //
    if (A.rows() != X.rows()) throw "multiplyt: Incompatible SparseMatrix4 and BitMatrix";
    AtX.row_.clear();
    AtX.row_.resize(A.cols());
    AtX.cols_ = X.cols_;
    BitMatrixRowIterator X_row_iter = X.row_.begin();
    for (SparseMatrix4::Point* p = A.set_points_; p != A.last_point_; ++p)
    {
        long int col = *p;
        if (col < 0)
        {
            if (p != A.set_points_)
            {
                ++X_row_iter;
            }
        }
        else
        {
            AtX.row_[col] ^= *X_row_iter;
        }
    }
}

void multiplyt(const SparseMatrix4& A, const BitMatrix64& X, BitMatrix64& AtX)
{
    //
    // Calculate A^tX = ((X^t) A)^t which is a BitMatrix having A.cols() rows and X.cols() columns
    // but without explicitly calculating A^t
    //
    //            m
    //   t       ---               t  t         t
    // (X A)   = \   X   A     = (X A)      = (A X)
    //      ji   /    kj  ki           ij          ij
    //           ---
    //           k=1
    //
    // where m = X.rows() = A.rows()
    // and all elements are in GF(2)
    //
    if (A.rows() != X.rows()) throw "multiplyt: Incompatible SparseMatrix4 and BitMatrix64";
    AtX.row_.clear();
    AtX.row_.resize(A.cols());
    AtX.cols_ = X.cols_;
    BitMatrix64RowIterator X_row_iter = X.row_.begin();
    for (SparseMatrix4::Point* p = A.set_points_; p != A.last_point_; ++p)
    {
        long int col = *p;
        if (col < 0)
        {
            if (p != A.set_points_)
            {
                ++X_row_iter;
            }
        }
        else
        {
            AtX.row_[col] ^= *X_row_iter;
        }
    }
}

void sym_multiply(const SparseMatrix4& B, const BitMatrix& X, BitMatrix& AX)
{
    BitMatrix BX;
    multiply(B, X, BX);
    multiplyt(B, BX, AX);
}

void sym_multiply(const SparseMatrix4& B, const BitMatrix64& X, BitMatrix64& AX)
{
    BitMatrix64 BX;
    multiply(B, X, BX);
    multiplyt(B, BX, AX);
}

void SparseMatrix4::write_dense_row(const std::string& str)
{
    if (!dense_file_)
    {
        dense_file_ = new std::fstream("dense_rows.txt", std::ios::out);
    }
    *dense_file_ << str << std::endl;
}

void SparseMatrix4::add_row(size_t row, size_t num_cols)
{
    if (next_point_ + num_cols + 1 > last_point_)
    {
        //std::cerr << "Extending: num_cols = " << num_cols << std::endl;
        extend(num_cols + 1);
    }
    *next_point_ = -static_cast<long int>(num_cols + 1);
    ++next_point_;
    while (char* s = strtok(0, " "))
    {
        long int n = std::atol(s);
        *next_point_ = n;
        ++next_point_;
        if (n > cols_)
        {
            cols_ = n;
        }
    }
    ++rows_added_;
}

void SparseMatrix4::add_row(size_t row, const std::string& str)
{
    long int num_cols = SparseRow::begin_parse__(str);

    add_row(row, num_cols);
}

bool SparseMatrix4::parse(const std::string& str, size_t row)
{
    long int num_cols = SparseRow::begin_parse__(str);

    if (split_ && dense_rows_ < max_very_dense_rows && (double)num_cols / (double)rows_ > very_dense_density_bound)
    {
        write_dense_row(str);
        ++dense_rows_;
        return false;
    }
    if (next_point_ + num_cols + 1 > last_point_)
    {
        //std::cerr << "Extending: num_cols = " << num_cols << std::endl;
        extend(num_cols + 1);
    }

    *next_point_ = -(num_cols + 1);
    ++next_point_;
    while (char* s = strtok(0, " "))
    {
        long int n = std::atol(s);
        *next_point_ = n;
        ++next_point_;
        if (n > cols_)
        {
            cols_ = n;
        }
    }
    ++rows_added_;
    return true;
}

void SparseMatrix4::extend(long int cols)
{
    size_t rows_left = rows_ - rows_added_;
    if (rows_ < rows_added_)
    {
        rows_left = 100UL;
    }
    //std::cerr << "SparseMatrix4::extend() : row = " << row << ", cols = " << cols << ", rows_left = " << rows_left << std::endl;
    //std::cerr << "allocated_points_ = " << allocated_points_ << std::endl;
    allocated_points_ += std::max<long int>(rows_left * 20, cols + rows_left);
    //std::cerr << "allocated_points_ = " << allocated_points_ << std::endl;
    size_t offset = next_point_ - set_points_;
    Point* new_set_points = new Point [ allocated_points_ ];
    next_point_ = new_set_points + offset;
    last_point_ = new_set_points + allocated_points_;
    memcpy(new_set_points, set_points_, sizeof(Point)*offset);
    delete [] set_points_;
    set_points_ = new_set_points;
}

void SparseMatrix4::clear()
{
    delete [] set_points_;
}

void SparseMatrix4::multiply_dense_part_by_bit_matrix(const BitMatrix& L, const BitMatrix& R, BitMatrix& BL, BitMatrix& BR) const
{
    if (!split_)
        return;

    delete [] set_points_;
    set_points_ = 0;

    BL.row_.clear();
    BL.row_.resize(dense_rows_);
    BL.cols_ = L.cols();
    BR.row_.clear();
    BR.row_.resize(dense_rows_);
    BR.cols_ = R.cols();

    MemoryMappedFile mmf("dense_rows.txt");
    std::string str;

    BitMatrixRowIterator BL_row_iter = BL.row_.begin();
    BitMatrixRowIterator BR_row_iter = BR.row_.begin();

    while (getline(mmf, str))
    {
        SparseRow::begin_parse__(str);
        uint32_t resultRowL = 0UL;
        uint32_t resultRowR = 0UL;

        while (char* s = strtok(0, " "))
        {
            long int n = std::atol(s);
            resultRowL ^= *(L.row_.vec_ + n);
            resultRowR ^= *(R.row_.vec_ + n);
        }
        *BL_row_iter = resultRowL;
        *BR_row_iter = resultRowR;
        ++BL_row_iter;
        ++BR_row_iter;
    }
}

std::ostream& operator<<(std::ostream& os, const SparseMatrix4& sm)
{
    os << static_cast<unsigned int>(sm.rows()) << std::endl;
    for (SparseMatrix4::Point* p = sm.set_points_; p != sm.last_point_; ++p)
    {
        long int col = *p;
        if (col < 0)
        {
            long int num_cols = -col - 1;
            if (p != sm.set_points_)
            {
                os << std::endl;
            }
            os << num_cols;
        }
        else
        {
            os << " " << col;
        }
    }
    os << std::endl;
    return os;
}

void SparseMatrix4::multiply_dense_part_by_bit_matrix(const BitMatrix& L, BitMatrix& BL) const
{
    if (!split_)
        return;

    delete [] set_points_;
    set_points_ = 0;

    BL.row_.clear();
    BL.row_.resize(dense_rows_);
    BL.cols_ = L.cols();

    MemoryMappedFile mmf("dense_rows.txt");
    std::string str;

    BitMatrixRowIterator BL_row_iter = BL.row_.begin();

    while (getline(mmf, str))
    {
        SparseRow::begin_parse__(str);
        uint32_t resultRow = 0UL;

        while (char* s = strtok(0, " "))
        {
            long int n = std::atol(s);
            resultRow ^= *(L.row_.vec_ + n);
        }
        *BL_row_iter = resultRow;
        ++BL_row_iter;
    }
}

