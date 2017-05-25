#ifndef __NEWMATRIX_H
#define __NEWMATRIX_H
#include <vector>
#include <sstream>
#include "Polynomial.h"

template <class X> void combine(X& x, const X& v, const X& q)
{
    // replace x by x - q * v
    x -= q * v;
}

template <class X> X divide(const X& a, const X& b)
{
    const X zero(0L);
    const X one(1L);
    X r = a / b;
    if (b * r != a &&
            ((a < zero && b > zero) ||
             (a > zero && b < zero)))
    {
        r -= one;
    }
    return r;
}
template <class X>
class Matrix
{
public:
    Matrix() : rows_(0), columns_(0), data_(0) {}
    Matrix(size_t rows, size_t columns) : rows_(rows), columns_(columns)
    {
        const X zero(0L);
        data_ = new X [ rows_ * columns_ ];
        X* p = data_;
        for (size_t i = 0; i < rows_; ++i)
        {
            for (size_t j = 0; j < columns_; ++j)
            {
                *p = zero;
                ++p;
            }
        }
    }
    Matrix(const Matrix& m) : rows_(m.rows()), columns_(m.columns())
    {
        data_ = new X [ rows_ * columns_ ];
        X* p = data_;
        X* q = m.data_;
        for (size_t i = 0; i < rows_; ++i)
        {
            for (size_t j = 0; j < columns_; ++j)
            {
                *p = *q;
                ++p;
                ++q;
            }
        }
    }
    Matrix(size_t rank, const X& value, size_t dummy) : rows_(rank), columns_(rank)
    {
        const X zero(0L);
        data_ = new X [ rows_ * columns_ ];
        X* p = data_;
        for (size_t i = 0; i < rows_; ++i)
        {
            for (size_t j = 0; j < columns_; ++j)
            {
                if (i == j) *p = value;
                else *p = zero;
                ++p;
            }
        }
    }
    ~Matrix()
    {
        delete [] data_;
    }
    void set_size(size_t rows, size_t columns)
    {
        const X zero(0L);
        size_t old_size = rows_ * columns_;
        size_t new_size = rows * columns;
        if (old_size < new_size)
        {
            delete [] data_;
            data_ = new X [ new_size ];
        }
        rows_ = rows;
        columns_ = columns;
        X* p = data_;
        for (size_t i = 0; i < rows_; ++i)
        {
            for (size_t j = 0; j < columns_; ++j)
            {
                *p = zero;
                ++p;
            }
        }
    }
    void add_column()
    {
        const X zero(0L);
        size_t new_size = rows_ * (columns_ + 1);
        X* newdata = new X [ new_size ];
        X* p = newdata;
        X* q = data_;
        for (size_t i = 0; i < rows_; ++i)
        {
            for (size_t j = 0; j < columns_; ++j)
            {
                *p = *q;
                ++p;
                ++q;
            }
            *p = zero;
            ++p;
        }
        ++columns_;
        delete [] data_;
        data_ = newdata;
    }
    Matrix& operator=(const Matrix& m)
    {
        if (&m == this)
            return *this;
        size_t new_size = m.rows() * m.columns();
        size_t old_size = rows() * columns();
        rows_ = m.rows();
        columns_ = m.columns();
        if (old_size < new_size)
        {
            delete [] data_;
            data_ = new X [ new_size ];
        }
        X* p = data_;
        X* q = m.data_;
        for (size_t i = 0; i < rows_; ++i)
        {
            for (size_t j = 0; j < columns_; ++j)
            {
                *p = *q;
                ++p;
                ++q;
            }
        }
        return *this;
    }
    friend Matrix operator+(const Matrix& m1, const Matrix& m2)
    {
        if (m1.rows() != m2.rows() ||
                m1.columns() != m2.columns())
        {
            throw std::string("operator+(Matrix, Matrix) : inconsistent sizes");
        }
        Matrix m = m1;
        X* p = m.data_;
        X* q = m2.data_;
        for (size_t i = 0; i < m.rows(); ++i)
        {
            for (size_t j = 0; j < m.columns(); ++j)
            {
                *p += *q;
                ++p;
                ++q;
            }
        }
        return m;
    }
    friend Matrix operator-(const Matrix& m1, const Matrix& m2)
    {
        if (m1.rows() != m2.rows() ||
                m1.columns() != m2.columns())
        {
            throw std::string("operator-(Matrix, Matrix) : inconsistent sizes");
        }
        Matrix m = m1;
        X* p = m.data_;
        X* q = m2.data_;
        for (size_t i = 0; i < m.rows(); ++i)
        {
            for (size_t j = 0; j < m.columns(); ++j)
            {
                *p -= *q;
                ++p;
                ++q;
            }
        }
        return m;
    }
    friend Matrix operator*(const X& x, const Matrix& m)
    {
        Matrix mm(m);
        X* p = mm.data_;
        for (size_t i = 0; i < m.rows(); ++i)
        {
            for (size_t j = 0; j < m.columns(); ++j)
            {
                *p *= x;
                ++p;
            }
        }
        return mm;
    }
    friend Matrix operator*(const Matrix& m, const X& x)
    {
        Matrix mm(m);
        X* p = mm.data_;
        for (size_t i = 0; i < m.rows(); ++i)
        {
            for (size_t j = 0; j < m.columns(); ++j)
            {
                *p *= x;
                ++p;
            }
        }
        return mm;
    }
    Matrix& operator/=(const X& x)
    {
        X* p = data_;
        for (size_t i = 0; i < rows(); ++i)
        {
            for (size_t j = 0; j < columns(); ++j)
            {
                *p /= x;
                ++p;
            }
        }
        return *this;
    }
    friend Matrix operator/(const Matrix& m, const X& x)
    {
        Matrix mm(m);
        X* p = mm.data_;
        for (size_t i = 0; i < m.rows(); ++i)
        {
            for (size_t j = 0; j < m.columns(); ++j)
            {
                *p /= x;
                ++p;
            }
        }
        return mm;
    }
    friend std::vector<X> operator*(const Matrix& m1, const std::vector<X>& v)
    {
        const X zero(0L);
        if (m1.columns() != v.size())
        {
            throw std::string("operator*(Matrix, vector) : inconsistent sizes");
        }
        std::vector<X> m;
        m.resize(v.size());
        X* p = m1.data_;
        for (size_t i = 0; i < m1.rows(); ++i)
        {
            m[i] = zero;
            for (size_t j = 0; j < m1.columns(); ++j)
            {
                m[i] += *p * v[j];
                ++p;
            }
        }
        return m;
    }
    template <class Y>
    //friend Matrix<X> operator* (const Matrix<X>& m1, const Matrix<Y>& m2);
    Matrix<X> operator* (const Matrix<Y>& m2) const
    {
        const Matrix<X>& m1(*this);
        const X zero(0L);
        if (m1.columns() != m2.rows())
        {
            throw std::string("operator*(Matrix, Matrix) : inconsistent sizes");
        }
        Matrix<X> m(m1.rows(), m2.columns());
        X* p = m.data_;
        for (size_t i = 0; i < m1.rows(); ++i)
        {
            for (size_t j = 0; j < m2.columns(); ++j)
            {
                *p = zero;
                for (size_t k = 0; k < m1.columns(); ++k)
                {
                    *p += m1(i, k) * m2(k, j);
                }
                ++p;
            }
        }
        return m;
    }

    friend std::ostream& operator<< (std::ostream& os, const Matrix& m)
    {
        X* p = m.data_;
        for (size_t i = 0; i < m.rows(); ++i)
        {
            for (size_t j = 0; j < m.columns(); ++j)
            {
                os << *p << " ";
                ++p;
            }
            os << std::endl;
        }
        return os;
    }
    friend std::istream& operator>> (std::istream& is, Matrix& m);
    Matrix& operator*=(const X& x)
    {
        X* p = data_;
        for (size_t i = 0; i < rows(); ++i)
        {
            for (size_t j = 0; j < columns(); ++j)
            {
                *p *= x;
                ++p;
            }
        }
        return *this;
    }
    template <class Y>
    Matrix<X>& operator*=(const Matrix<Y>& m)
    {
        const X zero(0L);
        if (columns() != m.rows())
        {
            throw std::string("Matrix::operator*=(Matrix) : inconsistent sizes");
        }
        Matrix mm(rows(), m.columns());
        X* p = mm.data_;
        for (size_t i = 0; i < rows(); ++i)
        {
            for (size_t j = 0; j < m.columns(); ++j)
            {
                *p = zero;
                for (size_t k = 0; k < columns(); ++k)
                {
                    *p += this->operator()(i, k) * m(k, j);
                }
                ++p;
            }
        }
        *this = mm;
        return *this;
    }
    Matrix& operator*=(const Matrix& m)
    {
        const X zero(0L);
        if (columns() != m.rows())
        {
            throw std::string("Matrix::operator*=(Matrix) : inconsistent sizes");
        }
        Matrix mm(rows(), m.columns());
        X* p = mm.data_;
        for (size_t i = 0; i < rows(); ++i)
        {
            for (size_t j = 0; j < m.columns(); ++j)
            {
                *p = zero;
                for (size_t k = 0; k < columns(); ++k)
                {
                    *p += this->operator()(i, k) * m(k, j);
                }
                ++p;
            }
        }
        *this = mm;
        return *this;
    }
    size_t rows() const
    {
        return rows_;
    }
    size_t columns() const
    {
        return columns_;
    }
    X trace() const
    {
        const X zero(0L);
        X tr(zero);
        size_t r = rows();
        size_t c = columns();
        size_t n = r;
        if (c < n) n = c;
        for (size_t i = 0; i < n; ++i)
        {
            tr += this->operator()(i, i);
        }
        return tr;
    }
    Matrix<X > transpose()
    {
        Matrix trans(columns(), rows());
        X* p = data_;
        for (size_t i = 0; i < rows(); ++i)
        {
            for (size_t j = 0; j < columns(); ++j)
            {
                trans(j,i) = *p;
                ++p;
            }
        }
        return trans;
    }
    bool operator==(const Matrix& m) const
    {
        if (rows() != m.rows()) return false;
        if (columns() != m.columns()) return false;
        X* p = data_;
        X* q = m.data_;
        for (size_t i = 0; i < rows(); ++i)
        {
            for (size_t j = 0; j < columns(); ++j)
            {
                if (*p != *q) return false;
                ++p;
                ++q;
            }
        }
        return true;
    }
    bool operator!=(const Matrix& m) const
    {
        return !(*this == m);
    }
    X& operator()(size_t i, size_t j)
    {
        return data_[i * columns_ + j];
    }
    const X& operator()(size_t i, size_t j) const
    {
        return data_[i * columns_ + j];
    }
    X& at(size_t i, size_t j)
    {
        if (i >= rows() || j >= columns())
        {
            throw std::string("Matrix::at() : out of range");
        }
        return operator()(i, j);
    }
    const X& at(size_t i, size_t j) const
    {
        if (i >= rows() || j >= columns())
        {
            throw std::string("Matrix::at() : out of range");
        }
        return operator()(i, j);
    }
    // swap rows r1 and r2
    void swap(size_t r1, size_t r2)
    {
        X* p1 = data_ + r1 * columns();
        X* p2 = data_ + r2 * columns();
        for (size_t j = 0; j < columns(); ++j)
        {
            std::swap(*p1, *p2);
            ++p1;
            ++p2;
        }
    }
    // swap columns c1 and c2
    void swap_columns(size_t c1, size_t c2)
    {
        X* p1 = data_ + c1;
        X* p2 = data_ + c2;
        for (size_t j = 0; j < rows(); ++j)
        {
            std::swap(*p1, *p2);
            p1 += columns();
            p2 += columns();
        }
    }
    void transform_row(size_t r1, size_t r2, const X& x, const X& y, const X& u, const X& v)
    {
        X t1;
        X t2;
        X t3;
        X t4;
        for (size_t i = 0; i < columns(); i++)
        {
            t1 = x * this->operator()(r1,i);
            t2 = y * this->operator()(r2,i);
            t1 += t2;
            t3 = u * this->operator()(r1,i);
            t4 = v * this->operator()(r2,i);
            t3 += t4;
            this->operator()(r1,i) = t1;
            this->operator()(r2,i) = t3;
        }
    }
    void transform_column(size_t c1, size_t c2, const X& x, const X& y, const X& u, const X& v)
    {
        X t1;
        X t2;
        X t3;
        X t4;
        for (size_t i = 0; i < rows(); i++)
        {
            t1 = x * this->operator()(i,c1);
            t2 = y * this->operator()(i,c2);
            t1 += t2;
            t3 = u * this->operator()(i,c1);
            t4 = v * this->operator()(i,c2);
            t3 += t4;
            this->operator()(i,c1) = t1;
            this->operator()(i,c2) = t3;
        }
    }
    X dot(int r1, int r2)
    {
        X result = 0L;
        for (size_t i = 0; i < columns(); i++)
        {
            result += this->operator()(r1,i) * this->operator()(r2,i);
        }
        return result;
    }
    X dot_on_columns(int c1, int c2)
    {
        X result = 0L;
        for (size_t i = 0; i < rows(); i++)
        {
            result += this->operator()(i, c1)* this->operator()(i, c2);
        }
        return result;
    }
    void combine_columns(size_t c1, size_t c2, const X& q)
    {
        // replace column c1 by column c1 - q * column c2
        X* p = data_;
        for (size_t r = 0; r < rows(); ++r)
        {
            //*(p + c1) -= q * (*(p + c2));
            combine(*(p + c1), *(p + c2), q);
            p += columns_;
        }
    }
    void combine_columns(size_t c1, size_t c2, const X& u, const X& v)
    {
        // replace column c1 by u * column c1 - v * column c2
        X* p1 = data_ + c1;
        X* p2 = data_ + c2;
        for (size_t r = 0; r < rows(); ++r)
        {
            //return data_[i * columns_ + j];
            //(*this)(r-1,c1) -= q * (*this)(r-1,c2);
            //data_[(r-1)*columns_ + c1] -= q * data_[(r-1)*columns_ + c2];
            *p1 *= u;
            *p1 -= v * (*p2);
            p1 += columns_;
            p2 += columns_;
        }
    }
    void combine_columns_mod_R(size_t c1, size_t c2, const X& u, const X& v, const X& R)
    {
        // replace column c1 by (u * column c1 - v * column c2) mod R
        X* p1 = data_ + c1;
        X* p2 = data_ + c2;
        for (size_t r = 0; r < rows(); ++r)
        {
            //return data_[i * columns_ + j];
            //(*this)(r-1,c1) -= q * (*this)(r-1,c2);
            //data_[(r-1)*columns_ + c1] -= q * data_[(r-1)*columns_ + c2];
            *p1 *= u;
            *p1 -= v * (*p2);
            *p1 %= R;
            p1 += columns_;
            p2 += columns_;
        }
    }
    void negate_column(size_t c)
    {
        X* p = data_;
        for (size_t r = 0; r < rows(); ++r)
        {
            *(p + c) = - *(p + c);
            p += columns_;
        }
    }
    friend void copy_column(const Matrix& source, size_t source_column,
                            Matrix& target, size_t target_column)
    {
        X* source_ptr = source.data_;
        X* target_ptr = target.data_;
        for (size_t r = 0; r < source.rows(); ++r)
        {
            *(target_ptr + target_column) = *(source_ptr + source_column);
            source_ptr += source.columns();
            target_ptr += target.columns();
        }

    }
private:
    size_t rows_;
    size_t columns_;
    X* data_;
};

// Algorithm 2.2.1 (Square Linear System using Gaussian Elimination)
template <class X> bool solve(const Matrix<X >& MM, const std::vector<X >& BB, std::vector<X >& XX)
{
    // only valid for n x n matrix M
    const X zero(0L);
    const X one(1L);
    if (MM.rows() != MM.columns())
    {
        return false;
    }
    int n = MM.rows();
    static Matrix<X > M(1,1);
    M = MM;
    std::vector<X > B = BB;
    int j = 0;
    while (j < n)
    {
        // step 3. [ Find non-zero entry ]
        int i = j;
        while (i < n && M(i, j) == zero) i++;
        if (i >= n)
        {
            throw std::string("M is not invertible");
        }
        //cout << "First non-zero entry = " << i << endl;
        // step 4. [Swap?]
        if (i > j)
        {
            for (int l = j; l < n; l++)
            {
                X tmp = M(j, l);
                M(j, l) = M(i, l);
                M(i, l) = tmp;
            }
            X tmp = B[j];
            B[j] = B[i];
            B[i] = tmp;
        }
        //cout << "M after Swap:" << endl;
        //cout << M;

        // step 5 [Eliminate]
        static std::vector<X > C;
        C.resize(n, zero);
        X d = one / M(j, j);
        //cout << "d = " << d << endl;
        int k = j + 1;
        for (k = j + 1; k < n; k++)
        {
            C[k] = d * M(k, j);
        }
        for (k = j + 1; k < n; k++)
        {
            M(k, j) = zero;
            for (int l = j + 1; l < n; l++)
            {
                M(k, l) = M(k, l) - C[k] * M(j, l);
            }
        }
        for (k = j + 1; k < n; k++)
        {
            B[k] = B[k] - C[k] * B[j];
        }
        //cout << "M after elimination:" << endl;
        //cout << M;
        // step 2. [Finished?]
        j++;
    }

    // step 6. [Solve triangular system]
    int i = n - 1;
    XX.resize(n);
    for (i = n - 1; i >= 0; --i)
    {
        XX[i] = B[i];
        for (int j = i + 1; j < n; j++)
        {
            XX[i] -= M(i, j) * XX[j];
        }
        XX[i] /= M(i, i);
    }

    return true;
}

// Algorithm 2.2.2 (Inverse of a Matrix)
template <class X> bool invert(const Matrix<X >& MM, Matrix<X >& XX)
{
    // only valid for n x n matrix M
    const X zero(0L);
    const X one(1L);
    if (MM.rows() != MM.columns())
    {
        return false;
    }
    int n = MM.rows();
    static Matrix<X > M(1,1);
    M = MM;
    int j = 0;
    Matrix<X > B(n, one, 0); // B is n x n identity matrix
    //cout << "B:" << endl;
    //cout << B;

    while (j < n)
    {
        // step 3. [ Find non-zero entry ]
        int i = j;
        while (i < n && M(i, j) == zero) i++;
        if (i >= n)
        {
            throw std::string("M is not invertible");
        }
        //cout << "First non-zero entry = " << i << endl;
        // step 4. [Swap?]
        if (i > j)
        {
            for (int l = j; l < n; l++)
            {
                X tmp = M(j,l);
                M(j,l) = M(i,l);
                M(i,l) = tmp;
            }
            for (int m = 0; m < n; m++)
            {
                X tmp = B(j,m);
                B(j,m) = B(i,m);
                B(i,m) = tmp;
            }
        }

        // step 5 [Eliminate]
        static std::vector<X > C;
        C.resize(n, zero);
        X d = one / M(j,j);
        //cout << "d = " << d << endl;
        int k = j + 1;
        for (k = j + 1; k < n; k++)
        {
            C[k] = d * M(k,j);
        }
        //cout << "C:" << endl;
        //cout << C;
        for (k = j + 1; k < n; k++)
        {
            M(k,j) = zero;
            for (int l = j + 1; l < n; l++)
            {
                M(k,l) -= C[k] * M(j,l);
            }
        }
        //cout << "M:" << endl;
        //cout << M;
        for (k = j + 1; k < n; k++)
        {
            //cout << "k = " << k << endl;
            for (int l = 0; l < n; l++)
            {
                //cout << "l = " << l << " " << flush;
                //cout << "B[" << k << "," << l << "] = " << B(k,l) << endl;
                //cout << "B[" << j << "," << l << "] = " << B(j,l) << endl;
                //cout << "C[" << k << "] = " << C[k] << endl;
                B(k,l) -= C[k] * B(j,l);
            }
            //cout << endl;
        }
        //cout << "B after elimination:" << endl;
        //cout << B;
        // step 2. [Finished?]

        j++;
    }

    // step 6. [Solve triangular system]
    //cout << "solve triangular system, M:" << endl;
    //cout << M;
    int i = n - 1;
    for (i = n - 1; i >= 0; --i)
    {
        int l = 0;
        for (l = 0; l < n; l++)
        {
            XX(i,l) = B(i,l);
        }
        for (int j = i + 1; j < n; j++)
        {
            for (int l = 0; l < n; l++)
            {
                XX(i,l) -= M(i,j) * XX(j,l);
            }
        }
        for (l = 0; l < n; l++)
        {
            XX(i,l) /= M(i,i);
        }
    }

    return true;
}

// Algorithm 2.2.3 (Determinant, Using Ordinary Elimination)
template <class X> X determinant(const Matrix<X >& MM)
{
    // only valid for n x n matrix M
    const X zero(0L);
    const X one(1L);
    if (MM.rows() != MM.columns())
    {
        return zero;
    }
    int n = MM.rows();
    static Matrix<X > M(1,1);
    M = MM;
    int j = 0;
    X x = one;
    while (j < n)
    {
        // step 3. [ Find non-zero entry ]
        int i = j;
        while (i < n && M(i,j) == zero) i++;
        if (i >= n)
        {
            return zero;
        }
        //cout << "First non-zero entry = " << i << endl;
        // step 4. [Swap?]
        if (i > j)
        {
            for (int l = j; l < n; l++)
            {
                X tmp = M(j,l);
                M(j,l) = M(i,l);
                M(i,l) = tmp;
            }
            x = zero - x;
        }

        // step 5 [Eliminate]
        //cout << "Step 5. [Eliminate]" << endl;
        static std::vector<X > C;
        C.resize(n, zero);
        X d = one / M(j,j);
        //cout << "d = " << d << endl;
        int k = j + 1;
        for (k = j + 1; k < n; k++)
        {
            C[k] = d * M(k,j);
        }
        //cout << "C:" << endl;
        //cout << C;
        for (k = j + 1; k < n; k++)
        {
            M(k,j) = zero;
            for (int l = j + 1; l < n; l++)
            {
                M(k,l) -= C[k] * M(j,l);
            }
        }
        //cout << "M after elimination:" << endl;
        //cout << M;
        x = x * M(j,j);
        // step 2. [Finished?]

        j++;
    }

    return x;
}

// Algorithm 2.2.6 (Determinant Using Gauss-Bareiss)
// MM has coefficients in an integral domain R
template <class X> X determinant_in_integral_domain(const Matrix<X >& MM)
{
    // only valid for n x n matrix M
    const X zero(0L);
    const X one(1L);
    const X minusone(-1L);
    if (MM.rows() != MM.columns())
    {
        return zero;
    }
    int n = MM.rows();
    static Matrix<X > M(1,1);
    M = MM;
    // 1. [Initialize]
    int k = 0;
    X c = one;
    X s = one;

    while (k < n - 1)
    {
        X p = M(k,k);
        // 3. [Is p = 0?]
        if (p == zero)
        {
            // look for first non-zero coefficient M(i,k)
            int i = k;
            while (i < n && M(i,k) == zero) i++;
            if (i >= n)
            {
                return zero;
            }
            //cout << "First non-zero entry = " << i << endl;

            for (int j = k; j < n; j++)
            {
                X tmp = M(i,j);
                M(i,j) = M(k,j);
                M(k,j) = tmp;
            }
            //cout << "M after swapping :" << endl;
            //cout << M;

            s = s * minusone;
//cout << "1" << endl;
//cout << "k = " << k << endl;
//cout << "M(k,k) = " << M(k,k) << endl;
            p = M(k,k);
//cout << "2" << endl;
        }

        //cout << "about to do step 4, p =" << p << endl;
        // 4. [Main step]
        for (int i = k + 1; i < n; i++)
        {
            for (int j = k + 1; j < n; j++)
            {
                //cout << "(i,j,k) = (" << i << "," << j << "," << k << ")" << endl;
                //cout << "t = " << p << "*" << M(i,j) << " - " << M(i,k) << "*" << M(k,j) << endl;
                X t = p * M(i,j) - M(i,k) * M(k,j);
                //cout << "t = " << t << endl;
                //cout << "c = " << c << endl;
                M(i,j) = t / c;
            }
            //cout << "M after main step :" << endl;
            //cout << M;
        }
        c = p;
        k++;
    }

    return s * M(n-1,n-1);
}

// Algorithm 2.2.7 (Characteristic Polynomial and Adjoint Matrix)
template <class X>
Polynomial<X > characteristic_polynomial(const Matrix<X >& MM, Matrix<X >& MM_adj)
{
    // only valid for n x n matrix M
    const X zero(0L);
    const X one(1L);
    const X minusone(-1L);
    if (MM.rows() != MM.columns())
    {
        return zero;
    }
    int n = MM.rows();
    static Matrix<X > M(1,1);
    M = MM;

    // 1. [Initialize]
    int i = 0;
    Matrix<X > I_n(n, one, 0); // n x n Identity matrix
    Matrix<X > C = I_n;
    static std::vector<X > a; // coefficients of characteristic polynomial
    a.resize(n + 1, zero);

    a[n] = one;
    //cout << "a[0] = " << a[n] << endl;

    //cout << "M: " << endl << M;
    //cout << "C: " << endl << C;
    //cout << "MC :" << endl << M * C;
    while (i < n - 1)
    {
        C = M * C;
        a[n - i - 1] = minusone * (C.trace() / X((long)(i + 1)));
        //cout << "a[" << i + 1 << "] = " << a[n - i - 1] << endl;
        //cout << "I_n = " << endl << I_n;
        C = C + a[n - i - 1] * I_n;
        //cout << "M: " << endl << M;
        //cout << "C: " << endl << C;
        //cout << "MC :" << endl << M * C;
        i++;
    }

    a[0] = minusone * ((M * C).trace() / X((long)n));
    //cout << "a[n] = " << a[0] << endl;
    MM_adj = C;
    if (n % 2 == 0)
    {
        MM_adj = minusone * MM_adj;
    }

    Polynomial<X > P(a);
    return P;
}

// Algorithm 2.2.9 (Hessenberg)
template <class X> Polynomial<X > characteristic_polynomial(const Matrix<X >& MM)
{
    // only valid for n x n matrix M
    const X zero(0L);
    const X one(1L);
    const X minusone(-1L);
    if (MM.rows() != MM.columns())
    {
        return zero;
    }
    int n = MM.rows();
    static Matrix<X > M(1,1);
    M = MM;

    // 1. [Initialize]
    static Matrix<X > H(1,1);
    H = MM;
    int m = 1;
    X t;

    while (m < n - 1)
    {
        // 2. [Search for non-zero]
        int i = m + 1;
        while (i < n && H(i, m - 1) == zero) i++;
        if (i < n)
        {
            if (H(m, m - 1) != zero) i = m;
            t = H(i, m - 1);
            if (i > m)
            {
                for (int j = m - 1; j < n; j++)
                {
                    X tmp = H(i,j);
                    H(i,j) = H(m,j);
                    H(m,j) = tmp;
                }
                for (int r = 0; r < n; r++)
                {
                    X tmp = H(r,i);
                    H(r,i) = H(r,m);
                    H(r,m) = tmp;
                }
            }
            // 3. [Eliminate]
            for (i = m + 1; i < n; i++)
            {
                if (H(i, m - 1) != zero)
                {
                    X u = H(i, m - 1) / t;
                    for (int j = m; j < n; j++)
                    {
                        H(i,j) = H(i,j) - u * H(m,j);
                    }
                    H(i, m - 1) = zero;
                    for (int r = 0; r < n; r++)
                    {
                        H(r,m) = H(r,m) + u * H(r,i);
                    }
                }
            }
        }

        // 4. [Hessenberg finished?]
        m++;
    }

    //cout << "Hessenberg form: " << endl << H;

    // 5. [Initialize characteristic polynomial]
    static std::vector<Polynomial<X > > p;
    p.resize(n + 1, zero);
    p[0] = one;
    m = 1;

    while (m <= n)
    {
        // 6. [Initialize computation]
        static std::vector<X > pp;
        pp.resize(2, zero);
        pp[1] = one;
        pp[0] = H(m - 1, m - 1) * minusone;
        p[m] = p[m - 1] * Polynomial<X >(pp);
        t = one;

        // 7. [Compute p_m]
        for (int i = 1; i <= m - 1; i++)
        {
            t = t * H(m - i, m - i - 1);
            p[m] = p[m] - t * H(m - i - 1, m - 1) * p[m - i - 1];
        }

        // 8. [Finished?]
        m++;
    }

    return p[n];
}

// Algorithm 2.3.1 (Kernel of a Matrix)
template <class X> Matrix<X > kernel(const Matrix<X >& MM)
{
    const X zero(0L);
    const X one(1L);
    const X minusone(-1L);
    int m = MM.rows();
    int n = MM.columns();
//   cout << "kernel: m = " << m << ", n = " << n << endl;
//cout << "0" << endl;
    static Matrix<X > M(1,1);
    M = MM;
//cout << "1" << endl;
    // 1. [Initialize]
    int r = 0;
    int k = 0;
    static std::vector<int > c;
    c.resize(m, 0);
    for (int i = 0; i < m; i++) c[i] = -1;
    static std::vector<int > d;
    d.resize(n, 0);

    while (k < n)
    {
//cout << "2, k = " << k << endl;
        // 2. [Scan column]
        int j = 0;
        while (j < m && (M(j,k) == zero || c[j] != -1)) j++;
        if (j >= m)
        {
            r++;
            d[k] = -1;
        }
        if (j < m)
        {
//cout << "3, j = " << j << endl;
            // 3. [Eliminate]
            X dd = minusone / M(j,k);
            M(j,k) = minusone;
            for (int s = k + 1; s < n; s++)
            {
                M(j,s) = dd * M(j,s);
            }
            for (int i = 0; i < m; i++)
            {
                if (i != j)
                {
                    dd = M(i,k);
                    M(i,k) = zero;
                    for (int s = k + 1; s < n; s++)
                    {
                        M(i,s) += dd * M(j,s);
                    }
                }
            }
            c[j] = k;
            d[k] = j;
        }

        // 4. [Finished?]
        k++;
//cout << "M = " << endl << M;
    }

//cout << "5" << endl;
    // 5. [Output kernel]

//cout << "r = " << r << endl;
    Matrix <X > ker(n, r);
    int col = 0;
    for (k = 0; k < n; k++)
    {
//cout << "6, k = " << k << endl;
        if (d[k] == -1)
        {
            for (int i = 0; i < n; i++)
            {
//cout << "7, i = " << i << ", d[i] = " << d[i] << endl;
                if (d[i] >= 0) ker(i, col) = M(d[i], k);
                else if (i == k) ker(i, col) = one;
                else ker(i, col) = zero;
            }
            col++;
        }
    }

    return ker;
}

// Algorithm 2.3.2 (Image of a Matrix)
template <class X> Matrix<X > image(const Matrix<X >& MM)
{
    const X zero(0L);
    const X minusone(-1L);
    int m = MM.rows();
    int n = MM.columns();
    static Matrix<X > M(1,1);
    M = MM;
    static Matrix<X > N(1,1);
    N = M;
//cout << "1" << endl;
    // 1. [Initialize]
    int r = 0;
    int k = 0;
    static std::vector<int > c;
    c.resize(m, 0);
    for (int i = 0; i < m; i++) c[i] = -1;
    static std::vector<int > d;
    d.resize(n, 0);

    while (k < n)
    {
//cout << "2, k = " << k << endl;
        // 2. [Scan column]
        int j = 0;
        while (j < m && (M(j,k) == zero || c[j] != -1)) j++;
        if (j >= m)
        {
            r++;
            d[k] = -1;
        }
        if (j < m)
        {
//cout << "3, j = " << j << endl;
            // 3. [Eliminate]
            X dd = minusone / M(j,k);
            M(j,k) = minusone;
            for (int s = k + 1; s < n; s++)
            {
                M(j,s) = dd * M(j,s);
            }
            for (int i = 0; i < m; i++)
            {
                if (i != j)
                {
                    dd = M(i,k);
                    M(i,k) = zero;
                    for (int s = k + 1; s < n; s++)
                    {
                        M(i,s) += dd * M(j,s);
                    }
                }
            }
            c[j] = k;
            d[k] = j;
        }

        // 4. [Finished?]
        k++;
//cout << "M = " << endl << M;
    }

//cout << "5" << endl;
    // 5. [Output image]

//cout << "r = " << r << endl;
    Matrix <X > im(m, n - r);
    int col = 0;
    for (int j = 0; j < m; j++)
    {
//cout << "6, j = " << j << endl;
        if (c[j] != -1)
        {
            for (int i = 0; i < m; i++)
            {
                im(i, col) = N(i, c[j]);
            }
            col++;
        }
    }

    return im;
}

// Algorithm 2.3.4 (Inverse Image)
template <class X> Matrix<X > inverse_image(const Matrix<X >& MM, const std::vector<X >& BB)
{
    const X zero(0L);
    const X minusone(-1L);
    size_t m = MM.rows();
    size_t n = MM.columns();
    Matrix<X > XX(n, 1);
    if (BB.size() != m)
    {
        return XX;
    }
    // 1. [Compute Kernel]
//   cout << "1. Compute Kernel" << endl;
    Matrix<X > M1(m, n+1);
    for (size_t row = 0; row < m; row++)
    {
        for (size_t col = 0; col < n; col++)
        {
            M1(row, col) = MM(row, col);
        }
        M1(row, n) = BB[row];
    }
    Matrix<X > V = kernel(M1);
    int r = V.columns();
//   cout << "kernel = " << endl << V;
//   cout << "r = " << r << endl;

    // 2. [Solution exists?]
//   cout << "2. Solution exists?" << endl;
    int j = 0;
    while (j < r && V(n, j) == zero) j++;
    if (j >= r)
    {
        throw std::string("inverse_image: No solution");
    }

    X d = minusone / V(n, j);

    // 3. [Output solution]
//   cout << "3. Output solution?" << endl;
    for (size_t i = 0; i < n; i++)
    {
//      cout << "i " << endl;
        XX(i, 0) = d * V(i, j);
    }

    return XX;
}

// Algorithm 2.3.5 (Inverse Image Matrix)
template <class X> Matrix<X > inverse_image_matrix(const Matrix<X >& MM, const Matrix<X >& VV)
{
    const X zero(0L);
    const X one(1L);
    int m = MM.rows();
    int n = MM.columns();
    int r = VV.columns();
    if (m != static_cast<int>(VV.rows()) || n > m)
    {
        std::ostringstream oss;
        oss << "inverse_image_matrix: M must be m x n with n <= m and V must be m x r, ";
        oss << "but M is " << m << " x " << n << " and V is " << VV.rows() << " x " << r;
        throw oss.str();
    }
    if (r == 0) return VV;
    static Matrix<X > M(1,1);
    M = MM;
    int j = 0;
    static Matrix<X > B(1,1);
    B = VV;
    Matrix<X > XX(n, r);

    while (j < n)
    {
        // step 3. [ Find non-zero entry ]
        //cout << "step 3, j = " << j << endl;
        int i = j;
        while (i < m && M(i, j) == zero) i++;
        if (i >= m)
        {
            throw std::string("columns of M are not linearly independent");
        }
        //cout << "First non-zero entry = " << i << endl;
        // step 4. [Swap?]
        if (i > j)
        {
            //cout << "step 4, i,j = " << i << "," << j << endl;
            for (int l = j; l < n; l++)
            {
                X tmp = M(j,l);
                M(j,l) = M(i,l);
                M(i,l) = tmp;
            }
            for (int m = 0; m < r; m++)
            {
                X tmp = B(j,m);
                B(j,m) = B(i,m);
                B(i,m) = tmp;
            }
        }

        // step 5 [Eliminate]
        //cout << "step 5" << endl;
        static std::vector<X > C;
        C.resize(m, zero);
        X d = one / M(j,j);
        //cout << "d = " << d << endl;
        int k = j + 1;
        for (k = j + 1; k < m; k++)
        {
            C[k] = d * M(k,j);
        }
        //cout << "C:" << endl;
        //cout << C;
        for (k = j + 1; k < m; k++)
        {
            M(k,j) = zero;
            for (int l = j + 1; l < n; l++)
            {
                M(k,l) -= C[k] * M(j,l);
            }
        }
        //cout << "M:" << endl;
        //cout << M;
        for (k = j + 1; k < m; k++)
        {
            //cout << "k = " << k << endl;
            for (int l = 0; l < r; l++)
            {
                //cout << "l = " << l << " " << flush;
                B(k,l) -= C[k] * B(j,l);
            }
            //cout << endl;
        }
        //cout << "B after elimination:" << endl;
        //cout << B;
        // step 2. [Finished?]

        j++;
    }

    // step 6. [Solve triangular system]
    //cout << "solve triangular system, M:" << endl;
    //cout << M;
    int i = n - 1;
    for (i = n - 1; i >= 0; --i)
    {
        int l = 0;
        for (l = 0; l < r; l++)
        {
            XX(i,l) = B(i,l);
        }
        for (int j = i + 1; j < n; j++)
        {
            for (int l = 0; l < r; l++)
            {
                XX(i,l) -= M(i,j) * XX(j,l);
            }
        }
        for (l = 0; l < r; l++)
        {
            XX(i,l) /= M(i,i);
        }
    }

    // step 7. [Check rest of matrix]
    //cout << "B = " << endl << B;
    int not_in_image = 0;
    for (int k = n; k < m; k++)
    {
        // calculate M'_k * X
        Matrix<X > M_k(1, n);
        for (int i = 0; i < n; i++) M_k(0,i) = M(k,i);
        //cout << "row " << k << " of M = " << endl << M_k;
        Matrix<X > B_k = M_k * XX;
        //cout << "M_k * XX = " << endl << B_k;
        for (int j = 0; j < r; j++)
        {
            if (B_k(0,j) != B(k,j))
            {
                not_in_image = 1;
            }
        }
    }
    if (not_in_image)
    {
        throw std::string("inverse_image_matrix: a column of V is not in the image of M");
    }

    return XX;
}

// Algorithm 2.3.6 (Supplement a Basis)
template <class X> Matrix<X > supplement(const Matrix<X >& MM)
{
    const X zero(0L);
    const X one(1L);
    static Matrix<X> M(1,1);
    M = MM;
    int n = M.rows();
    int k = M.columns();

    // Step 1. [Initialize]
    int s = 0;
    Matrix<X> B(n, one, 0);

    // Step 2. [Finished?]
    while (s < k)
    {
        // Step 3. [Search for non-zero]
        s++;
        //cout << "Step 3, B = " << endl << B;
        //cout << "M = " << endl << M;
        //cout << "B * M = " << endl << B * M;
        int t = s;
        while (t <= n && M(t-1, s-1) == zero) t++;
        if (t <= n)
        {
            X d = one / M(t-1, s-1);

            // Step 4. [Modify basis and eliminate]
            if (t != s)
            {
                for (int r = 0; r < n; r++)
                {
                    B(r, t-1) = B(r, s-1);
                }
            }
            for (int r = 0; r < n; r++)
            {
                B(r, s-1) = MM(r, s-1);
            }

            for (int j = s + 1; j <= k; j++)
            {
                if (t != s)
                {
                    X tmp = M(s-1,j-1);
                    M(s-1,j-1) = M(t-1,j-1);
                    M(t-1,j-1) = tmp;
                }
                M(s-1,j-1) = d * M(s-1,j-1);
                for (int i = 1; i <= n; i++)
                {
                    if (i != t && i != s)
                    {
                        M(i-1,j-1) -= M(i-1,s-1) * M(s-1,j-1);
                    }
                }
            }
            for (int i = 1; i <= n; i++)
            {
                M(i-1,s-1) = zero;
            }
            M(s-1,s-1) = one;
        }
        else
        {
            std::cout << "Supplement: rank of M = " << std::endl << M;
            std::cout << "is less than k = " << k << std::endl;
            return Matrix<X>(1,1);
        }
    }

    return B;
}

// Algorithm 2.3.7 (Supplement a Subspace in Another)
template <class X> Matrix<X > supplement_subspace(const Matrix<X >& VV, const Matrix<X >& MM)
{
    int r = VV.columns();
    int n = MM.columns();
    static Matrix<X> V(1,1);
    V = VV;
    static Matrix<X> M(1,1);
    M = MM;

    // Step 1. [Find new coordinates]
    Matrix<X> XX = inverse_image_matrix(M, V);
    if (XX.rows() == 1 && XX.columns() == 1)
    {
        throw std::string("F is not a subspace of E");
    }

    // Step 2. [Supplement X]
    Matrix<X> B = supplement(XX);

    // Step 3. [Supplement F in E]
    Matrix<X> C(n, n-r);
    for (int row = 0; row < n; row++)
    {
        for (int column = r; column < n; column++)
        {
            C(row, column - r) = B(row, column);
        }
    }

    return M * C;
}

// Algorithm 2.4.4 (Hermite Normal Form)
template <class X> Matrix<X > HNF(const Matrix<X >& AA)
{
    const X zero(0L);
    int m = AA.rows();
    int n = AA.columns();
    static Matrix <X > A(1,1);
    A = AA;
    int w_rows = m;
    int w_cols = n;
    if (m < n)
    {
        w_cols = m;
    }
    Matrix <X > W(w_rows, w_cols);

    // 1. [Initialize]
    int i = m;
    int k = n;
    int l = 1;
    if (m > n) l = m - n + 1;

    int finished = 0;
    while (!finished)
    {
        int done = 0;
        while (!done)
        {
            // 2. [Row Finished]
            //cout << "2. [Row Finished], (i, k) = (" << i << ", " << k << ")" << endl;
            int j = 1;
            while (j < k && A(i-1,j-1) == zero) j++;
            if (j >= k)
            {
                // replace column A_k by -A_k
                //std::cout << "i = " << i << ", k = " << k << ", A(i-1,k-1) = " << A(i-1,k-1) << std::endl;
                if (A(i-1,k-1) < zero)
                {
                    //std::cout << "HNF : i = " << i << ", k = " << k << std::endl;
                    for (int r = 1; r <= m; r++)
                    {
                        //A(r-1,k-1) *= minusone;
                        //std::cout << "r = " << r << "A(r-1,k-1) = " << A(r-1,k-1) << std::endl;
                        A(r-1,k-1) = -A(r-1,k-1);
                        //std::cout << "r = " << r << "A(r-1,k-1) = " << A(r-1,k-1) << std::endl;
                    }
                }
                //std::cout << "A = " << std::endl << A;
                done = 1;
            }
            else
            {
                // j is first non-zero entry
                // 3. [Choose non-zero entry]
                //cout << "3. [Choose non-zero entry]" << endl;
                int j0 = j;
                X min_val = A(i-1,j-1);
                //if (min_val < zero) min_val *= minusone;
                if (min_val < zero) min_val = -min_val;

                j++;
                while (j <= k)
                {
                    X val = A(i-1,j-1);
                    //if (val < zero) val *= minusone;
                    if (val < zero) val = -val;
                    if (val > zero && val < min_val)
                    {
                        min_val = val;
                        j0 = j;
                    }
                    j++;
                }
                //cout << "j0 = " << j0 << ", k = " << k << endl;
                if (j0 < k)
                {
                    // exchange column k with column j0
                    for (int r = 1; r <= m; r++)
                    {
                        X tmp = A(r-1,k-1);
                        A(r-1,k-1) = A(r-1, j0-1);
                        A(r-1, j0-1) = tmp;
                    }
                }
                if (A(i-1,k-1) < zero)
                {
                    for (int r = 1; r <= m; r++)
                    {
                        //A(r-1,k-1) *= minusone;
                        A(r-1,k-1) = -A(r-1,k-1);
                    }
                }
                //std::cout << "A = " << std::endl << A;
                X b = A(i-1,k-1);
                //cout << "b = " << b << endl;

                // 4. [Reduce]
                //cout << "4. [Reduce]" << endl;
                for (j = 1; j < k; j++)
                {
                    X q = A(i-1,j-1) / b;
                    for (int r = 1; r <= m; r++)
                    {
                        A(r-1,j-1) -= q * A(r-1,k-1);
                    }
                }
                //std::cout << "A = " << std::endl << A;
            }
        }

        // 5. [Final Reductions]
        //cout << "5. [Final Reductions]" << endl;
        X b = A(i-1,k-1);
        if (b == zero)
        {
            k++;
        }
        else
        {
            for (int j = k + 1; j <= n; j++)
            {
                X q = divide(A(i-1,j-1), b);
                for (int r = 1; r <= m; r++)
                {
                    A(r-1,j-1) -= q * A(r-1,k-1);
                }
            }
        }

        // 6. [Finished?]
        //cout << "6. [Finished?], i = " << i << endl;
        if (i == l)
        {
            for (int j = 1; j <= n - k + 1; j++)
            {
                for (int r = 1; r <= m; r++)
                {
                    W(r-1,j-1) = A(r-1, j+k-2);
                }
            }
            finished = 1;
        }
        else
        {
            i--;
            k--;
        }
    }
    return W;
}

// Algorithm 2.4.5 (Hermite Normal Form)
template <class X> Matrix<X > HNF1(const Matrix<X >& AA)
{
    const X zero(0L);

    Matrix<X > A(AA);

    int m = A.rows();
    int n = A.columns();
    int w_rows = m;
    int w_cols = n;
    if (m < n)
    {
        w_cols = m;
    }
    Matrix <X > W(w_rows, w_cols);
    std::vector<X > B;
    B.resize(m);

    // 1. [Initialize]
    int i = m;
    int j = n;
    int k = n;
    int l = 1;
    if (m > n) l = m - n + 1;

    bool done = false;
    while (!done)
    {
        // 2. [Check zero]
        while (j != 1)
        {
            --j;
            if (A(i-1,j-1) != zero)
            {
                // 3. [Euclidean step]
                // "Using Euclid's extended algorithm, compute (u,v,d) such that
                // u A(i-1,k-1) + v A(i-1,j-1) = d = gcd(A(i-1,k-1),A(i-1,j-1)),
                // with |u| and |v| minimal, in the sense that
                // -|a|/d < v sign(b) <= 0 and 1 <= u sign(a) <= |b|/d"
                X u;
                X v;
                X d = extended_gcd(A(i-1,k-1), A(i-1,j-1), u, v);
                X A_ik = A(i-1,k-1) / d;
                X A_ij = A(i-1,j-1) / d;

                for (int r = 1; r <= m; ++r)
                {
                    if (r == i) continue;
                    B[r-1] = u;
                    B[r-1] *= A(r-1,k-1);
                    B[r-1] += v * A(r-1,j-1);
                }

                A.combine_columns(j-1, k-1, A_ik, A_ij);
                A(i-1,k-1) = d;
                for (int r = 1; r <= m; ++r)
                {
                    if (r == i) continue;
                    A(r-1,k-1) = B[r-1];
                }
            }
        }

        // 4. [Final reductions]
        X b = A(i-1,k-1);
        if (b < zero)
        {
            A.negate_column(k-1);
            b = -b;
        }
        if (b == zero)
        {
            ++k;
        }
        else
        {
            for (int j = k + 1; j <= n; ++j)
            {
                X q = divide(A(i-1,j-1), b);

                A.combine_columns(j-1, k-1, q);
            }
        }

        // 5. [Finished?]
        if (i == l)
        {
            for (int j = 1; j <= n - k + 1; ++j)
            {
                copy_column(A, j+k-2, W, j-1);
            }
            done = true;
        }
        else
        {
            --i;
            --k;
            j = k;
        }
    }
    return W;
}

// Algorithm 2.4.8 (HNF modulo D)
template <class X, class INTEGER> Matrix<X > HNF_mod_D(const Matrix<X >& AA, const INTEGER& D)
{
    const X zero(0L);
    Matrix<X > A(AA);

    int m = A.rows();
    int n = A.columns();
    if (n < m)  // must have n >= m
    {
        throw std::string("HNF_mod_D : must have at least as many columns as rows");
    }

    // Step 1. [Initialize]
    int i = m;
    int j = n;
    int k = n;
    INTEGER R = D;

    std::vector<X> B;
    B.resize(m, zero);

    Matrix<X> W(m, m);

    while (1)
    {
        // Step 2. [Check zero]
        if (j != 1)
        {
            --j;
            if (A(i-1, j-1) != zero)
            {
                // Step 3. [Euclidean step]
                X u;
                X v;
                X d = extended_gcd(A(i-1,k-1), A(i-1,j-1), u, v);
                X A_ik = A(i-1,k-1) / d;
                X A_ij = A(i-1,j-1) / d;
                //std::cout << "i = " << i << ", j = " << j << ", k = " << k << ", A(i,k) = " << A(i-1,k-1) << ", A(i,j) = " << A(i-1,j-1) << ", d = " << d << ", u = " << u << ", v = " << v << std::endl;
                for (int r = 1; r <= m; ++r)
                {
                    if (r == i) continue;
                    B[r-1] = u;
                    B[r-1] *= A(r-1,k-1);
                    B[r-1] += v * A(r-1,j-1);
                    B[r-1] %= R;
                }
                A.combine_columns_mod_R(j-1, k-1, A_ik, A_ij, R);
                A(i-1,k-1) = d;
                for (int r = 1; r <= m; ++r)
                {
                    if (r == i) continue;
                    A(r-1,k-1) = B[r-1];
                }
            }
        }
        else
        {
            // Step 4. [Next row]
            X u;
            X v;
            X d = extended_gcd(A(i-1,k-1), R, u, v);
            W(i-1,i-1) = d;
            for (int r = 1; r <= m; ++r)
            {
                if (r == i) continue;
                W(r-1,i-1) = (u * A(r-1,k-1)) % R;
            }
            for (j = i + 1; j <= m; j++)
            {
                //X q = W(i-1,j-1) / d;
                X q = divide(W(i-1,j-1), d);
                W.combine_columns(j-1, i-1, q);
            }
            if (i == 1)
            {
                return W;
            }
            R /= d;
            i--;
            k--;
            j = k;
            if (A(i-1,k-1) == zero)
            {
                A(i-1,k-1) = R;
            }
        }
    }
}

#endif
