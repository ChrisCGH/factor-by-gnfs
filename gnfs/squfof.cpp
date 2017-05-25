#include <math.h>
#include "VeryLong.h"
#include "gcd.h"
#include "mod.h"
#include "timings.h"
#include <string>
#include <unordered_map>
#include <iostream>
#include <iomanip>

class QuadraticForm
{
public:
    QuadraticForm() : A_(0), B_(0), C_(0), D_(0), root_D_(0), root_D_d_(0.0), distance_(0.0)
    {}
    QuadraticForm(long int A, long int B, long int C, double distance = 0.0) : A_(A), B_(B), C_(C), distance_(distance)
    {
        calculate_discriminant();
    }

    long int A() const
    {
        return A_;
    }
    long int B() const
    {
        return B_;
    }
    long int C() const
    {
        return C_;
    }
    long long int discriminant() const
    {
        return D_;
    }
    double distance() const
    {
        return distance_;
    }

    void reduce()
    {
        long int a = C_;
        long int b = r(-B_, C_);
        // we know that b = -B_ + 2uC_, so
        //
        // b^2 - D_ = B_^2 - 4uC_B_ + 4u^2C_^2 - D_
        //          = 4A_C_ - 4uC_B_ + 4u^2C_^2
        // so
        // (b^2 - D_) / 4C_ = A_ - B_ u + C_ u^2
        // i.e.
        // c = A_ - B_ u + C_ u^2
        // c = A_ - B_ (b + B_) / 2C_ + (b + B_)^2 / 4C_
        //   = A_ + (b + B_)(b - B_) / 4C_
        //   = A_ + u * (b - B_) / 2

        long int u = (b + B_) / C_;
        u >>= 1;
        //long int c = A_ - B_ * u + C_ * u * u;
        //long int c = u * C_;
        //c -= B_;
        //c *= u;
        //c += A_;
        long int v = b - B_;
        long int c = u * v;
        c >>= 1;
        c += A_;
        //long int c = ((long long)b * b - D_) / ((long long)C_ * 4L);
        distance_ += log(fabs((B_ + root_D_d_) / (B_ - root_D_d_))) / 2.0;
        A_ = a;
        B_ = b;
        C_ = c;
#ifdef DEBUG
        long long int saveD_ = D_;
        calculate_discriminant();
        if (D_ != saveD_)
        {
            std::cerr<<"Reduce: Problem!" << std::endl;
        }
#endif
    }

    QuadraticForm rho()
    {
        long int a = C_;
        long int b = r(-B_, C_);
        // we know that b = -B_ + 2uC_, so
        //
        // b^2 - D_ = B_^2 - 4uC_B_ + 4u^2C_^2 - D_
        //          = 4A_C_ - 4uC_B_ + 4u^2C_^2
        // so
        // (b^2 - D_) / 4C_ = A_ - B_ u + C_ u^2
        // i.e.
        // c = A_ - B_ u + C_ u^2
        // c = A_ - B_ (b + B_) / 2C_ + (b + B_)^2 / 4C_
        //   = A_ + (b + B_)(b - B_) / 4C_
        //   = A_ + u * (b - B_) / 2
        long int c = (b + B_) / C_;
        c >>= 1;
        c *= (b - B_);
        c >>= 1;
        c += A_;
#ifdef DEBUG
        long int c1 = ((long long)b * b - D_) / ((long long)C_ * 4L);
        if (c1 != c)
        {
            std::cerr << "Problem : c = " << c << ", c1 = " << c1 << std::endl;
        }
#endif
        return QuadraticForm(a,b,c);
    }

    bool reduced() const
    {
        long int limit1 = 0L;
        long int limit2 = root_D_;
        if (A_ < 0) limit1 = root_D_ + 2L * A_;
        else limit1 = root_D_ - 2L * A_;
        if (limit1 < 0) limit1 = -limit1;
        if (B_ > limit1 && B_ < limit2) return true;
        return false;
    }

    friend ostream& operator<< (ostream& os, const QuadraticForm& qf)
    {
        os << "(" << qf.A_ << "," << qf.B_ << "," << qf.C_ << "), disc = " << qf.D_ << ", distance = " << qf.distance_;
        return os;
    }

    void full_reduce()
    {
        while (!reduced())
        {
            //std::cout << "3. About to call reduce(), C_ = " << C_ << std::endl;
            reduce();
#ifdef DEBUG
            std::cout << "reduces -> " << *this << std::endl;
#endif
        }
    }


private:
    long int A_;
    long int B_;
    long int C_;
    long long int D_;
    long int root_D_;
    double root_D_d_;
    double distance_;
    void calculate_discriminant()
    {
        D_ = (long long)B_ * B_ - 4L * (long long)A_ * C_;
        root_D_d_ = sqrt((double)D_);
        root_D_ = (long int)root_D_d_;
    }
    long int r(long int b, long int a)
    {
        if (a < 0) a = -a;
        long int two_a = a << 1;
        if (a > root_D_)
        {
            // trying to find rr = b mod 2a
            // with -a < rr <= a
            if (b > -a && b <= a) return b;
            long int rr = b % two_a;
            if (rr > a) rr -= two_a;
            if (rr <= -a) rr += two_a;
            return rr;
        }
        else
        {
            // trying to find rr = b mod 2a
            // with root_D_ - 2a < rr < root_D_
#if 0
            long int k = (root_D_ - b) / two_a;
            //long int rr = b + k * a + k * a;
            long int rr = b + k * two_a;
            if (rr <= root_D_ - two_a) rr += two_a;
            if (rr > root_D_) rr -= two_a;
#else
            long int rr = root_D_ - (root_D_ - b) % two_a;
#endif
            return rr;
        }
    }

};

class QuadraticFormL
{
public:
    QuadraticFormL() : A_(0), B_(0), C_(0), D_(0), root_D_(0), root_D_d_(0.0), distance_(0.0)
    {}
    QuadraticFormL(long long int A, long long int B, long long int C, double distance = 0.0) : A_(A), B_(B), C_(C), distance_(distance)
    {
        calculate_discriminant();
    }

    long long int A() const
    {
        return A_;
    }
    long long int B() const
    {
        return B_;
    }
    long long int C() const
    {
        return C_;
    }
    long long int discriminant() const
    {
        return D_;
    }
    double distance() const
    {
        return distance_;
    }

    void reduce()
    {
        long long int a = C_;
        long long int b = r(-B_, C_);
        /*
        	 // we know that b = -B_ + 2uC_, so
        	 //
        	 // b^2 - D_ = B_^2 - 4uC_B_ + 4u^2C_^2 - D_
        	 //          = 4A_C_ - 4uC_B_ + 4u^2C_^2
        	 // so
        	 // (b^2 - D_) / 4C_ = A_ - B_ u + C_ u^2
        	 // i.e.
        	 // c = A_ - B_ u + C_ u^2
        */
        long long int u = (b + B_) / C_;
        u >>= 1;
        long long int c = A_ - B_ * u + C_ * u * u;
        //long long int c = ((long long)b * b - D_) / ((long long)C_ * 4);
        distance_ += log(fabs(((double)B_ + root_D_d_) / ((double)B_ - root_D_d_))) / 2.0;
        A_ = a;
        B_ = b;
        C_ = c;
#ifdef DEBUG
        long long int saveD_ = D_;
        calculate_discriminant();
        if (D_ != saveD_)
        {
            std::cerr<<"ReduceL: Problem!" << std::endl;
        }
#endif
    }

    void inverse_reduce()
    {
        long long int c = A_;
        long long int b = r(-B_, A_);
        /*
        	 // we know that b = -B_ + 2uA_, so
        	 //
        	 // b^2 - D_ = B_^2 - 4uA_B_ + 4u^2A_^2 - D_
        	 //          = 4A_C_ - 4uA_B_ + 4u^2A_^2
        	 // so
        	 // (b^2 - D_) / 4A_ = C_ - B_ u + A_ u^2
        	 // i.e.
        	 // a = C_ - B_ u + A_ u^2
        */
        long long int u = (b + B_) / A_;
        u >>= 1;
        long long int a = C_ - B_ * u + A_ * u * u;
        distance_ += log(fabs(((double)b - root_D_d_) / ((double)b + root_D_d_))) / 2.0;
        A_ = a;
        B_ = b;
        C_ = c;
#ifdef DEBUG
        long long int saveD_ = D_;
        calculate_discriminant();
        if (D_ != saveD_)
        {
            std::cerr << "Inverse ReduceL: Problem!" << std::endl;
        }
#endif
    }

    QuadraticFormL rho()
    {
        long long int a = C_;
        long long int b = r(-B_, C_);
        long long int c = ((long long)b * b - D_) / ((long long)C_ * 4);
        return QuadraticFormL(a,b,c);
    }

    bool reduced() const
    {
        double limit1 = 0.0;
        double limit2 = root_D_;
        if (A_ < 0) limit1 = static_cast<double>(root_D_ + 2*A_);
        else limit1 = static_cast<double>(root_D_ - 2*A_);
        if (limit1 < 0) limit1 = -limit1;
        if (B_ > limit1 && B_ < limit2) return true;
        return false;
    }

    void full_reduce()
    {
        while (!reduced())
        {
            reduce();
#ifdef DEBUG
            std::cout << "reduces -> " << *this << std::endl;
#endif
        }
    }

    friend ostream& operator<< (ostream& os, const QuadraticFormL& qf)
    {
        os << "(" << qf.A_ << "," << qf.B_ << "," << qf.C_ << "), disc = " << qf.D_ << ", distance = " << qf.distance_;
        return os;
    }

    void compose(const QuadraticFormL& f2)
    {
        // composition of quadratic forms, assuming same discriminant
        if (D_ == f2.D_)
        {
            long long int beta = (B_ + f2.B_) / 2L;
            long long int x = 0L;
            long long int y = 0L;
            long long int m = extended_gcd<long long int>(A_, beta, x, y);
            long long int n = gcd<long long int>(m, f2.A_);
            long long int modulus = f2.A_ / n;
            if (modulus < 0) modulus = -modulus;
            long long int m_div_n_inv = 0L;
            if (m / n > 0) m_div_n_inv = inverse(m / n, modulus);
            else m_div_n_inv = -inverse(-m / n, modulus);
            long long int tmp_ll = (f2.B_ - B_) / 2L;
            tmp_ll *= x;
            tmp_ll -= C_ * y;
            tmp_ll %= modulus;
            tmp_ll *= m_div_n_inv;
            long long int z = tmp_ll % modulus;
            long long int a = (A_ / n) * (f2.A_ / n);
            long long int b = B_ + (A_ / n) * z * 2L;
            // c = (b^2 - D_) / 4a, but this might cause overflow
            // let b = q 2^32 + r
            // b^2 = q^2 2^64 + q r 2^33 + r^2
            //     = (q^2 2^32 + 2 q r) 2^32 + r^2 % 2^32 + (r^2 / 2^32) * 2^32
            //     = (q^2 2^32 + 2 q r + r^2 / 2^32) 2^32 + r^2 % 2^32
            // which overflows if
            // q^2 2^32 + 2 q r + r^2 / 2^32 >= 2^32
            // a = xy
            // b = B_ + 2xz
            // where x = A_ / n, y = r2.A_ / n
            // so b^2 = B_^2 + 4 xzB_ + 4x^2z^2
            // b^2 - D_ = B_^2 - D_ + 4x(zB_ + xz^2)
            //          = 4A_C_ + 4x(zB_ + xz^2)
            // (b^2 - D_) / 4a = (A_C_ + x(zB_ + xz^2)) / xy
            //                 = (nxC_ + x(zB_ + xz^2)) / xy
            //                 = (nC_ + z(B_ + xz)) / y
            //                 = (nC_ + zB_ + z^2x) / y
            // Let z^2 = q y + r
            //                 = (nC_ + zB_ + rx) / y + qx
            long long int q = (z * z) / (f2.A_ / n);
            long long int r = (z * z) % (f2.A_ / n);
            //long long int c = (n * C_ + z * (B_ + z * (A_ / n))) / (f2.A_ / n);
            long long int c = (n * C_ + z * B_ + r * (A_ / n)) / (f2.A_ / n) + q * (A_ / n);
            //long long int c = ((long long)b * b - D_) / a;
            //c >>= 2;
            A_ = a;
            B_ = b;
            C_ = c;
            distance_ += f2.distance_;
#ifdef DEBUG
            {
                long long int saveD_ = D_;
                calculate_discriminant();
                if (D_ != saveD_)
                {
                    std::cerr<<"ComposeL: Problem!" << std::endl;
                }
            }
#endif
        }
    }

private:
    long long int A_;
    long long int B_;
    long long int C_;
    long long int D_;
    long int root_D_;
    double root_D_d_;
    double distance_;

    void calculate_discriminant()
    {
        D_ = (long long)B_ * B_ - 4L * (long long)A_ * C_;
        root_D_d_ = sqrt((double)D_);
        root_D_ = (long int)root_D_d_;
    }
    long long int r(long long int b, long long int aa)
    {
        long long int rr = 0L;
        long long int a = aa;
        if (a < 0) a = -a;
        if (a > root_D_)
        {
            long int k = static_cast<long int>((a - b) / (2L * a));
            rr = b + k * a + k * a;
#if 0
            if (rr > a) rr -= a + a;
            if (rr <= -a) rr += a + a;
#endif
        }
        else
        {
            long int k = static_cast<long int>((root_D_ - b) / (2L * a));
            rr = b + k * a + k * a;
#if 0
            if (rr <= root_D_ - a - a) rr += a + a;
            if (rr >= root_D_) rr -= a + a;
#endif
        }
        return rr;
    }

};

namespace
{
// Algorithn 1.7.1 (Integer Square Root)

long int integer_square_root(long int n)
{
    if (n == 1L) return 1L;
    // 1. [Initialize]
    long int x = (long int)(sqrt((float)n) + 0.5);
    long int y = n / x;
    y += x;
    y >>= 1;
    while (y < x)
    {
        x = y;
        y = n / x;
        y += x;
        y >>= 1;
    }
    return x;
}

long int integer_square_root(long long int n)
{
    // 1. [Initialize]
    long long int x = (long long int)(sqrt((double)n) + 0.5);
    long long int y = n / x;
    //std::cout << "integer_square_root : n = " << n << ", x = " << x << ", y = " << y << std::endl;
    y += x;
    y >>= 1;
    while (y < x)
    {
        x = y;
        y = n / x;
        y += x;
        y >>= 1;
    }
    //std::cout << "integer_square_root : n = " << n << ", x = " << x << ", y = " << y << std::endl;
    return static_cast<long int>(x);
}

// Algorithm 1.7.3 (Square Test).
//
#if 0
int q64[64];
int q65[65];
int q11[11];
int q63[63];

void square_test_precomputations()
{
    for (int k = 0; k < 11; ++k) q11[k] = 0;
    for (int k = 0; k < 6; ++k) q11[(k*k)%11] = 1;

    for (int k = 0; k < 63; ++k) q63[k] = 0;
    for (int k = 0; k < 32; ++k) q63[(k*k)%63] = 1;

    for (int k = 0; k < 64; ++k) q64[k] = 0;
    for (int k = 0; k < 32; ++k) q64[(k*k)%64] = 1;

    for (int k = 0; k < 65; ++k) q65[k] = 0;
    for (int k = 0; k < 33; ++k) q65[(k*k)%65] = 1;
}

long int square_test(long int n)
{
#if 0
    static int first_time = 1;
    if (first_time)
    {
        first_time = 0;
        square_test_precomputations();
    }
    if (n < 0) return 0;

    if (q64[(n & 0x3F)] == 0) return 0; // not a square
    long int r = n % 45045;
    if (q63[(r % 63)] == 0) return 0; // not a square
    if (q65[(r % 65)] == 0) return 0; // not a square
    if (q11[(r % 11)] == 0) return 0; // not a square
    // [Compute square root]
    long int q = integer_square_root(n);
    if (n != q*q) return 0;
    return q;
#else
    long int q = static_cast<long int>(sqrt(static_cast<double>(n)));
    if (n != q*q) return 0;
    return q;
#endif
}
#endif
long int square_test(long long n)
{
    static int q64[64];
    static int q65[65];
    static int q11[11];
    static int q63[63];
    static int first_time = 1;
    if (first_time)
    {
        first_time = 0;
        memset(q11, 0, 11 * sizeof(int));
        //for (int k = 0; k < 11; ++k) q11[k] = 0;
        for (int k = 0; k < 6; ++k) q11[(k*k)%11] = 1;

        memset(q63, 0, 63 * sizeof(int));
        //for (int k = 0; k < 63; ++k) q63[k] = 0;
        for (int k = 0; k < 32; ++k) q63[(k*k)%63] = 1;

        memset(q64, 0, 64 * sizeof(int));
        //for (int k = 0; k < 64; ++k) q64[k] = 0;
        for (int k = 0; k < 32; ++k) q64[(k*k)%64] = 1;

        memset(q65, 0, 65 * sizeof(int));
        //for (int k = 0; k < 65; ++k) q65[k] = 0;
        for (int k = 0; k < 33; ++k) q65[(k*k)%65] = 1;
        // square_test_precomputations();
    }

    if (n < 0)
    {
        //std::cout << "1. not a square" << std::endl;
        return 0;
    }
    if (q64[(n & 0x3f)] == 0)
    {
        //std::cout << "2. not a square" << std::endl;
        return 0L; // not a square
    }
    long int r = static_cast<long int>(n % 45045L);
    if (q63[(r % 63)] == 0)
    {
        //std::cout << "3. not a square" << std::endl;
        return 0L; // not a square
    }
    if (q65[(r % 65)] == 0)
    {
        //std::cout << "4. not a square" << std::endl;
        return 0L; // not a square
    }
    if (q11[(r % 11)] == 0)
    {
        //std::cout << "5. not a square" << std::endl;
        return 0L; // not a square
    }
    // [Compute square root]
    long int q = integer_square_root(n);
    if (n != (long long)q*q)
    {
        //std::cout << "6. not a square" << std::endl;
        //std::cout << "n = " << n << std::endl;
        //std::cout << "q = " << q << std::endl;
        //std::cout << "q*q = " << (long long)q*q << std::endl;
        return 0L;
    }
    //std::cout << "7. a square" << std::endl;
    return q;
}
const long long int max_D_0 = 1152921504606846976LL;
//const long long int max_D_0 = 2305843009213693952LL;
//const long long int max_D = 4611686018427387904LL;
//const long long int max_D = 6000000000000000000LL;
//const long long int max_D =  9223372036854775808LL;

};

struct squfof_info
{
    long long int D_;
    long int multiplier_;
    long int d_;
    long int b_;
    long int L_;
    long int a_;
    std::unordered_map<long int, bool> Q_;
    QuadraticForm f_;
    QuadraticFormL g_;
    long int reductions_;
    long int mask_;
    std::vector<QuadraticForm> fp0_;

    squfof_info() : D_(0), multiplier_(0), d_(0), b_(0), L_(0), a_(0), reductions_(0), mask_(0x1)
    {}

    squfof_info(long long int N, long int multiplier) : multiplier_(multiplier), reductions_(0), mask_(0x1)
    {
        D_ = N * multiplier_;
        //if (D_ > max_D || D_ < 0 || D_ / multiplier_ != N)
        if (D_ > max_D_0 || D_ < 0 || D_ / multiplier_ != N)
        {
            D_ = 0LL;
            multiplier_ = 0L;
            d_ = 0L;
            b_ = 0L;
            return;
        }
        d_ = integer_square_root(D_);
        if (multiplier_ == 1L)
        {
            b_ = d_ - 1L;
            b_ >>= 1;
            b_ <<= 1;
            ++b_;
        }
        else
        {
            //switch (multiplier_)
            switch (multiplier_ % 4L)
            {
            case 3L:
            {
                b_ = d_ - 1L;
                b_ >>= 1;
                b_ <<= 1;
                ++b_;
            }
            break;
            case 0L:
            {
                b_ = d_;
                b_ >>= 1;
                b_ <<= 1;
            }
            break;
            case 1L:
            {
                b_ = d_;
                b_ >>= 1;
                b_ <<= 1;
                --b_;
            }
            break;
            }
        }

        long long int tmp(b_);
        tmp *= tmp;
        tmp -= D_;
        tmp /= 4L;
        if (D_ < max_D_0)
        {
            f_ = QuadraticForm(1L, b_, static_cast<long int>(tmp));
        }
        Q_.clear();
        L_ = integer_square_root(d_);
        a_ = 0L;
    }

    void reduce_f()
    {
        reductions_++;
        //std::cout << "5. About to call reduce(), f_.C_ = " << f_.C() << std::endl;
        f_.reduce();
        if (reductions_ == mask_)
        {
            mask_ <<= 1;
            fp0_.push_back(f_);
        }
    }

    double f_distance() const
    {
        return f_.distance();
    }

    bool f_is_identity()
    {
        return (f_.A() == 1L);
    }

    void fill_queue()
    {
        long int a = f_.A();
        if (a < 0) a = -a;
        if (a <= L_) Q_[a] = true;
    }

    long int check_square_factor()
    {
        long int s = gcd<long int>(a_, f_.B());
        return static_cast<long int>(gcd<long long int>(s, D_));
    }

    void make_g()
    {
        g_ = QuadraticFormL(a_, -f_.B(), (long long)a_*(f_.C()));
    }

    long int square_test_f()
    {
        return square_test(f_.A());
    }

    friend ostream& operator<< (ostream& os, squfof_info& si)
    {
        os << si.f_;
        return os;
    }
};

#ifdef DEBUG
#define QFOUT(s1, qf, s2) std::cout << s1 << qf << s2 << std::endl;
#else
#define QFOUT(s1, qf, s2)
#endif

// Shanks' SQUFOF algorithm [Algorithm 8.7.2]
bool SQUFOF(long long int N, long int& factor)
{
    if (N > max_D_0) return false;
    // 1. [Is N prime?]
    VeryLong N_vl(N);
    if (N_vl.is_probable_prime()) return false;

    // 2. [Is N square?]
    long int n = square_test(N);
    if (n != 0L)
    {
        factor = n;
        return true;
    }

    // 3. [Initializations]
    const int info_list_size = 5;
    squfof_info info_list[info_list_size];
    int info_count = 3;

    if (N % 4L != 1L)
    {
        // N == 3 mod 4
        // N*3 = 1 mod 4
        // N*4 = 0 mod 4
        // N*7 = 1 mod 4
        // N*15 = 1 mod 4
        info_list[0] = squfof_info(N, 3L);
        info_list[1] = squfof_info(N, 4L);
        info_list[2] = squfof_info(N, 7L);
#if 0
        info_list[3] = squfof_info(N, 15L);
        info_list[4] = squfof_info(N, 19L);
#endif
    }
    else
    {
        info_list[0] = squfof_info(N, 1L);
        info_list[1] = squfof_info(N, 5L);
        info_list[2] = squfof_info(N, 9L);
        //info_count = 1;
#if 0
        info_list[3] = squfof_info(N, 13L);
        info_list[4] = squfof_info(N, 17L);
#endif
    }

    for (int ic = info_count - 1; ic >= 0; --ic)
    {
        if (info_list[ic].multiplier_ == 0) --info_count;
    }

    if (info_count <= 0) return false;

    int i = 0;
    const int max_i = 20000;
#ifdef DEBUG
    for (int ic = 0; ic < info_count; ++ic)
    {
        std::cout << "ic = " << ic << ", f = " << info_list[ic] << std::endl;
    }
#endif

    squfof_info* infop = 0;
    squfof_info* found_infop = 0;
    bool done_outer = false;
    while (!done_outer)
    {
        bool done = false;
        while (!done)
        {
            // 4. [Apply rho]
            for (infop = info_list; infop != info_list + info_count; ++infop)
            {
                infop->reduce_f();
            }
            ++i;
            if (i >= max_i) return false;
#ifdef DEBUG
            for (int ic = 0; ic < info_count; ++ic)
            {
                std::cout << "ic = " << ic << ", f = " << info_list[ic] << std::endl;
            }
#endif
            if (i % 2 == 0)
            {
                // 5. [Squareform?]
                squfof_info* clear_info = 0;
                for (infop = info_list; !done && !clear_info && infop != info_list + info_count; ++infop)
                {
                    infop->a_ = infop->square_test_f();
                    if (infop->a_ > 1 && infop->Q_.find(infop->a_) == infop->Q_.end())
                    {
                        found_infop = infop;
                        done = true;
                    }
                    else
                    {
                        // 6. [Short period?]
                        if (infop->f_is_identity())
                        {
                            clear_info = infop;
                        }
                    }
                }
                if (clear_info)
                {
                    if (info_count == 1)
                    {
                        std::cerr << "SQUFOF: Ran through the i elements of the principal cycle without finding a non-trivial squareform : N = " << N << std::endl;
                        return false;
                    }
                    squfof_info* tmp_infop1 = info_list;
                    squfof_info* tmp_infop2 = info_list;
                    int new_info_count = info_count;
                    while (tmp_infop1 != info_list + info_count)
                    {
                        if (!tmp_infop1->f_is_identity())
                        {
                            *tmp_infop2 = *(tmp_infop1);
                            ++tmp_infop2;
                        }
                        else
                        {
                            --new_info_count;
                        }
                        ++tmp_infop1;
                    }
                    info_count = new_info_count;
                    if (info_count == 0) return false;
                }
            }
            // 7. [Fill queue and cycle]
            for (infop = info_list; !done && infop != info_list + info_count; ++infop)
            {
                infop->fill_queue();
            }
        }

        infop = found_infop;

        // 8. [Initialize back-cycle]
        long int s = infop->check_square_factor();
        bool skip = false;

        if (s > 1)
        {
            factor = s*s;
            if (factor != infop->multiplier_) return true;
            else skip = true;
        }
        if (!skip)
        {
            int g_count = 0;
            infop->make_g();
#ifdef DEBUG
            std::cout << "g = " << infop->g_ << std::endl;
#endif
            while (!infop->g_.reduced())
            {
                //std::cout << "6. About to call reduce(), infop->g_.C_ = " << infop->g_.C() << std::endl;
                infop->g_.reduce();
                ++g_count;
#ifdef DEBUG
                std::cout << g_count << " : g = " << infop->g_ << " (reducing)" << std::endl;
#endif
            }

#ifdef DEBUG
            std::cout << g_count << " : g = " << infop->g_ << " (reduced)" << std::endl;
#endif
            // 9. [Back-cycle]
#define FASTRETURN 1
#ifdef FASTRETURN
            // Implement "Fast Return" using Shanks' distance to help find a form f
            // that is distance 1/2 infop->distance(), and compose it with infop->g_;
            double delta = infop->f_distance() / 2.0;
            double save_delta = delta;
#ifdef DEBUG
            std::cout << "Aiming for delta = " << delta << std::endl;
#endif
#ifdef DEBUG
            for (int i = 0; i < infop->fp0_.size(); ++i)
            {
                std::cout << "rho^" << i << "(I) = " << infop->fp0_[i] << std::endl;
            }
#endif
            auto it = infop->fp0_.rbegin();
//	 while (it->distance() >= delta && it != infop->fp0_.rend()) ++it;
            while (it != infop->fp0_.rend() && it->distance() >= delta) ++it;
            if (it == infop->fp0_.rend()) return false;

            QuadraticFormL F(it->A(), it->B(), it->C(), it->distance());
            delta -= it->distance();
            ++it;
            while (it != infop->fp0_.rend())
            {
                if (it->distance() < delta)
                {
                    delta -= it->distance();
#ifdef DEBUG
                    std::cout << "composing with " << *it << std::endl;
#endif
                    F.compose(QuadraticFormL(it->A(), it->B(), it->C(), it->distance()));
#ifdef DEBUG
                    std::cout << "-> " << F << ", delta = " << delta << std::endl;
#endif
                    F.full_reduce();
                    delta = save_delta - F.distance();
#ifdef DEBUG
                    std::cout << "reduced -> " << F << ", delta = " << delta << std::endl;
#endif
                }
                ++it;
            }
#ifdef DEBUG
            std::cout << "composing with " << infop->g_ << std::endl;
#endif
            infop->g_.compose(F);
#ifdef DEBUG
            std::cout << "-> " << infop->g_ << std::endl;
#endif

            while (!infop->g_.reduced())
            {
                //std::cout << "7. About to call reduce(), infop->g_.C_ = " << infop->g_.C() << std::endl;
                infop->g_.reduce();
#ifdef DEBUG
                std::cout << "reduces -> " << infop->g_ << std::endl;
#endif
            }
#ifdef DEBUG
            std::cout << "reduced -> " << infop->g_ << std::endl;
#endif
            infop->g_.inverse_reduce();
#ifdef DEBUG
            std::cout << "inverse reduced -> " << std::setprecision(20) << infop->g_ << std::endl;
#endif
            while (infop->g_.distance() > save_delta)
            {
                // we've overshot, apply inverse of rho
                infop->g_.inverse_reduce();
#ifdef DEBUG
                std::cout << "inverse reduced -> " << std::setprecision(20) << infop->g_ << std::endl;
#endif
            }

            long int b1 = 0;
            while (fabs(infop->g_.distance() - save_delta) > 1e-3 && b1 != infop->g_.B())
            {
                b1 = static_cast<long int>(infop->g_.B());
                //std::cout << "8. About to call reduce(), infop->g_.C_ = " << infop->g_.C() << std::endl;
                infop->g_.reduce();
#ifdef DEBUG
                std::cout << "reduces -> " << std::setprecision(20) << infop->g_ << std::endl;
#endif
            }
#else
            done = false;
            while (!done)
            {
                long int b1 = infop->g_.B();
                //std::cout << "9. About to call reduce(), infop->g_.C_ = " << infop->g_.C() << std::endl;
                infop->g_.reduce();
                ++g_count;
#ifdef DEBUG
                std::cout << g_count << " : g = " << infop->g_ << std::endl;
#endif
                if (b1 == infop->g_.B()) done = true;
            }
#endif

            infop->a_ = static_cast<long int>(infop->g_.A());
            if (infop->a_ < 0) infop->a_ = -infop->a_;
            if (infop->a_ % 2L) factor = infop->a_;
            else factor = infop->a_/2;
            if (factor % infop->multiplier_ == 0) factor /= infop->multiplier_;
            else
            {
                int g = gcd<long int>(factor, infop->multiplier_);
                if (g > 1L) factor /= g;
            }

            if (factor != 1) return true;
        }
        if (info_count == 1) return false;
#ifdef DEBUG
        std::cout << "Retry!" << std::endl;
#endif
        squfof_info* tmp_infop = infop;
        while (tmp_infop != info_list + info_count - 1)
        {
            *tmp_infop = *(tmp_infop + 1);
            ++tmp_infop;
        }
        --info_count;
    } // done_outer
    return true;
}

