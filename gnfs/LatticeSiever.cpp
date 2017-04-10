#include "LatticeSiever.h"
#include "Matrix.h"
#include "LongModular.h"
#include "mod.h"
#include "lll.h"
#include <math.h>
#include <set>
#include <cstdlib>

//#define RANGE_CHECK 1
#ifdef RANGE_CHECK
class SieveError
{
public:
    SieveError(long int p, long int e, long int f, long int c, long int d,
               std::pair<long int, long int> e1,
               std::pair<long int, long int> e2)
            : p_(p), e_(e), f_(f), c_(c), d_(d), e1_(e1), e2_(e2)
    {}
    friend std::ostream& operator<<(std::ostream& os, SieveError& se)
    {
        os << "Attempted to hit sieve at (c,d) = (" << se.c_ << "," << se.d_ << ") <- (e,f) = (" << se.e_ << "," << se.f_ << ") for p = " << se.p_ << std::endl;
        os << "e1 = (" << se.e1_.first << "," << se.e1_.second << ")" << std::endl;
        os << "e2 = (" << se.e2_.first << "," << se.e2_.second << ")" << std::endl;
        return os;
    }
private:
    long int p_;
    long int e_;
    long int f_;
    long int c_;
    long int d_;
    std::pair<long int, long int> e1_;
    std::pair<long int, long int> e2_;
};
#endif

// assume n > 0
namespace
{
//template <> long int inverse(long int a, long int modulus)
#if 0
long int inverse(long int a, long int modulus)
{
    long int q = 1L;
    long int rem = a;
    long int dividend = modulus;
    long int divisor = a;
    long int ps1 = 1L;
    long int ps2 = 0;
    long int parity = 0;

    while (divisor > 1)
    {
        rem = dividend - divisor;
        long int t = rem - divisor;
        if (t >= 0) //1
        {
            q += ps1;
            rem = t;
            t -= divisor;
            if (t >= 0) //2
            {
                q += ps1;
                rem = t;
                t -= divisor;
                if (t >= 0) //3
                {
                    q += ps1;
                    rem = t;
                    t -= divisor;
                    if (t >= 0) //4
                    {
                        q += ps1;
                        rem = t;
                        t -= divisor;
                        if (t >= 0) //5
                        {
                            q += ps1;
                            rem = t;
                            t -= divisor;
                            if (t >= 0) //6
                            {
                                q += ps1;
                                rem = t;
                                t -= divisor;
                                if (t >= 0) //7
                                {
                                    q += ps1;
                                    rem = t;
                                    t -= divisor;
                                    if (t >= 0) //8
                                    {
                                        q += ps1;
                                        rem = t;
                                        if (rem >= divisor)
                                        {
                                            q = dividend / divisor;
                                            rem = dividend - q * divisor;
                                            q *= ps1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        q += ps2;
        parity = ~parity;
        dividend = divisor;
        divisor = rem;
        ps2 = ps1;
        ps1 = q;
    }

    if (parity == 0)
        return (ps1);
    else
        return (modulus - ps1);
}
#endif
    bool verbose()
    {
        if (std::getenv("LATTICE_SIEVER_VERBOSE"))
        {
            return true;
        }
        return false;
    }
}
//LatticeSiever::IPrimeFactorList* LatticeSiever::SieveCacheItem::pf_list_ = 0;
LatticeSiever::PrimeFactorList* LatticeSiever::SieveCacheItem::pf_list_ = 0;
long int LatticeSiever::total_relations_ = 0;
double LatticeSiever::total_sieving_time_ = 0;

namespace
{
// calculate log of x to base q
double logq(double x, int q)
{
    if (q == 10) return log10(x);
    return log10(x) / log10((double)q);
}

VeryLong evaluate_on_lattice(const Polynomial<VeryLong>& f,
                             long int c, long int d,
                             const std::pair<long int, long int>& c1,
                             const std::pair<long int, long int>& c2)
{
    const VeryLong zero(0L);
    if (c == 0L && d == 0L) return zero;
    VeryLong a(c);
    a *= c1.first;
    VeryLong tmp(d);
    tmp *= c2.first;
    a += tmp;
    VeryLong b(c);
    b *= c1.second;
    tmp = d;
    tmp *= c2.second;
    b += tmp;
    VeryLong res = f.evaluate_homogeneous(a, b);
    //std::cerr << "F(" << a << "," << b << ") = " << res << std::endl;
    return res;
}

double evaluate_on_lattice(const Polynomial<double>& f,
                           long int c, long int d,
                           const std::pair<long int, long int>& c1,
                           const std::pair<long int, long int>& c2)
{
    if (c == 0L && d == 0L) return 0.0;
    double a(c);
    a *= c1.first;
    double tmp(d);
    tmp *= c2.first;
    a += tmp;
    double b(c);
    b *= c1.second;
    tmp = d;
    tmp *= c2.second;
    b += tmp;
    double res = f.evaluate_homogeneous(a, b);
    return res;
}

inline long int nearest_integer(double d)
{
    if (d < 0) return( -(long int)(-d + 0.5));
    return((long int)(d + 0.5));
    /*   if (d < 0)
            return std::ceil(d - 0.5);
       else
            return std::floor(d + 0.5);*/
}

long int max_total_relations = 25000L;
};

LatticeSiever::LatticeSiever(const std::string& config_file)
        : relfile_(0), alg_factor_base_(0), rat_factor_base_(0),
        potentially_smooth_point_(0), number_potentially_smooth_(0),
        c1_(std::make_pair<long int, long int>(0L, 0L)),
        c2_(std::make_pair<long int, long int>(0L, 0L)),
        c_region_(Point<double>(min_c, min_d), Point<double>(max_c, min_d), Point<double>(min_c, max_d)),
        rat_pf_list_(rat_pf_list_size, fixed_sieve_array_),
        alg_pf_list_(alg_pf_list_size, fixed_sieve_array_),
        sieveCache_(fixed_sieve_array_, sieve_bit_array_),
        timer_("sieve.tim")
{
    SieveConfig config(config_file.c_str());
    if (debug_)
    {
        config.display();
    }
    f1_ = config.f1();
    f1d_ = Polynomial<VeryLong>::convert_to_double<double>(f1_);
    f2_ = config.f2();
    f2d_ = Polynomial<VeryLong>::convert_to_double<double>(f2_);
    B1_ = config.B1();
    L1_ = config.L1();
    LP1_ = config.LP1();
    B2_ = config.B2();
    L2_ = config.L2();
    LP2_ = config.LP2();
    L_LP_1_ = L1_;
    for (int i = 0; i < LP1_; i++) L_LP_1_ *= L1_;
    L_LP_2_ = L2_;
    for (int i = 0; i < LP2_; i++) L_LP_2_ *= L2_;
    relation_file_ = config.RELATION_FILE();
    SIEVE_BOUND_ADJUSTMENT1_ = config.SIEVE_BOUND_ADJUSTMENT1();
    SIEVE_BOUND_ADJUSTMENT2_ = config.SIEVE_BOUND_ADJUSTMENT2();
    SMALL_PRIME_BOUND1_ = config.SMALL_PRIME_BOUND1();
    SMALL_PRIME_BOUND2_ = config.SMALL_PRIME_BOUND2();
    INITIAL_CUTOFF_ = config.INITIAL_CUTOFF();
    MIN_A_ = config.MIN_A();
    MAX_A_ = config.MAX_A();
    MIN_B_ = config.MIN_B();
    MAX_B_ = config.MAX_B();
    debug_ = config.DEBUG();
    fixed_sieve_region_ = config.FIXED_SIEVE_REGION();

    SKEWEDNESS_ = config.SKEWEDNESS();

    std::string tmp = config.SIEVE_ID() + ".fb.dat";
    const char* alg_base_file = tmp.c_str();
    std::string tmp1 = config.SIEVE_ID() + ".rb.dat";
    const char* rat_base_file = tmp1.c_str();
    try
    {
        alg_factor_base_ = new FactorBase(alg_base_file);
    }
    catch (...)
    {
        delete alg_factor_base_;
        alg_factor_base_ = new FactorBase(f1_, B1_, alg_base_file);
    }

    try
    {
        rat_factor_base_ = new FactorBase(rat_base_file);
    }
    catch (...)
    {
        delete rat_factor_base_;
        rat_factor_base_ = new FactorBase(f2_, B2_, rat_base_file);
    }
}

LatticeSiever::~LatticeSiever()
{
    delete alg_factor_base_;
    delete rat_factor_base_;
    delete [] potentially_smooth_point_;
}
// The offset of (c,d) from start of sieve_array_ is
//
//    c - min_c + (max_c - min_c + 1) * d
//
std::pair<long int, long int> LatticeSiever::offset_to_c_d(size_t offset)
{
    std::pair<long int, long int> cd;

    cd.first = offset % (max_c - min_c + 1) + min_c;
    cd.second = offset / (max_c - min_c + 1);
    return cd;
}

size_t LatticeSiever::c_d_to_offset(const std::pair<long int, long int>& cd)
{
    return cd.first - min_c + (max_c - min_c + 1) * cd.second;
}

long int LatticeSiever::check_interval1(long int q)
{
    double L1d = L1_;
    double L1d2 = L1d;
    for (int i = 0; i < LP1_ - 1; i++) L1d2 *= L1d;
    L1d2 /= static_cast<double>(q);
    double log_L1d2 = logq(L1d2, LOGQ_BASE);

    int adjustment = SIEVE_BOUND_ADJUSTMENT1_;
    int cutoff0 = INITIAL_CUTOFF_;
    long int potential = 0;

    SIEVE_TYPE* sieve_ptr = fixed_sieve_array_;
    SIEVE_TYPE* sieve_end_ptr = fixed_sieve_array_ + fixed_sieve_array_size;
#if 0
    long int c = min_c;
    long int d = min_d;
#endif

    while (sieve_ptr < sieve_end_ptr)
    {
        SIEVE_TYPE val = -128;
        if (debug_)
        {
            std::pair<long int, long int> cd = offset_to_c_d(sieve_ptr - fixed_sieve_array_);
            std::cerr << "1. (c,d) = (" << cd.first << "," << cd.second << "), *sieve_ptr = " << (int)(*sieve_ptr) << ", cutoff0 = " << cutoff0 << std::endl;
        }
        if (!sieve_bit_array_.isSet(sieve_ptr - fixed_sieve_array_)) 
        {
            if ((int)(*sieve_ptr) >= cutoff0)
            {
#if 1
                std::pair<long int, long int> cd = offset_to_c_d(sieve_ptr - fixed_sieve_array_);
                double value1 = evaluate_on_lattice(f1d_, cd.first, cd.second, c1_, c2_);
#else
                double value1 = evaluate_on_lattice(f1d_, c, d, c1_, c2_);
#endif
                int cutoff = static_cast<int>(logq(fabs(value1), LOGQ_BASE) - log_L1d2);
                cutoff -= adjustment;
                if (debug_)
                {
                    std::pair<long int, long int> cd = offset_to_c_d(sieve_ptr - fixed_sieve_array_);
                    std::cerr << "2. (c,d) = (" << cd.first << "," << cd.second << "), *sieve_ptr = " << (int)(*sieve_ptr) << ", cutoff = " << cutoff << std::endl;
                }
    
                if ((int)(*sieve_ptr) > cutoff)
                {
                    if (debug_)
                    {
                        size_t offset = sieve_ptr - fixed_sieve_array_;
                        std::pair<long int, long int> cd = offset_to_c_d(offset);
                        std::cerr << "2.1 (c,d) = (" << cd.first << "," << cd.second << "), sieve_ptr = " << std::hex << (size_t)sieve_ptr << std::dec << ", offset = " << offset << ", cd = (" << cd.first << "," << cd.second << "), *sieve_ptr = " << (int)(*sieve_ptr) << ", cutoff = " << cutoff << std::endl;
                    }
                    val = 0;
                    ++potential;
                }
            }
    
            if (val < 0)
            {
                sieve_bit_array_.set(sieve_ptr - fixed_sieve_array_);
            }
            else
            {
                *sieve_ptr = val;
            }
        }
        ++sieve_ptr;
#if 0
        ++c;
        if (c > max_c)
        {
            c = min_c;
            ++d;
        }
#endif
    }
    return potential;
}

void LatticeSiever::check_interval2()
{
    double L1d = L2_;
    double L1d2 = L1d;
    for (int i = 0; i < LP2_ - 1; i++) L1d2 *= L1d;
    double log_L1d2 = logq(L1d2, LOGQ_BASE);

    int adjustment = SIEVE_BOUND_ADJUSTMENT2_;
    SIEVE_TYPE* sieve_ptr = fixed_sieve_array_;
    SIEVE_TYPE* sieve_end_ptr = fixed_sieve_array_ + fixed_sieve_array_size;
    long int c = min_c;
    long int d = min_d;

    int cutoff0 = 0;
    while (sieve_ptr < sieve_end_ptr)
    {
        if (!sieve_bit_array_.isSet(sieve_ptr - fixed_sieve_array_)) 
        {
            if ((int)(*sieve_ptr) >= cutoff0)
            {
                double value1 = evaluate_on_lattice(f2d_, c, d, c1_, c2_);
                int cutoff = static_cast<int>(logq(fabs(value1), LOGQ_BASE) - log_L1d2);
                cutoff -= adjustment;
                if (debug_)
                {
                    std::cerr << "3. (c,d) = (" << c << "," << d << "), *sieve_ptr = " << (int)(*sieve_ptr) << ", cutoff = " << cutoff << std::endl;
                }
                if ((int)(*sieve_ptr) > cutoff)
                {
                    if (debug_)
                    {
                        std::cerr << "3.1 (c,d) = (" << c << "," << d << "), *sieve_ptr = " << (int)(*sieve_ptr) << ", cutoff = " << cutoff << std::endl;
                    }
                    VeryLong v = abs(evaluate_on_lattice(f2_, c, d, c1_, c2_));
                    potentially_smooth_point_[number_potentially_smooth_] = PotentiallySmoothPoint(c, d, sieve_ptr, v);
                    if (number_potentially_smooth_ > 0)
                    {
                        potentially_smooth_point_[number_potentially_smooth_ - 1].next_ =
                            potentially_smooth_point_ + number_potentially_smooth_;
                    }
                    ++number_potentially_smooth_;
                }
                else
                {
                    //*sieve_ptr = -128;
                    sieve_bit_array_.set(sieve_ptr - fixed_sieve_array_);
                }
            }
            else
            {
                //*sieve_ptr = -128;
                sieve_bit_array_.set(sieve_ptr - fixed_sieve_array_);
            }
        }
        ++c;
        ++sieve_ptr;
        if (sieve_ptr < sieve_end_ptr && c > max_c)
        {
            c = min_c;
            ++d;
        }
    }
}

void LatticeSiever::divide_by_small_primes1(PotentiallySmoothPoint* smooth_iter)
{
    if (smooth_iter->c_ == 4994L && smooth_iter->d_ == 3578L)
    {
        std::cerr << "divide_by_small_primes1 : remaining_quotient = " << smooth_iter->partial1_.remaining_quotient_ << std::endl;
    }
    VeryLong::generate_prime_table();
    long int p = VeryLong::firstPrime();
    long int B = std::max(SMALL_PRIME_BOUND1_, 1000L);
    while (p < B)
    {
        if (smooth_iter->c_ == 4994L && smooth_iter->d_ == 3578L)
        {
            std::cerr << "divide_by_small_primes1 : p = " << p << ", remaining_quotient = " << smooth_iter->partial1_.remaining_quotient_ << std::endl;
        }
        while (smooth_iter->partial1_.remaining_quotient_ % p == 0L)
        {
            smooth_iter->add_factor1(p);
        }
        p = VeryLong::nextPrime();
    }
    if (smooth_iter->c_ == 4994L && smooth_iter->d_ == 3578L)
    {
        std::cerr << "divide_by_small_primes1 : final remaining_quotient = " << smooth_iter->partial1_.remaining_quotient_ << std::endl;
    }
}

void LatticeSiever::divide_by_small_primes2(PotentiallySmoothPoint* smooth_iter)
{
    //std::cerr << "divide_by_small_primes2 : (c,d) = (" << smooth_iter->c_ << "," << smooth_iter->d_ << ")" << std::endl;
    if (smooth_iter->partial2_.remaining_quotient_ == 0L)
    {
        throw std::string("Problem: remaining_quotient is zero");
    }
    if (debug_)
    {
        std::cerr << "divide_by_small_primes2 : remaining_quotient = " << smooth_iter->partial2_.remaining_quotient_ << std::endl;
    }
    VeryLong::generate_prime_table();
    long int p = VeryLong::firstPrime();
    long int B = std::max(SMALL_PRIME_BOUND2_, 1000L);
    while (p < B)
    {
        while (smooth_iter->partial2_.remaining_quotient_ % p == 0L)
        {
            if (debug_)
            {
                std::cerr << "divide_by_small_primes2 : p = " << p << ", remaining_quotient = " << smooth_iter->partial2_.remaining_quotient_ << std::endl;
            }
            smooth_iter->add_factor2(p);
        }
        p = VeryLong::nextPrime();
    }
}

bool LatticeSiever::is_likely_to_be_smooth(const FastVeryLong& remaining_quotient, long int B, long int L, long int LP, const VeryLong& L_LP)
{
    if (remaining_quotient <= B) return true;
    bool is_prime = remaining_quotient.is_probable_prime();
    if (LP > 0 && remaining_quotient <= L && is_prime) return true;
    // Here if
    //    remaining_quotient is > B and
    //    LP == 0 or remaining_quotient > L or remaining_quotient is composite
    // But remember that we have already removed factors <= B, so all factors
    // of remaining_quotient must be > B, which means ...
    if (LP == 0) return false;
    // Here if
    //    remaining_quotient is > B and LP > 0 and
    //    remaining_quotient > L or remaining_quotient is composite
    // But if remaining_quotient is prime, then it's a prime > L, so ...
    if (is_prime) return false;
    // Here if
    //    remaining_quotient is > B and LP > 0 and remaining_quotient is composite
    // note: remaining_quotient could still be <= L
    if (remaining_quotient <= L) return true;
    // Check if it's possible to be smooth still
    //if (remaining_quotient / L_LP > B) return false;
    if (remaining_quotient / B > L_LP) return false;
    // remaining_quotient could still be non-smooth, but it's quite likely to be smooth
    // Don't do the final factorisation here, but at the end.
    return true;
}

void LatticeSiever::eliminate1(long int q)
{
    // This function tries to do the final factorisation
    // and if the partial factorisation isn't smooth it eliminates it
    PotentiallySmoothPoint* smooth_iter = head_psp_;
    PotentiallySmoothPoint* prev_smooth_iter = 0;
    int nonzeros = 0;
    while (smooth_iter)
    {
        divide_by_small_primes1(smooth_iter);
        smooth_iter->add_factor1(q);
        if (is_likely_to_be_smooth(smooth_iter->partial1_.remaining_quotient_, B1_, L1_, LP1_, L_LP_1_))
        {
            if (debug_)
            {
                std::cerr << "(c,d) = (" << smooth_iter->c_ << "," << smooth_iter->d_ << "), partial1_.remaining_quotient = " << smooth_iter->partial1_.remaining_quotient_ << ", is likely to be smooth" << std::endl;
            }
            ++nonzeros;
            prev_smooth_iter = smooth_iter;
        }
        else
        {
            if (debug_)
            {
                std::cerr << "(c,d) = (" << smooth_iter->c_ << "," << smooth_iter->d_ << "), partial1_.remaining_quotient = " << smooth_iter->partial1_.remaining_quotient_ << ", is not likely to be smooth" << std::endl;
            }
            if (prev_smooth_iter)
            {
                prev_smooth_iter->next_ = smooth_iter->next_;
            }
            else
            {
                head_psp_ = smooth_iter->next_;
            }
        }
        smooth_iter = smooth_iter->next_;
    }
    number_potentially_smooth_ = nonzeros;
}

void LatticeSiever::eliminate2()
{
    // This function tries to do the final factorisation
    // and if the partial factorisation isn't smooth it eliminates it
    PotentiallySmoothPoint* smooth_iter = head_psp_;
    PotentiallySmoothPoint* prev_smooth_iter = 0;
    int zeros = 0;
    while (smooth_iter)
    {
        divide_by_small_primes2(smooth_iter);
        if (is_likely_to_be_smooth(smooth_iter->partial2_.remaining_quotient_, B2_, L2_, LP2_, L_LP_2_))
        {
            if (debug_)
            {
                std::cerr << "(c,d) = (" << smooth_iter->c_ << "," << smooth_iter->d_ << "), partial2_.remaining_quotient = " << smooth_iter->partial2_.remaining_quotient_ << ", is likely to be smooth" << std::endl;
            }
            VeryLong v = abs(evaluate_on_lattice(f1_, smooth_iter->c_, smooth_iter->d_, c1_, c2_));
            smooth_iter->partial1_.remaining_quotient_ = v;
            smooth_iter->partial1_.factor_.clear();
            --zeros;
            prev_smooth_iter = smooth_iter;
        }
        else
        {
            if (debug_)
            {
                std::cerr << "(c,d) = (" << smooth_iter->c_ << "," << smooth_iter->d_ << "), partial2_.remaining_quotient = " << smooth_iter->partial2_.remaining_quotient_ << ", is not likely to be smooth" << std::endl;
            }
            sieve_bit_array_.set(smooth_iter->ptr_ - fixed_sieve_array_);
            if (prev_smooth_iter)
            {
                prev_smooth_iter->next_ = smooth_iter->next_;
            }
            else
            {
                head_psp_ = smooth_iter->next_;
            }
        }
        smooth_iter = smooth_iter->next_;
        ++zeros;
    }
    number_potentially_smooth_ -= zeros;
}

void LatticeSiever::remove_sieved_factors1()
{
    //std::cout << "Removing sieved factors for ... f = " << f1_ << std::endl;
#ifdef RESIEVE1
    alg_pf_list_.set_end();
    alg_pf_list_.sort();
    PrimeFactor* pf_iter = alg_pf_list_.begin();

    PotentiallySmoothPoint* smooth_iter = head_psp_;

    while (pf_iter != alg_pf_list_.end() && smooth_iter)
    {
        if (pf_iter->offset_ + fixed_sieve_array_ > smooth_iter->ptr_)
        {
            smooth_iter = smooth_iter->next_;
        }
        else
        {
            if (pf_iter->offset_ + fixed_sieve_array_ == smooth_iter->ptr_)
            {
                smooth_iter->add_factor1(pf_iter->p_);
            }
            ++pf_iter;
        }
    }
    alg_pf_list_.reset();
#else
    PSPHashTable psp_hash_table(PotentiallySmoothPoint::get_sieve_ptr, PotentiallySmoothPoint::compare);
    psp_hash_table.load(head_psp_);
    alg_pf_list_.set_end();
    PrimeFactor* pf_iter = alg_pf_list_.begin();
    while (pf_iter != alg_pf_list_.end())
    {
        if (!sieve_bit_array_.isSet(pf_iter->offset_)) 
        {
            PotentiallySmoothPoint* psp = psp_hash_table.find((pf_iter->offset_ + fixed_sieve_array_));
            if (psp)
            {
                psp->add_factor1((pf_iter)->p_);
            }
        }
        ++pf_iter;
    }
#endif
}

void LatticeSiever::remove_sieved_factors2()
{
    //std::cout << "Removing sieved factors for ... f = " << f2_ << std::endl;
    rat_pf_list_.set_end();
    rat_pf_list_.sort();
    PrimeFactor* pf_iter = rat_pf_list_.begin();

    PotentiallySmoothPoint* smooth_iter = head_psp_;

    while (pf_iter != rat_pf_list_.end() && smooth_iter)
    {
        if (pf_iter->offset_ + fixed_sieve_array_ > smooth_iter->ptr_)
        {
            smooth_iter = smooth_iter->next_;
        }
        else
        {
            if (pf_iter->offset_ + fixed_sieve_array_ == smooth_iter->ptr_)
            {
                smooth_iter->add_factor2(pf_iter->p_);
            }
            ++pf_iter;
        }
    }
    rat_pf_list_.reset();
}

void LatticeSiever::generate_lattice(long int q, long int s)
{
    // from Lenstra's Lattice siever, code by Peter Montgomery, ratio <=> s ("skewedness")
    //double ratio = 23677.0; // Hardcoded for the moment
    double ratio = SKEWEDNESS_;
    //ratio = 1.0;

    long int x1 = q;
    long int y1 = 0L;
    long int x2 = s;
    long int y2 = 1L;
    long int x3 = 0L;
    long int y3 = 0L;

    while (1)
    {
        double dx1 = x1;
        double dx2 = x2;
        double dy1 = y1;
        double dy2 = y2;
        double dot12 = dx1*dx2 + ratio*ratio*dy1*dy2;
        double dot22 = dx2*dx2 + ratio*ratio*dy2*dy2;
        double dquot = dot12/dot22;
        if (fabs(dquot) < 0.500001) break;
        long int squot = nearest_integer(dquot);
        x3 = x1 - squot*x2;
        y3 = y1 - squot*y2;
        x1 = x2;
        y1 = y2;
        x2 = x3;
        y2 = y3;
    }
    if (x1*x1 + y1*y1 < x2*x2 + y2*y2)
    {
        x3 = x1;
        y3 = y1;
        x1 = x2;
        y1 = y2;
        x2 = x3;
        y2 = y3;
    }
    c1_.first = x1;
    c1_.second = y1;
    c2_.first = x2;
    c2_.second = y2;

    if (c1_.first > c2_.first)
    {
        std::swap(c1_, c2_);
    }

    if (c2_.first < 0)
    {
        c2_.first = -c2_.first;
        c2_.second = -c2_.second;
    }

    long int det = static_cast<long int>((long long)(c1_.first) * (long long)(c2_.second)
                                         - (long long)(c1_.second) * (long long)(c2_.first));

    if (det < 0)
    {
        c1_.first = -c1_.first;
        c1_.second = -c1_.second;
    }
}

/*
       A = a.a, B = b.b
       N = a.b
       r = | N / B | (nearest integer)
       tmp = b
       b = a - r b
       a = b
       T = A - 2r N + r^2 B
       A = B
       B = T
       new N is (a - rb).b = old N - r old B
*/
void LatticeSiever::generate_ef_lattice(int32_t p, int32_t r1,
                                        std::pair<int32_t, int32_t>& e2,
                                        std::pair<int32_t, int32_t>& e1)
{
    e1.second = p;
    e1.first = 0L;
    e2.second = r1;
    e2.first = 1L;
    double dx1 = e1.second;
    double dx2 = e2.second;
    double dy1 = e1.first;
    double dy2 = e2.first;
    double A = dx1 * dx1 + dy1 * dy1;
    double B = dx2 * dx2 + dy2 * dy2;
    double N = dx1 * dx2 + dy1 * dy2;
    long int r = nearest_integer(N / B);
    //long int r = (N / B < 0 ? -(long int)(0.5 - N / B) : (long int)(0.5 + N / B));
    //long int r = ::round(N / B);

    while (1)
    {
        double T = r * B;
        T -= N;
        T -= N;
        T *= r;
        T += A;
        std::pair<int32_t, int32_t> tmp = e2;
        e2.first = e1.first - r * e2.first;
        e2.second = e1.second - r * e2.second;
        e1 = tmp;
        if (T >= B)
        {
            if (debug_)
            {
                std::cerr << "generate_ef_lattice: p = " << p << ", r1 = " << r1 << ", e1 = (" << e1.first << "," << e1.second << "), e2 = (" << e2.first << "," << e2.second << ")" << std::endl;
                double A = (double)e1.first * e1.first + (double)e1.second * e1.second;
                double B = (double)e2.first * e2.first + (double)e2.second * e2.second;
                std::cerr << "A = " << A << ", B = " << B << std::endl;
            }
            return;
        }

        N -= r * B;
        A = B;
        B = T;
        r = nearest_integer(N / B);
        //r = (N / B < 0 ? -(long int)(0.5 - N / B) : (long int)(0.5 + N / B));
        //r = ::round(N / B);
    }
}

bool LatticeSiever::is_actually_smooth(FastVeryLong& remaining_quotient,
                                       FactorList& factor,
                                       long int B, long int L, long int LP, bool debug)
{
    if (remaining_quotient == 1L) return true;
    if (remaining_quotient.is_very_long()) return false;
    if (debug)
    {
        std::cerr << "is_actually_smooth : remaining_quotient = " << remaining_quotient << ", " << " L = " << L << std::endl;
    }
    const VeryLong one(1L);
    bool is_prime = remaining_quotient.is_probable_prime();
    if (is_prime &&
            (remaining_quotient <= B ||
             (LP > 0 && remaining_quotient <= L)))
    {
        if (debug)
        {
            std::cerr << "1. factor.push_back(" << remaining_quotient << "), [" << remaining_quotient.get_long() << "]" << std::endl;
        }
        factor.push_back(remaining_quotient.get_long());
        return true;
    }
    // Here if
    //    remaining_quotient is > B and LP > 0 and
    //    remaining_quotient > L or remaining_quotient is composite
    // But if remaining_quotient is prime, then it's a prime > L, but we've already excluded
    // this case so we must have
    //    remaining_quotient is > B and LP > 0 and remaining_quotient is composite
    // note: remaining_quotient could still be <= L
    if (remaining_quotient <= L)
    {
        if (debug)
        {
            std::cerr << "1. Calling factorise on " << remaining_quotient << std::endl;
        }
        std::vector<long int> fac;
        long int rq_l = remaining_quotient.get_long();
        if (!VeryLong::factorise_no_trial_l(rq_l, &fac)) return false;
        if (debug)
        {
            std::cerr << "1. done" << std::endl;
        }
        for (auto& f: fac)
        {
            if (debug)
            {
                std::cerr << "2. factor.push_back(f=" << f << ")" << std::endl;
            }
            factor.push_back(f);
            rq_l /= f;
        }
        if (rq_l != 1L)
        {
            if (debug)
            {
                std::cerr << "3. factor.push_back(rq_ll=" << rq_l << ")" << std::endl;
            }
            factor.push_back((long int)rq_l);
        }
        return true;
    }

    // Check if it's possible to be smooth still
    if (remaining_quotient.get_long_long() > (long long)L * (long long)L)
        return false;
    std::vector<long int> fac;
    if (debug)
    {
        std::cerr << "2. Calling factorise on " << remaining_quotient << std::endl;
    }
    if (!remaining_quotient.factorise_no_trial(&fac)) return false;

    for (auto& f: fac)
    {
        if (f > L)
        {
            if (debug)
            {
                std::cerr << std::endl;
                std::cerr << "not smooth : " << L << " < " << f << std::endl;
            }
            return false;
        }
        if (debug)
        {
            std::cerr << "4. factor.push_back(f=" << f << ")" << std::endl;
        }
        factor.push_back(f);
    }
    return true;
}

void LatticeSiever::print_relation(long int c, long int d,
                                   FactorList& factors1,
                                   FactorList& factors2)
{
    const VeryLong zero(0L);
    const VeryLong one(1L);
    if (debug_)
    {
        std::cerr << "(c,d) = (" << c << "," << d << ") => (a,b) = ";
    }
    VeryLong a(c);
    a *= c1_.first;
    VeryLong tmp(d);
    tmp *= c2_.first;
    a += tmp;
    VeryLong b(c);
    b *= c1_.second;
    tmp = d;
    tmp *= c2_.second;
    b += tmp;
    if (b < zero)
    {
        b = -b;
        a = -a;
    }
    if (debug_)
    {
        std::cerr << "(" << a << "," << b << ")" << std::endl;
    }
    if (gcd<VeryLong>(a, b) == one)
    {
        if (!relfile_)
        {
            //cout << "Opening RELATION_FILE (" << relation_file_ << ") ...";
            relfile_ = new std::fstream(relation_file_.c_str(), std::ios::out|std::ios::app);
            //cout << " opened" << endl;
        }
        *relfile_ << a << " " << b << " :";
        std::sort(factors1.begin(), factors1.end());
        for (size_t i = 0; i < factors1.size(); i++)
        {
            long int p = factors1[i];
            if (p >= B1_ && !alg_factor_base_->exists_extra(p))
            {
                std::vector<LongModular> roots;
                find_roots_mod_p<VeryLong, long int, LongModular>(f1_, p, roots);
                alg_factor_base_->add_extra(p, roots);
            }
            int done = 0;
            for (FactorBase::a_const_root_iterator iter = alg_factor_base_->begin(p);
                    !done && iter != alg_factor_base_->end(p);
                    ++iter)
            {
                // find the r such that a = br mod p or r = p
                long int r = *iter;
                if (r == p) done = 1;
                else
                {
                    long long int xll = b % p;
                    xll *= r;
                    //long int x = xll % p;  replace by modasm
                    long int x = modasm(xll, p);
                    //if (x < 0) x += p; replace by modasm
                    long int y = a % p;
                    if (y < 0) y += p;
                    if (x == y) done = 1;
                }
                if (done) *relfile_ << " " << p << "/" << r;
            }
        }
        *relfile_ << " :";
        std::sort(factors2.begin(), factors2.end());
        for (size_t i = 0; i < factors2.size(); i++)
        {
            *relfile_ << " " << factors2[i];
        }
        if (factors2.size() == 0) *relfile_ << " 1";
        *relfile_ << " :";
        *relfile_ << std::endl;
    }
}

int LatticeSiever::check_for_remaining_relations()
{
    int relations = 0;
    for (PotentiallySmoothPoint* smooth_iter = head_psp_;
            smooth_iter;
            smooth_iter = smooth_iter->next_)
    {
        if (is_actually_smooth(smooth_iter->partial2_.remaining_quotient_,
                               smooth_iter->partial2_.factor_, B2_, L2_, LP2_, debug_))
        {
            if (debug_)
            {
                std::cerr << "(c,d) = (" << smooth_iter->c_ << "," << smooth_iter->d_ << "), partial2_.remaining_quotient = " << smooth_iter->partial2_.remaining_quotient_ << ", is actually smooth" << std::endl;

            }
            if (is_actually_smooth(smooth_iter->partial1_.remaining_quotient_,
                                   smooth_iter->partial1_.factor_, B1_, L1_, LP1_, debug_))
            {
                if (debug_)
                {
                    std::cerr << "(c,d) = (" << smooth_iter->c_ << "," << smooth_iter->d_ << "), partial1_.remaining_quotient = " << smooth_iter->partial1_.remaining_quotient_ << ", is actually smooth" << std::endl;
                }
                print_relation(smooth_iter->c_, smooth_iter->d_, smooth_iter->partial1_.factor_, smooth_iter->partial2_.factor_);
                relations++;
            }
            else
            {
                if (debug_)
                {
                    std::cerr << "(c,d) = (" << smooth_iter->c_ << "," << smooth_iter->d_ << "), partial1_.remaining_quotient = " << smooth_iter->partial1_.remaining_quotient_ << ", is not actually smooth" << std::endl;
                }
            }
        }
        else
        {
            if (debug_)
            {
                std::cerr << "(c,d) = (" << smooth_iter->c_ << "," << smooth_iter->d_ << "), partial2_.remaining_quotient = " << smooth_iter->partial2_.remaining_quotient_ << ", is not actually smooth" << std::endl;
            }
        }

    }
    return relations;
}

inline void LatticeSiever::sieve1(FactorBase::a_iterator iter, long int r1)
{
    std::pair<int32_t, int32_t> e1;
    std::pair<int32_t, int32_t> e2;
    generate_ef_lattice(iter->get_p(), r1, e1, e2);
    int32_t e12 = e1.first + (e1.second << c_span_bits);
    int32_t e22 = e2.first + (e2.second << c_span_bits);

    Parallelogram E_region(c_region_, e1, e2);
    int32_t e_min = static_cast<int32_t>(std::ceil(E_region.min_x()));
    int32_t e_max = static_cast<int32_t>(std::floor(E_region.max_x()));
    int32_t f_min = 0L;
    int32_t f_max = 0L;
    int32_t e = e_min;

    while (e <= e_max)
    {
        if (E_region.y_limits1(e, f_min, f_max))
        {
            uint32_t ptr = e * e12 + (f_min - 1) * e22 - min_c;
            int32_t f_span = f_max - f_min + 1;
#ifdef RESIEVE1
            sieveCache_.add1(ptr, f_span, e22, iter);
#else
            sieveCache_.add(ptr, f_span, e22, iter);
#endif
        }
        ++e;
    }
}

#ifdef RESIEVE1
inline void LatticeSiever::sieve1_again(FactorBase::a_iterator iter, long int r1)
{
    std::pair<int32_t, int32_t> e1;
    std::pair<int32_t, int32_t> e2;
    generate_ef_lattice(iter->get_p(), r1, e1, e2);
    int32_t e12 = e1.first + (e1.second << c_span_bits);
    int32_t e22 = e2.first + (e2.second << c_span_bits);

    Parallelogram E_region(c_region_, e1, e2);
    int32_t e_min = static_cast<int32_t>(std::ceil(E_region.min_x()));
    int32_t e_max = static_cast<int32_t>(std::floor(E_region.max_x()));
    int32_t f_min = 0L;
    int32_t f_max = 0L;
    int32_t e = e_min;

    while (e <= e_max)
    {
        if (E_region.y_limits1(e, f_min, f_max))
        {
            uint32_t ptr = e * e12 + (f_min - 1) * e22 - min_c;
            int32_t f_span = f_max - f_min + 1;
            sieveCache_.add1_again(ptr, f_span, e22, iter);
        }
        ++e;
    }
}
#endif

inline void LatticeSiever::sieve2(FactorBase::a_iterator iter, long int r1)
{
    std::pair<int32_t, int32_t> e1;
    std::pair<int32_t, int32_t> e2;
    generate_ef_lattice(iter->get_p(), r1, e1, e2);
    int32_t e12 = e1.first + (e1.second << c_span_bits);
    int32_t e22 = e2.first + (e2.second << c_span_bits); 

    Parallelogram E_region(c_region_, e1, e2);
    int32_t e_min = static_cast<int32_t>(std::ceil(E_region.min_x()));
    int32_t e_max = static_cast<int32_t>(std::floor(E_region.max_x()));
    int32_t f_min = 0L;
    int32_t f_max = 0L;
    int32_t e = e_min;

    while (e <= e_max)
    {
        if (E_region.y_limits1(e, f_min, f_max))
        {
            uint32_t ptr = e * e12 + (f_min - 1) * e22 - min_c;
            int32_t f_span = f_max - f_min + 1;
            sieveCache_.add(ptr, f_span, e22, iter);
        }
        ++e;
    }
}

/*
 *   (a,b) = c c1 + d c2
 *         = (c c1.first + d c2.first, c c1.second + d c2.second)
 * If (c,d) gives a point in the (q,s) lattice then for it also to be in the (p,r) lattice we must have
 *   c c1.first + d c2.first = r (c c1.second + d c2.second) mod p     (1)
 * i.e.
 *   c (c1.first - r c1.second) = d (r c2.second - c2.first) mod p     (2)
 * If
 *   c1.first = r c1.second mod p
 * then c1 is also in the (p,r) lattice so can be one of the basis vectors, and we have
 *   d (r c2.second - c2.first) = 0 mod p
 * If
 *   r c2.second != c2.first mod p
 * then we must have d = 0 mod p
 * Similarly, if
 *   c2.first = r c2.second mod p
 * then c2 can be one of the basis vectors, and we must have c = 0 mod p.
 * Otherwise, we have
 *   d = c (r c2.second - c2.first)^-1 (c1.first - r c1.second) mod p
 * or
 *   d = c r' mod p
 * where
 *   r' = (r c2.second - c2.first)^-1 (c1.first - r c1.second) mod p
 * If we calculate a basis for the (p, r') lattice, {e1, e2} say, then for any
 *   (c,d) = e e1 + f e2
 * we will have
 *   c = d r' mod p
 * by definition of e1 and d2.
 * If
 *   (a,b) = c c1 + d c2                                               (3)
 * then (a,b) will satisfy
 *   a = b s mod q
 * by definition of c1 and c2 and in addition
 *   a = c c1.first + d c2.first
 *     = r ( c c1.second + d c2.second) mod p [from (1)]
 *     = r b mod p [ from (3) ]
 */

void LatticeSiever::sieve_by_vectors1()
{
    if (debug_)
    {
        std::cerr << "sieve_by_vectors1 : min_c = " << min_c << std::endl;
        std::cerr << "c1 = (" << c1_.first << "," << c1_.second << "), c2 = (" << c2_.first << "," << c2_.second << ")" << std::endl;
    }
#ifndef RESIEVE1
    SieveCacheItem::set_pf_list(&alg_pf_list_);
#endif
    FactorBase::a_iterator iter = alg_factor_base_->begin();
    FactorBase::a_iterator enditer = alg_factor_base_->end();

    for (; iter != enditer; ++iter)
    {
        if (iter->get_p() < SMALL_PRIME_BOUND1_) continue;
        if (iter->get_p() > B1_) break;
        for (FactorBase::a_const_root_iterator root_info_iter = alg_factor_base_->begin(iter);
                root_info_iter != alg_factor_base_->end(iter);
                ++root_info_iter)
        {
            long int p = iter->get_p();
            long int r = *root_info_iter;
            if (p == r) continue;
            // now find short vectors in the sub-lattice of the (q,s) lattice
            // which intersects the (p,r) lattice.

            long long int Q_ll = c1_.first - (long long)c1_.second * r;
            long int Q = modasm(Q_ll, p);
            if (Q)
            {
                long long int R_ll = (long long)c2_.second * r - c2_.first;
                long int R = modasm(R_ll, p);
                // we are interested in (c,d) such that c Q = d R mod p
                // or c Q R^-1 = d mod p or cr' = d mod p
                if (R)
                {
                    long int R_inv = inverse<long int>((long int)R, p);
                    // r' = (r c2_.second - c2_.first)^-1 (c1_.first - r c1_.second) mod p
                    long int r1 = 0L;
                    mulmodasm2(Q, R_inv, p, r1);
                    sieve1(iter, r1);
                }
            }
        }
    }
#ifdef RESIEVE1
    sieveCache_.dump(false);
#else
    sieveCache_.dump(true);
#endif
}

#ifdef RESIEVE1
void LatticeSiever::sieve_by_vectors1_again()
{
    if (debug_)
    {
        std::cerr << "sieve_by_vectors1_again : min_c = " << min_c << std::endl;
        std::cerr << "c1 = (" << c1_.first << "," << c1_.second << "), c2 = (" << c2_.first << "," << c2_.second << ")" << std::endl;
    }
    SieveCacheItem::set_pf_list(&alg_pf_list_);
    FactorBase::a_iterator iter = alg_factor_base_->begin();
    FactorBase::a_iterator enditer = alg_factor_base_->end();

    for (; iter != enditer; ++iter)
    {
        if (iter->get_p() < SMALL_PRIME_BOUND1_) continue;
        for (FactorBase::a_const_root_iterator root_info_iter = alg_factor_base_->begin(iter);
                root_info_iter != alg_factor_base_->end(iter);
                ++root_info_iter)
        {
            long int p = iter->get_p();
            long int r = *root_info_iter;
            if (p == r) continue;
            // now find short vectors in the sub-lattice of the (q,s) lattice
            // which intersects the (p,r) lattice.

            long long int Q_ll = c1_.first - (long long)c1_.second * r;
            long int Q = modasm(Q_ll, p);
            if (Q)
            {
                long long int R_ll = (long long)c2_.second * r - c2_.first;
                long int R = modasm(R_ll, p);
                // we are interested in (c,d) such that c Q = d R mod p
                // or c Q R^-1 = d mod p or cr' = d mod p
                if (R)
                {
                    long int R_inv = inverse<long int>((long int)R, p);
                    // r' = (r c2_.second - c2_.first)^-1 (c1_.first - r c1_.second) mod p
                    long int r1 = 0L;
                    mulmodasm2(Q, R_inv, p, r1);
                    sieve1_again(iter, r1);
                }
            }
        }
    }
    sieveCache_.dump(true);
}
#endif

void LatticeSiever::sieve_by_vectors2()
{
    if (debug_)
    {
        std::cerr << "sieve_by_vectors2 : min_c = " << min_c << std::endl;
        std::cerr << "c1 = (" << c1_.first << "," << c1_.second << "), c2 = (" << c2_.first << "," << c2_.second << ")" << std::endl;
    }
    SieveCacheItem::set_pf_list(&rat_pf_list_);
    FactorBase::a_iterator iter = rat_factor_base_->begin();
    FactorBase::a_iterator enditer = rat_factor_base_->end();

    for (; iter != enditer; ++iter)
    {
        long int p = iter->get_p();
        if (p < SMALL_PRIME_BOUND2_) continue;
        if (p > B2_) break;
#if 0
        for (FactorBase::a_const_root_iterator* root_info_iter = rat_factor_base_->begin(iter);
                root_info_iter != rat_factor_base_->end(iter);
                ++root_info_iter)
        {
            long int r = root_info_iter->r;
#else
        {
            long int r = *(rat_factor_base_->begin(iter));
#endif
            if (p == r) continue;
            // now find short vectors in the sub-lattice of the (q,s) lattice
            // which intersects the (p,r) lattice.

            long long int Q_ll = c1_.first - (long long)c1_.second * r;
            long int Q = modasm(Q_ll, p);
            if (Q)
            {
                long long int R_ll = (long long)c2_.second * r - c2_.first;
                long int R = modasm(R_ll, p);
                // we are interested in (c,d) such that c Q = d R mod p
                // or c Q R^-1 = d mod p or cr' = d mod p
                if (R)
                {
                    long int R_inv = inverse<long int>((long int)R, p);
                    // r' = (r c2_.second - c2_.first)^-1 (c1_.first - r c1_.second) mod p
                    long int r1 = 0L;
                    mulmodasm2(Q, R_inv, p, r1);
                    sieve1(iter, r1);
                }
            }
        }
    }
    sieveCache_.dump();
}

bool LatticeSiever::allocate_c_d_region()
{
    memset(fixed_sieve_array_, (SIEVE_TYPE)0, fixed_sieve_array_size * sizeof(SIEVE_TYPE));
    sieve_bit_array_.clear();
    // The offset of (c,d) from start of fixed_sieve_array_ is
    //
    //    c - min_c + (max_c - min_c + 1) * d
    //

    // we are only interested in points (c,d) with gcd(c,d) = 1
    // so for each prime which could divide a c in the region
    // mark each point (c,d) with the prime dividing c and d.
    VeryLong::generate_prime_table();
    long int p = VeryLong::firstPrime();
    long int B = std::max(std::abs(min_c), std::abs(max_c));
    while (p < B)
    {
        long int z = min_d % p;
        if (z <= 0) z += p;
        long int d = min_d - z + p;
        while (d <= max_d)
        {
            z = min_c % p;
            if (z <= 0) z += p;
            long int c = min_c - z + p;
            SIEVE_TYPE* ptr = fixed_sieve_array_ + c - min_c + d * (max_c - min_c + 1);
            while (c <= max_c)
            {
                //*ptr = -128;
                sieve_bit_array_.set(ptr - fixed_sieve_array_);
                ptr += p;
                c += p;
            }
            d += p;
        }

        p = VeryLong::nextPrime();
    }

    return true;
}

//#define RANGE_CHECK 1
#ifdef RANGE_CHECK
bool LatticeSiever::in_range(long int c, long int d)
{
    if (c < min_c || c > max_c)
    {
        std::cerr << "in_range : c = " << c << ", min_c = " << min_c << ", max_c = " << max_c << std::endl;
        return false;
    }
    if (d < min_d || d > max_d)
    {
        std::cerr << "in_range : c = " << c << ",d = " << d << ", min_d = " << min_d << ", max_d = " << max_d << std::endl;
        return false;
    }
    return true;
}
#endif

void LatticeSiever::sieve_by_vectors(long int q, long int s)
{
    if (!allocate_c_d_region())
        return;
    Timer sieving_timer;
    sieving_timer.start();
    int relations = 0;
    rat_pf_list_.reset();
    alg_pf_list_.reset();
    if (!potentially_smooth_point_)
    {
        potentially_smooth_point_ = new PotentiallySmoothPoint[max_potentially_smooth];
    }
    head_psp_ = potentially_smooth_point_;
    number_potentially_smooth_ = 0;
    timer_.start("sieve by vectors 1");
    sieve_by_vectors1();
    timer_.stop();

    timer_.start("check interval 1");
    int number_of_potentially_alg_smooth = check_interval1(q);
    timer_.stop();
    if (verbose())
    {
        std::cout << number_of_potentially_alg_smooth << " -> " << std::flush;
    }

    timer_.start("sieve by vectors 2");
    sieve_by_vectors2();
    timer_.stop();

    timer_.start("check interval 2");
    check_interval2();
    timer_.stop();

    timer_.start("remove factors for rational");
    remove_sieved_factors2();
    timer_.stop();
    timer_.start("eliminate rational");
    eliminate2();
    timer_.stop();
    if (verbose())
    {
        std::cout << number_potentially_smooth_ << " -> " << std::flush;
    }

    if (number_potentially_smooth_ < max_potentially_smooth)
    {
#ifdef RESIEVE1
        sieve_by_vectors1_again();
#endif
        // Now remove the known factors for the algebraic norms
        timer_.start("remove factors for algebraic");
        remove_sieved_factors1();
        timer_.stop();
        timer_.start("eliminate algebraic");
        eliminate1(q);
        timer_.stop();

        // We should now have the smooth relations which remain
        timer_.start("check relations");
        relations = check_for_remaining_relations();
        if (verbose())
        {
            std::cout << relations;
        }
        timer_.stop();
    }
    double sieving_time = sieving_timer.stop();
    //total_sieving_time_ += sieving_time;
    total_sieving_time_ = sieving_timer.total_elapsed();
    total_relations_ += relations;
    double average_relations_per_second = (double)relations / sieving_time;
    double average_relations_per_hour = average_relations_per_second * 60 * 60;
    double average_relations_per_day = average_relations_per_hour * 24;
    double running_average_relations_per_second = (double)total_relations_ / total_sieving_time_;
    double running_average_relations_per_hour = running_average_relations_per_second * 60 * 60;
    double running_average_relations_per_day = running_average_relations_per_hour * 24;
    static double prev_running_average_relations_per_day = 0;
    static double prev_total_sieving_time = 0;
    if (verbose())
    {
        std::cout << " (" << average_relations_per_second << "," << (int)average_relations_per_hour << "," << (int)average_relations_per_day << "),";
        std::cout << " (" << running_average_relations_per_second << "," << (int)running_average_relations_per_hour << "," << (int)running_average_relations_per_day << "),";
    }
    if (prev_running_average_relations_per_day > 0 && 
        running_average_relations_per_day / prev_running_average_relations_per_day < 0.5)
    {
        std::cerr << "!!!! relation rate just reduced" << std::endl;
        std::cerr << "prev_total_sieving_time = " << prev_total_sieving_time << std::endl;
        std::cerr << "total_sieving_time_ = " << total_sieving_time_ << std::endl;
        sieving_timer.reset();
        total_relations_ = 0;
    }
    prev_running_average_relations_per_day = running_average_relations_per_day;
    prev_total_sieving_time = total_sieving_time_;
   
   if (verbose()) 
   {
       std::cout << std::endl;
   }
}

// Lattice sieving for NFS
bool LatticeSiever::sieve(long int q)
{
    static int firsttime = 1;
    if (firsttime)
    {
        firsttime = 0;
        VeryLong::generate_prime_table();
    }
    // Lattice sieve for NFS, based on description in
    // "The Development of the Number Field Sieve", Lenstra & Lenstra (eds.) 1993

//   cout << "NFS Lattice Sieve:" << endl;
//   cout << "f1 = " << f1_ << endl;
//   cout << "f2 = " << f2_ << endl;

    VeryLong qq(q);
    const VeryLong one(1L);

    // find roots of f1_ mod q
    std::vector<LongModular> q_roots;
    find_roots_mod_p<VeryLong, long int, LongModular>(f1_, q, q_roots);
    try
    {
        for (size_t i = 0; i < q_roots.size(); i++)
        {
            long int s = q_roots[i].get_long();
            if (debug_)
            {
                std::cerr << "q = " << q << std::endl;
                std::cerr << "s = " << s << std::endl;
            }
            // Generate lattice of points (a,b) such that a = bs mod q
            // Note that f1_(a,b) = 0 mod q for all points in this lattice
            // thus making it more likely that f1_ will be smooth.
            // To generate the lattice, we look for a pair of short vectors
            // which form a basis for the lattice.
            generate_lattice(q, s);
            if (debug_)
            {
                std::cerr << "c1 = (" << c1_.first << "," << c1_.second << ")" << std::endl;
                std::cerr << "c2 = (" << c2_.first << "," << c2_.second << ")" << std::endl;
            }
            sieve_by_vectors(q, s);
        }
        if (q_roots.size()) timer_.summary();
    }
    catch (...)
    {
        std::cerr << "Unknown exception";
    }
    if (total_relations_ > max_total_relations)
    {
        return false;
    }
    else
    {
        return true;
    }
}
