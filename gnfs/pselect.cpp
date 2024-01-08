#include "Polynomial.h"
#include "Config.h"
#include "Combinations.h"
#include "LongModular.h"
#include "VeryLongModular.h"
#include "discriminant.h"
#include "dickman.h"
#include "PriorityQueue.h"
#include "MPFloat.h"
#include "PolynomialOptimizer.h"
#include <time.h>
#include <iomanip>
#include <fstream>
#include "timings.h"
#include "lip.h"
#include "pselect.h"
#include "gcd.h"
#include <math.h>
#include <set>
#include <cmath>

#define CHECK 1

namespace
{
Skewed_selection_config Skewed_config("skewed.cfg");

struct Kleinjung_poly_info
{
    Kleinjung_poly_info(const VeryLong& a, const VeryLong& b, const Polynomial<VeryLong>& fm)
        : a_(a), b_(b), fm_(fm)
    {
        //als_ = average_log_size(fm, 1000);
    }
    Kleinjung_poly_info(const Kleinjung_poly_info& pi)
    {
        a_ = pi.a_;
        b_ = pi.b_;
        fm_ = pi.fm_;
        als_ = pi.als_;
    }
    Kleinjung_poly_info& operator=(const Kleinjung_poly_info& pi)
    {
        if (&pi != this)
        {
            a_ = pi.a_;
            b_ = pi.b_;
            fm_ = pi.fm_;
            als_ = pi.als_;
        }
        return *this;
    }
    bool operator<(const Kleinjung_poly_info& pi) const
    {
        //return (als_ < pi.als_);
        int d = fm_.deg();
        //VeryLong x = fm_.coefficient(d - 1) * fm_.coefficient(d - 2);
        //VeryLong y = pi.fm_.coefficient(d - 1) * pi.fm_.coefficient(d - 2);
        VeryLong x = fm_.coefficient(d - 2);
        VeryLong y = pi.fm_.coefficient(d - 2);
        return (abs(x) < abs(y));
//      return (abs(fm_.coefficient(d - 2)) < abs(pi.fm_.coefficient(d - 2)));
    }

    VeryLong a_;
    VeryLong b_;
    Polynomial<VeryLong> fm_;
    double als_;
};
//---------------------------------------------------------

Polynomial<VeryLong> adjust_base_m_polynomial(const Polynomial<VeryLong>& poly, const VeryLong& m)
{
    const VeryLong one(1L);
    const VeryLong two(2L);
    // adjust coefficients of poly, so that f(m) = 0 mod N still holds
    // but coefficients are smaller
    Polynomial<VeryLong> adjusted_poly = poly;
    int i = 0;
    VeryLong half_m = m / two;
    for (i = 0; i < poly.deg(); i++)
    {
        VeryLong a = adjusted_poly.coefficient(i);
        while (a > half_m)
        {
            VeryLong k = a / m + one;
            adjusted_poly.set_coefficient(i, a - k * m);
            adjusted_poly.set_coefficient(i + 1, adjusted_poly.coefficient(i + 1) + k);
            a = adjusted_poly.coefficient(i);
        }
    }
//   std::cout << "Original polynomial = " << poly << std::endl;
//   std::cout << "Adjusted polynomial = " << adjusted_poly << std::endl;
    return adjusted_poly;
}
};
namespace
{
Polynomial<VeryLong> adjust_root_properties_orig(const Polynomial<VeryLong>& min_poly,
        const VeryLong& m,
        VeryLong& s,
        double average_log_size,
        std::vector<PolynomialOptimizer::Poly_info>& poly_list)
{
    const VeryLong one(1L);
    std::cout << "Sieving for best root properties ..." << std::endl;
    const long int MAX_J0 = Skewed_config.MAX_J0();
    const long int MAX_J1 = Skewed_config.MAX_J1();
    const long int MAX_SMALL_PRIME = Skewed_config.MAX_SMALL_PRIME();
    const unsigned long int MAX_P_K = MAX_SMALL_PRIME;
    const double target_E_F = Skewed_config.PRINTING_BOUND();
    double alpha_cutoff = target_E_F - average_log_size;
    if (alpha_cutoff < Skewed_config.ALPHA_CUTOFF())
    {
        alpha_cutoff = Skewed_config.ALPHA_CUTOFF();
    }
    std::cout << "alpha cutoff = " << alpha_cutoff << std::endl;

    static short* cont_array_data = 0;
    static long int cont_array_data_size = 0;
    if (cont_array_data == 0)
    {
        cont_array_data_size = (2 * MAX_J1 + 1) * (2 * MAX_J0 + 1) ;
        cont_array_data = new short [ cont_array_data_size ];
    }

    static short** cont_array = 0;
    if (cont_array == 0)
    {
        cont_array = new short* [ 2 * MAX_J1 + 1 ];
        for (long int j1 = 0; j1 < 2 * MAX_J1 + 1; j1++)
        {
            cont_array[j1] = cont_array_data + j1 * (2 * MAX_J0 + 1);
        }
    }
    static short* j0_array = 0;
    if (j0_array == 0)
    {
        j0_array = new short [ MAX_J0 ];
    }
    unsigned long int p_k;  // p^k
    //long int k = 0;

    memset((char*)cont_array_data, 0, cont_array_data_size * sizeof(short));
    long int p = zpnextb(2);
    VeryLong leading_coefficient = min_poly.coefficient(min_poly.deg());
    time_t start = time(0);
    while (p < MAX_SMALL_PRIME)
    {
        p_k = p;
        //k = 1;
        int projective_root = (leading_coefficient % p == 0L);
        double logp = log((double)p);
        while (p_k < MAX_P_K)
        {
            p_k *= p;
            //k++;
        }
        p_k /= p;
        //k--;
        //cout << "p^k = " << p_k << std::endl;

        double prob_here = (double)p / (p_k * (p + 1.0));

        LongModular::set_default_modulus(p_k);
        VeryLong mm = m % (long)p_k;
        VeryLong m_mod_p = m % p;
        LongModular m_lm(mm.get_long());
        Polynomial<LongModular> f = convert_to_F_p<VeryLong, long int, LongModular>(min_poly, p_k);
        int degree = f.deg();

        unsigned long int max_l = p_k;
        if (projective_root) max_l += p_k / p;

        double initial_evaluation = logp / (p - 1.0);
        short initial_evaluation_s = (short)floor(initial_evaluation * 1000.0 + 0.5);

        for (long int j1 = -MAX_J1; j1 <= MAX_J1; j1++)
        {
            memset((char*)j0_array, 0, sizeof(short)*MAX_J0);
            for (unsigned long int jj0 = 0; jj0 < p_k; jj0++) j0_array[jj0] = -initial_evaluation_s;
            // use finite differences to calculate f(l) for consecutive values of l
            // assume at least a cubic f
            LongModular f_value = f.coefficient(0);
            LongModular d1 = f.coefficient(1) + f.coefficient(2) + f.coefficient(3);
            LongModular d2 = 2L*f.coefficient(2) + 6L*f.coefficient(3);
            LongModular d3 = 6L*f.coefficient(3);
            LongModular d4 = 0L;
            if (degree > 3)
            {
                d1 += f.coefficient(4);
                d2 += 14L*f.coefficient(4);
                d3 += 36L*f.coefficient(4);
                d4 += 24L*f.coefficient(4);
            }
            if (degree > 4)
            {
                d1 += f.coefficient(5);
                d2 += 30L*f.coefficient(5);
                d3 += 150L*f.coefficient(5);
                d4 += 240L*f.coefficient(5);
            }

            for (unsigned long int l = 0; l < max_l; l++)
            {
                long int con0;
                long int con1;
                long int l_minus_m = (p_k + l - m_lm.get_long()) % p_k; // l - m (mod p^k)
                if (l < p_k)
                {
                    // possible non-projective
                    con0 = (f_value.get_long() + (p_k + j1) * l * l_minus_m) % p_k;
                    con1 = l_minus_m;
                    // next value of f
                    if (l >= 4 && degree > 4) d4 += 120L*f.coefficient(5);
                    if (l >= 3) d3 += d4;
                    if (l >= 2) d2 += d3;
                    if (l >= 1) d1 += d2;
                    f_value += d1;
                }
                else
                {
                    // possible projective
                    // evaluate homogeneous polys
                    LongModular y(static_cast<long int>(p_k + p * (l - p_k)));
                    LongModular con0_lm = f.evaluate_homogeneous(LongModular(1L), y);
                    LongModular con1_lm = LongModular(1L) - m_lm * y;
                    for (int i = 0; i < degree - 2; i++) con1_lm *= y;
                    con0_lm += LongModular(static_cast<long int>(p_k + j1)) * con1_lm;
                    con1_lm *= y;
                    con0 = con0_lm.get_long();
                    con1 = con1_lm.get_long();
                }

                if (con0 < 0) con0 += p_k;
                if (con1 < 0) con1 += p_k;

                long int con1gcd = gcd<long int>(con1, p_k);
                if (con1gcd == 0L) con1gcd = p_k;
                long int con10gcd = gcd<long int>(con0, con1gcd);
                if (con10gcd == 0L) con10gcd = gcd<long int>(p_k, con1gcd);
                long int con0left = con0;
                long int con1left = con1;
                long int p_k_left = p_k;

                if (con10gcd > 1)
                {
                    con0left /= con10gcd;
                    con1left /= con10gcd;
                    p_k_left /= con10gcd;
                    short addon = (short)floor(log((double)con10gcd) * prob_here * 1000.0 + 0.5);
                    for (unsigned long int j0 = 0; j0 < p_k; j0++)
                    {
                        j0_array[j0] += addon;
                    }
                }
                if (p_k_left > 1 && con1gcd == con10gcd)
                {
                    LongModular::set_default_modulus(p_k_left);
                    LongModular j0_start = LongModular(con0left) / LongModular(con1left);
                    double addon = logp * prob_here;
                    long int p_q = 1;
                    for (p_q = 1; p_q <= p_k_left; p_q *= p)
                    {
                        double addonj = 0.0;
                        if (p_q > 1) addonj += addon;
                        if (p_q == p_k_left) addonj += addon/(p - 1.0);
                        if (addonj > 0.0)
                        {
                            short addonj_s = (short)floor(addon * 1000.0 + 0.5);
                            unsigned long int j0 = j0_start.get_long() % p_q;
                            //if (j0 < 0) j0 += p_q;
                            while (j0 < p_k)
                            {
                                j0_array[j0] += addonj_s;
                                j0 += p_q;
                            }
                        }
                    }
                    LongModular::set_default_modulus(p_k);
                }

            }
            // replicate through J0-space

            short* src_ptr = &j0_array [0];
            short* src_end_ptr = &j0_array[p_k];
            short* tar_ptr = &cont_array [MAX_J1 + j1] [MAX_J0]; // j1, j0 = 0
            short* tar_end_ptr = &cont_array[MAX_J1 + j1][MAX_J0 + MAX_J0];
            short* ptr = src_ptr;
            while (tar_ptr != tar_end_ptr + 1)
            {
                // scale, round and store in short
                *tar_ptr -= *ptr;
                ptr++;
                if (ptr == src_end_ptr) ptr = src_ptr; // wrap
                tar_ptr++;
            }
            tar_ptr = &cont_array [MAX_J1 + j1] [MAX_J0 - 1]; // j1, j0 = -1
            tar_end_ptr = &cont_array [MAX_J1 + j1] [0];
            ptr = src_end_ptr - 1;
            while (tar_ptr != tar_end_ptr - 1)
            {
                *tar_ptr -= *ptr;
                ptr--;
                if (ptr == src_ptr - 1) ptr = src_end_ptr - 1; // wrap
                tar_ptr--;
            }
        }
        p = zpnext();
    }
    time_t end = time(0);
    std::cout << "Time for sieving (J0 = " << MAX_J0 << ", J1 = " << MAX_J1 << ") = " << end - start << std::endl;

    double best = 0.0;
    long int j0_best = 0;
    long int j1_best = 0;
    Polynomial<VeryLong> best_F;
    double best_alpha = 0.0;
    double best_E = Skewed_config.PRINTING_BOUND() + 1;
    double printing_bound = Skewed_config.PRINTING_BOUND() + 2;
    // calculate scaling for alpha estimate
    long int j0 = 0;
    long int j1 = 0;

    for (j1 = -MAX_J1; j1 <= MAX_J1; j1++)
    {
        for (j0 = -MAX_J0; j0 <= MAX_J0; j0++)
        {
            double cont = ((double)(cont_array [MAX_J1 + j1] [MAX_J0 + j0])*0.001);
            if (cont < best)
            {
                best = cont;
                j0_best = j0;
                j1_best = j1;
            }
            if (cont < alpha_cutoff)
            {
                // worth a look

                std::vector<VeryLong> c;
                c.resize(3);
                c[0] = VeryLong(j0) * m;
                c[1] = -one * (VeryLong(j0) + m * VeryLong(j1));
                c[2] = j1;
                Polynomial<VeryLong> F = min_poly + Polynomial<VeryLong>(c);

                s = PolynomialOptimizer::minimize_I_over_s(Polynomial<VeryLong>::convert_to_double<double>(F), 10000000.0);
                double I_F_S = PolynomialOptimizer::average_log_size(F, s);
                if (I_F_S + cont < printing_bound)
                {
                    poly_list.push_back(PolynomialOptimizer::Poly_info(F, s, I_F_S, cont));
                }

                // since calculating alpha is fairly slow, only
                // do it if we have a hope of beating the best E so far
                if (I_F_S + cont < best_E)
                {
                    //double alpha = alpha_F(F, 100, 50);
                    double alpha = cont;
                    double E = I_F_S + alpha;
                    //cerr << alpha << "   " << cont << "   " << alpha - cont << std::endl;
                    if (E < best_E)
                    {
                        best_alpha = alpha;
                        best_F = F;
                        best_E = E;
                        std::cout << "j0 = " << j0 << " and j1 = " << j1 << ", cont = " << cont << std::endl;
                        std::cout << F.evaluate(m) << std::endl;
                        std::cout << "F = " << F << std::endl;
                        std::cout << "m = " << m << std::endl;
                        std::cout << "s = " << s << std::endl;
                        std::cout << "alpha = " << alpha << std::endl;
                        std::cout << "E(F) = " << E << std::endl;
                    }
                }
            }
        }
    }

    if (best_alpha == 0.0)
    {
        std::cout << "j0 = " << j0_best << " and j1 = " << j1_best << ", cont = " << best << std::endl;
        std::vector<VeryLong> c;
        c.resize(3);
        c[0] = VeryLong(j0_best) * m;
        c[1] = -one * (VeryLong(j0_best) + m * VeryLong(j1_best));
        c[2] = j1_best;
        Polynomial<VeryLong> F = min_poly + Polynomial<VeryLong>(c);
        std::cout << F.evaluate(m) << std::endl;
        std::cout << "F = " << F << std::endl;
        double alpha = PolynomialOptimizer::alpha_F(F, 2000, 200);
        if (alpha < best_alpha)
        {
            best_alpha = alpha;
            best_F = F;
        }
        std::cout << "alpha = " << alpha << std::endl;
    }

    s = PolynomialOptimizer::minimize_I_over_s(Polynomial<VeryLong>::convert_to_double<double>(best_F), 1000000.0);
    std::cout << "Best F = " << best_F << std::endl;
    std::cout << "m = " << m << std::endl;
    std::cout << "alpha = " << best_alpha << std::endl;
    std::cout << "New s = " << s << std::endl;

    return best_F;
}

Polynomial<VeryLong> base_m_polynomial(const VeryLong& N, const VeryLong& a, const VeryLong& b,
                                       const VeryLong& c_d, const VeryLong& c_d_1, int d)
{
    const VeryLong zero(0L);
    std::vector<VeryLong> r(d + 1);
    std::vector<VeryLong> c(d + 1);

    c[d] = c_d;
    r[d] = N;

    c[d-1] = c_d_1;
    r[d-1] = (r[d] - c[d] * pow<VeryLong, long int>(b, d)) / a;

    for (int i = d - 2; i >= 0; --i)
    {
        r[i] = (r[i+1] - c[i+1] * pow<VeryLong, long int>(b, i+1)) / a;
        VeryLong b_power = pow<VeryLong, long int>(b, i);
        VeryLong b_power_inv = b_power.inverse(a);
        VeryLong c_i_mod_a = (r[i] * b_power_inv) % a;
        VeryLong c_i = r[i] / pow<VeryLong, long int>(b, i);
        c[i] = c_i - c_i % a + c_i_mod_a;
        VeryLong delta = c[i] - c_i;
        //std::cout << "base_m_polynomial() : delta[" << i << "] = " << delta << std::endl;
        while (delta < zero)
        {
            delta += a;
            c[i] += a;
        }

        while (delta >= a)
        {
            delta += a;
            c[i] += a;
        }


    }
    Polynomial<VeryLong> poly(c);
    return poly;
}

#if 0
Polynomial<VeryLong> base_m_polynomial(const VeryLong& N, const VeryLong& a, const VeryLong& b,
                                       const VeryLong& c_d, int d)
{
    const VeryLong zero(0L);
    std::vector<VeryLong> r(d + 1);
    std::vector<VeryLong> c(d + 1);

    c[d] = c_d;
    r[d] = N;

    for (int i = d - 1; i >= 0; --i)
    {
        r[i] = (r[i+1] - c[i+1] * pow<VeryLong, long int>(b, i+1)) / a;
        VeryLong b_power = pow<VeryLong, long int>(b, i);
        VeryLong b_power_inv = b_power.inverse(a);
        VeryLong c_i_mod_a = (r[i] * b_power_inv) % a;
        VeryLong c_i = r[i] / pow<VeryLong, long int>(b, i);
        c[i] = c_i - c_i % a + c_i_mod_a;
        VeryLong delta = c[i] - c_i;
        //std::cout << "base_m_polynomial() : delta[" << i << "] = " << delta << std::endl;
        while (delta < zero)
        {
            delta += a;
            c[i] += a;
        }

        while (delta >= a)
        {
            delta += a;
            c[i] += a;
        }


    }
    Polynomial<VeryLong> poly(c);
    return poly;
}

void mpz_set_ll(mpz_t& vl_, long long l)
{
    int negative = 0;
    if (l < 0)
    {
        l = -l;
        negative = 1;
    }
    unsigned long int q = (unsigned long int)(l >> 32);
    unsigned long int r = (unsigned long int)(l & 0x00000000ffffffff);
    mpz_set_ui(vl_, q);
    mpz_mul_2exp(vl_, vl_, 32);
    mpz_add_ui(vl_, vl_, r);
    if (negative) mpz_neg(vl_, vl_);
}

void mpz_set_ull(mpz_t& vl_, unsigned long long l)
{
    unsigned long int q = (unsigned long int)(l >> 32);
    unsigned long int r = (unsigned long int)(l & 0x00000000ffffffff);
    mpz_set_ui(vl_, q);
    mpz_mul_2exp(vl_, vl_, 32);
    mpz_add_ui(vl_, vl_, r);
}
#endif
};
namespace
{
#if 0
Polynomial<VeryLong> base_m_polynomial(const VeryLong& N, const long int degree, VeryLong* m)
{
    // create a Polynomial of degree d from N, by
    // expanding N base m where m ~ N^(1/d) or should it be N^(1/(d+1)) ?
    // if m ~ N^(1/d) then we get a monic polynomial,
    // but m ~ N^(1/(d+1)) gives a polynomial with more evenly matched coefficients
    *m = N.nth_root(degree+1);
    VeryLong mp = exp(*m, VeryLong(degree));
    VeryLong r = N;

    vector<VeryLong> coefficients;
    coefficients.resize(degree + 1);
    int d = degree;
    while (d >= 0)
    {
        VeryLong c = r / mp;
        coefficients[d] = c;
        std::cout << d << " : " << c << std::endl;
        r -= c * mp;
        mp /= *m;
        d--;
    }

    Polynomial<VeryLong> poly(coefficients);
    std::cout << "m = " << *m << std::endl;
    std::cout << poly.evaluate(*m) << std::endl;
    return poly;
}
#endif

Polynomial<VeryLong> base_m_polynomial_1(const VeryLong& N, const long int degree, const VeryLong& m)
{
    VeryLong mp = exp(m, VeryLong(degree));
    VeryLong r = N;

    vector<VeryLong> coefficients;
    coefficients.resize(degree + 1);
    int d = degree;
    while (d >= 0)
    {
        VeryLong c = r / mp;
        coefficients[d] = c;
        r -= c * mp;
        mp /= m;
        d--;
    }

    Polynomial<VeryLong> poly(coefficients);
    return poly;
}
#if 0
// Polynomial selection for NFS

double rate_polynomial(const Polynomial<VeryLong>& f1,
                       const Polynomial<VeryLong>& f2,
                       double alpha,
                       long int B1, long int B2, long int d)
{
    const int K = 1000;

    Polynomial<double> F1 = Polynomial<VeryLong>::convert_to_double<double>(f1);
    Polynomial<double> F2 = Polynomial<VeryLong>::convert_to_double<double>(f2);
    double rating = 0.0;
    double logB1 = log((double)B1);
    double logB2 = log((double)B2);
    for (int i = 1; i <= K; i++)
    {
        double theta = M_PI/K*(i - 0.5);
        double tmp = F1.evaluate_homogeneous(cos(theta), sin(theta));
        if (tmp < 0) tmp = -tmp;
        double u1 = (log(tmp) + alpha) / logB1 + d;
        tmp = F2.evaluate_homogeneous(cos(theta), sin(theta));
        if (tmp < 0) tmp = -tmp;
        double u2 = log(tmp) / logB2 + 1.0;
        rating += dickman_rho(u1) * dickman_rho(u2);
    }
    return rating / K;
}

double rate_skewed_polynomial(const Polynomial<VeryLong>& f1,
                              const Polynomial<VeryLong>& f2,
                              double alpha, double alpha2,
                              long int B1, long int B2, long int d, long int s)
{
    const int K = 1000;

    Polynomial<double> F1 = Polynomial<VeryLong>::convert_to_double<double>(f1);
    Polynomial<double> F2 = Polynomial<VeryLong>::convert_to_double<double>(f2);
    if (s == 0L)
    {
        std::cerr << "Bad s passed to rate_skewed_polynomial, s set to 1000" << std::endl;
        s = 1000L;
    }
    double s1 = sqrt((double)s);
    double s2 = 1 / s1;
    double rating = 0.0;
    double logB1 = log((double)B1);
    double logB2 = log((double)B2);
    for (int i = 1; i <= K; i++)
    {
        double theta = M_PI/K*(i - 0.5);
        double tmp = fabs(F1.evaluate_homogeneous(s1 * cos(theta), s2 * sin(theta)));
        double u1 = (log(tmp) + alpha) / logB1;
        tmp = fabs(F2.evaluate_homogeneous(s1 * cos(theta), s2 * sin(theta)));
        double u2 = (log(tmp) + alpha2) / logB2;
        rating += dickman_rho(u1) * dickman_rho(u2);
    }
    return rating;
}
#endif
}

void display_mu(const std::vector<long int>& mu)
{
    std::cout << "mu = (";
    for (size_t i = 0; i < mu.size() - 1; ++i)
    {
        std::cout << mu[i] << ",";
    }
    std::cout << mu[mu.size() - 1] << ")" << std::endl;
}


VeryLong calculate_m_mu(const std::vector<std::vector<VeryLong> >& m,
                        const std::vector<long int>& mu)
{
    //
    //       ---
    // m   = \   m
    //  mu   /    i,mu
    //       ---      i
    //        i
    //
    VeryLong m_mu(0L);
    for (size_t i = 0; i < m.size(); ++i)
    {
        m_mu += m[i][mu[i]];
    }

    return m_mu;
}

bool check_c_d(const VeryLong& c_d, const VeryLong& N, const VeryLong& a, const VeryLong& b, long int d)
{
    const VeryLong zero(0L);
    // check that
    //     d
    // c  b = N mod a
    //  d
    VeryLong check(N);
    check -= c_d * pow<VeryLong, long int>(b, d);
    check %= a;
    if (check != zero)
    {
        std::cout << "Problem: c   = " << c_d << " does not work" << std::endl;
        std::cout << "          d" << std::endl;
    }
    return (check == zero);
}

bool check_c_d_1(const VeryLong& c_d_1, const VeryLong& N, const VeryLong& c_d,
                 const VeryLong& m_mu, const VeryLong& a, long int d)
{
    const VeryLong zero(0L);
    // check that
    //          d-1            d
    // c       m    = (N - c  m  ) / a (mod a)
    //  d-1,mu  mu          d  mu

    VeryLong check(N);
    check -= c_d * pow<VeryLong, long int>(m_mu, d);
    if (check % a != zero)
    {
        std::cout << "Problem: check % a != 0 : " << check % a << std::endl;
    }
    check /= a;
    check -= c_d_1 * pow<VeryLong, long int>(m_mu, d-1);
    check %= a;
    if (check != zero)
    {
        std::cout << "Problem: c       = " << c_d_1 << " does not work" << std::endl;
        std::cout << "          d-1,mu" << std::endl;
    }
    return (check == zero);
}


class PolynomialPairCalculator
{
public:
    PolynomialPairCalculator(const VeryLong& N, const VeryLong& c_d, std::fstream& output_file, bool debug = false)
        : N_(N), c_d_(c_d), output_file_(output_file), debug_(debug)
    {}
    ~PolynomialPairCalculator() {}

    bool generate(long int degree);

private:
    /*
     * Each flist_item records a sum, and a mu index
     * The mu_index encodes the permutation mu as digits in base d
     */
    struct flist_item
    {
        flist_item(double sum, size_t mu_index)
            : sum_(sum), mu_index_(mu_index)
        {
        }
        flist_item()
            : sum_(0.0), mu_index_(0L)
        {
        }

        double sum_;
        size_t mu_index_;
        bool operator<(const flist_item& fi) const
        {
            return (sum_ < fi.sum_);
        }

        void display() const
        {
            std::cout << "flist_item : sum_ = " << sum_ << ", mu_index_ = " << mu_index_ << std::endl;
        }

        void display(long int degree, size_t mu_length) const
        {
            std::vector<long int> mu;
            get_mu(degree, mu, mu_length);
            std::cout << "flist_item : sum_ = " << sum_ << ", mu_index_ = " << mu_index_ << " : (";
            for (size_t i = 0; i < mu.size(); ++i)
            {
                std::cout << mu[i];
                if (i < mu.size() - 1)
                    std::cout << ",";
                else
                    std::cout << ")" << std::endl;
            }
        }

        // Compute a permutation mu from the index in mu_index_
        void get_mu(long int degree, std::vector<long int>& mu, size_t mu_length) const
        {
            const long int d(degree);
            size_t i = 0;
            long int digit = 0L;
            size_t index = mu_index_;
            while (index > 0L)
            {
                if (i >= mu_length)
                    throw "Index out of range in get_mu";

                digit = index % d;
                mu.push_back(digit);

                index -= digit;
                index /= d;
                ++i;
            }

            for (; i < mu_length; ++i)
            {
                mu.push_back(0L);
            }
        }
    };

    template <size_t S, size_t M>
    class flist_hash_table_
    {
    public:
        void add(const flist_item& item)
        {
            size_t slot = (item.sum_ + 0.5) * S;
            if (slot >= S)
                slot = S - 1;
            flist_list& fl = hash_table_[slot];
            if (fl.item_count_ < M)
            {
                fl.item_[fl.item_count_] = item;
                ++fl.item_count_;
            }
        }
        void check_for_match(long int degree,
                             const flist_item& item,
                             double epsilon,
                             size_t flist1_mu_length,
                             size_t flist2_mu_length,
                             std::vector<std::vector<long int> >& good_mu,
                             bool debug = false)
        {
            const long int d(degree);
            // we are looking for it such that
            // -(it->sum_) - epsilon <= item.sum_ <= -(it->sum_) + epsilon)
            // <=>
            // -epsilon <= it->sum_ + item.sum_ <= epsilon)
            // i.e. it->sum_ + item.sum_ is very close to zero
            //
            // Given item, we want to find the slots in hash_table which could contain such *it
            // i.e. -item.sum_ is close to it->sum_
            double lb = -item.sum_ - epsilon;
            double ub = -item.sum_ + epsilon;
            size_t slot = (-item.sum_ + 0.5) * S;
            if (slot >= S)
                slot = S - 1;
            if (debug)
            {
                std::cout << "Looking for match in slot " << slot << " for item : ";
                item.display(degree, flist1_mu_length);
            }
            flist_list& fl = hash_table_[slot];
            if (debug)
            {
                if (fl.item_count_ > 0)
                {
                    std::cout << "Slot " << slot << std::endl;
                    fl.display();
                }
            }
            for (flist_item* it = fl.item_;
                    it != fl.item_ + fl.item_count_;
                    ++it)
            {
                if (it->sum_ >= lb && it->sum_ <= ub)
                {
                    std::vector<long int> mu;
                    mu.reserve(flist1_mu_length + flist2_mu_length);
                    item.get_mu(d, mu, flist1_mu_length);
                    it->get_mu(d, mu, flist2_mu_length);
                    good_mu.push_back(mu);
                    if (debug)
                    {
                        std::cout << std::setprecision(16) << "Match found: epsilon = " << epsilon << ", it->sum_ = " << it->sum_ << ", item.sum_ = " << item.sum_ << ", ";
                        display_mu(mu);
                    }
                }
            }
        }
        void display() const
        {
            for (size_t i = 0; i < S; ++i)
            {
                const flist_list& fll = hash_table_[i];
                if (fll.item_count_ > 0)
                {
                    std::cout << "Slot " << i << " : " << std::endl;
                }
                fll.display();
            }
        }
        void clear()
        {
            for (size_t i = 0; i < S; ++i)
            {
                flist_list& fll = hash_table_[i];
                fll.item_count_ = 0;
            }
        }
    private:
        struct flist_list
        {
            size_t item_count_;
            flist_item item_[M];
            flist_list() : item_count_(0) {}
            void display() const
            {
                for (size_t j = 0; j < item_count_; ++j)
                {
                    item_[j].display();
                }
            }
        };
        flist_list& find(double sum)
        {
            size_t slot = (sum + 0.5) * S;
            if (slot >= S)
                slot = S - 1;
            return hash_table_[slot];
        }
        flist_list hash_table_[S];
    };

    enum { number_of_slots = 512L, max_items = 8L };

    typedef flist_hash_table_<number_of_slots, max_items> flist_hash_table;

    typedef std::vector<flist_item> flist_collection;

    static double modZ(double d)
    {
        double f = ::fmod(d, 1.0);
        if (::fabs(f) < 0.5)
        {
            return f;
        }
        if (f >= 0.5)
        {
            return f - 1.0;
        }
        return f + 1.0;
    }

    static VeryLong calculate_c_d_1(const VeryLong& N,
                                    const VeryLong& c_d,
                                    long int degree,
                                    long int primes_to_combine,
                                    const VeryLong& a,
                                    const VeryLong& N_inv_c_d,
                                    const VeryLong& m_mu_0,
                                    const VeryLong& c_d_1_0,
                                    const std::vector<std::vector<VeryLong> >& m,
                                    const std::vector<long int>& mu);

    struct XYZ
    {
        XYZ(const PolynomialPairCalculator& ppc, const Combinations& combination, long int primes_to_combine)
            : ppc_(ppc), combination_(combination), primes_to_combine_(primes_to_combine), p_(primes_to_combine_), a_(1L), x_(primes_to_combine_),
              m_(primes_to_combine_), f_(primes_to_combine_) {}
        XYZ(const PolynomialPairCalculator::XYZ& xyz, const VeryLong& q, const VeryLong& s);
        XYZ(const PolynomialPairCalculator::XYZ& xyz, const VeryLong& m_adjustment);
        void generate(std::vector<Kleinjung_poly_info>& top_polys);
        void generate(long int iterations, const VeryLong& q, const VeryLong& s, std::vector<Kleinjung_poly_info>& top_polys);
        bool make_a();
        void make_x();
        void make_f();
        void make_flists();
        void find_good_mu();
        void process_good_mu(std::vector<Kleinjung_poly_info>& top_polys);

        const PolynomialPairCalculator& ppc_;
        const Combinations& combination_;
        long int primes_to_combine_;
        std::vector<long int> p_;
        VeryLong a_;
        VeryLong m0_;
        std::vector<std::vector<double> > x_;
        std::vector<std::vector<VeryLong> > m_;
        VeryLong m_mu_0_;
        VeryLong c_d_1_0_;
        std::vector<std::vector<double> > f_;
        double f0_;
        VeryLong N_inv_c_d_;
        std::vector<std::vector<long int> > good_mu_;
        flist_collection flist1_;
        size_t flist1_mu_length_;
        flist_hash_table flist2_;
        size_t flist2_mu_length_;
    };

    VeryLong N_;
    long int d_;
    std::vector<long int> d_powers_;
    Polynomial<VeryLong> poly_;
    VeryLong c_d_;
    VeryLong c_d_1_max_;
    VeryLong c_d_2_max_;
    double minus_c_d_d_d_;
    std::vector<Kleinjung_poly_info> top_polys_;
    mutable VeryLong prev_top_c_d_1_;
    VeryLong m_;
    std::vector<long int> primes_;
    std::vector<long int> extra_primes_;
    std::vector<std::vector<LongModular> > roots_mod_p_;
    std::vector<LongModular > extra_roots_mod_p_;
    std::fstream& output_file_;
    mutable int count_;
    mutable bool debug_;
};

VeryLong PolynomialPairCalculator::calculate_c_d_1(const VeryLong& N,
        const VeryLong& c_d,
        long int degree,
        long int primes_to_combine,
        const VeryLong& a,
        const VeryLong& N_inv_c_d,
        const VeryLong& m_mu_0,
        const VeryLong& c_d_1_0,
        const std::vector<std::vector<VeryLong> >& m,
        const std::vector<long int>& mu)
{
    //std::cout << "DEBUG : ";
    //display_mu(mu);
    const VeryLong zero(0L);
    //
    //           ---
    // c       = \   e
    //  d-1,mu   /    i,mu
    //           ---      i
    //            i
    //
    VeryLong c_d_1(zero);
    //const int d(5L);
    const int d(degree);
    // i = 0
    {
        if (mu[0] == 0)
        {
            VeryLong m_mu(m_mu_0);
            VeryLong e_i_j(N);
            e_i_j -= c_d * pow<VeryLong, long int>(m_mu, d);
            e_i_j *= m_mu;
            c_d_1 += e_i_j;
        }
        else
        {
            VeryLong m_mu = m_mu_0 + m[0][mu[0]] - m[0][0];
            VeryLong e_i_j(N);
            e_i_j -= c_d * pow<VeryLong, long int>(m_mu, d);
            e_i_j *= m_mu;
            c_d_1 += e_i_j;
        }
    }

    // i > 0
    for (long int i = 1; i < primes_to_combine; ++i)
    {
        if (mu[i] > 0)
        {
            VeryLong m_mu = m_mu_0 + m[i][mu[i]] - m[i][0];
            VeryLong e_i_j(N);
            e_i_j -= c_d * pow<VeryLong, long int>(m_mu, d);
            e_i_j *= m_mu;
            e_i_j -= c_d_1_0;
            c_d_1 += e_i_j;
        }
    }
    if (c_d_1 % a != zero)
    {
        std::cout << "Problem: c_d_1 % a = " << c_d_1 % a << std::endl;
    }
    c_d_1 /= a;
    c_d_1 *= N_inv_c_d;
    c_d_1 %= a;

    return c_d_1;
}

bool PolynomialPairCalculator::XYZ::make_a()
{
    const VeryLong zero(0L);
    for (int i = 0; i < primes_to_combine_; ++i)
    {
        p_[i] = ppc_.primes_[combination_(i)];
        a_ *= p_[i];
    }

    if (ppc_.debug_)
    {
        std::cout << "a        = " << a_ << std::endl;
        std::cout << "c        = " << ppc_.c_d_1_max_ << std::endl;
        std::cout << " d-1,max" << std::endl;
    }

    if (a_ > ppc_.c_d_1_max_)
    {
        return false;
    }

    if (ppc_.debug_)
    {
        std::cout << "a        = " << a_ << std::endl;
    }

    //
    // m  is the smallest integer bigger than m_ and divisible by a
    //  0
    //
    m0_ = ppc_.m_;
    if (m0_ % a_ != zero)
    {
        m0_ -= m0_ % a_;
        m0_ += a_;
    }
    if (ppc_.debug_)
    {
        std::cout << "m        = " << m0_ << std::endl;
        std::cout << " 0" << std::endl;
    }

    return true;
}

void PolynomialPairCalculator::XYZ::make_x()
{
    const VeryLong zero(0L);
    //
    // find x    satisfying
    //       i,j
    //
    // x    = r    mod p
    //  i,j    i,j      i
    //
    //                                    d
    // where r    are the roots of N = c X  mod p
    //        i,j                       d        i
    //
    // and a / p  | x
    //          i    i,j
    //

    // i = 0
    {
        VeryLong a_p_i = a_ / p_[0];
        VeryLong a_p_i_inv = a_p_i.inverse(p_[0]);
        for (size_t j = 0; j < ppc_.roots_mod_p_[combination_(0)].size(); ++j)
        {
            const LongModular& r = ppc_.roots_mod_p_[combination_(0)][j];
            m_[0].push_back(r.get_long());
            VeryLong& x_i_j = m_[0][m_[0].size() - 1];
            x_i_j *= a_p_i_inv;
            x_i_j *= a_p_i;
            x_i_j %= a_;
            x_[0].push_back(x_i_j.get_double());
            x_i_j += m0_;
#ifdef CHECK
            {
                // check
                VeryLong check = x_i_j - VeryLong(r.get_long());
                if (check % p_[0] != zero)
                {
                    std::cout << "Problem : x   != r    mod p" << std::endl;
                    std::cout << "           0,j    i,j      0" << std::endl;
                }
                VeryLong check1 = ppc_.c_d_ * pow<VeryLong, long int>(x_i_j, ppc_.d_) - ppc_.N_;
                VeryLong check2 = ppc_.poly_.evaluate(x_i_j);
                if (check1 % p_[0] != zero)
                {
                    std::cout << "Problem : x    doesn't work" << std::endl;
                    std::cout << "           0," << j << std::endl;
                    std::cout << "check1 = " << check1 << std::endl;
                    std::cout << "check2 = " << check2 << std::endl;
                    std::cout << "check2 % p[0] = " << check2 % p_[0] << std::endl;
                }
            }
#endif
        }
    }

    for (int i = 1; i < primes_to_combine_; ++i)
    {
        VeryLong a_p_i = a_ / p_[i];
        VeryLong a_p_i_inv = a_p_i.inverse(p_[i]);
        for (size_t j = 0; j < ppc_.roots_mod_p_[combination_(i)].size(); ++j)
        {
            const LongModular& r = ppc_.roots_mod_p_[combination_(i)][j];
            m_[i].push_back(r.get_long());
            VeryLong& x_i_j = m_[i][m_[i].size() - 1];
            x_i_j *= a_p_i_inv;
            x_i_j *= a_p_i;
            x_i_j %= a_;
            x_[i].push_back(x_i_j.get_double());
#ifdef CHECK
            {
                // check
                VeryLong check = x_i_j - VeryLong(r.get_long());
                if (check % p_[i] != zero)
                {
                    std::cout << "Problem : x  != r    mod p" << std::endl;
                    std::cout << "           i,j   i,j      i" << std::endl;
                }
                VeryLong check1 = ppc_.c_d_ * pow<VeryLong, long int>(x_i_j, ppc_.d_) - ppc_.N_;
                if (check1 % p_[i] != zero)
                {
                    std::cout << "Problem : x    doesn't work" << std::endl;
                    std::cout << "          " << i << "," << j << std::endl;
                }
            }
#endif
        }
    }

    //
    // Now x_[i] has ppc_.roots_mod_p_[combination_(i)].size() elements in it
    // Similarly m_[i]
    //

    for (size_t i = 0; i < m_.size(); ++i)
    {
        m_mu_0_ += m_[i][0];
    }
#ifdef CHECK
    {
        // check
        VeryLong check(ppc_.N_);
        check -= ppc_.c_d_ * pow<VeryLong, long int>(m_mu_0_, ppc_.d_);
        if (check % a_ != zero)
        {
            std::cout << "Problem: m_mu_0 is not working!" << std::endl;
            std::cout << "m_mu_0 = " << m_mu_0_ << std::endl;
            std::cout << "m_mu_0 - m0 = " << m_mu_0_ - m0_ << std::endl;
            std::cout << "check = " << check << std::endl;
            std::cout << "check % a = " << check % a_ << std::endl;
        }
    }
#endif
    c_d_1_0_ = ppc_.N_;
    c_d_1_0_ -= ppc_.c_d_ * pow<VeryLong, long int>(m_mu_0_, ppc_.d_);
    c_d_1_0_ *= m_mu_0_;
}

void PolynomialPairCalculator::XYZ::make_f()
{
    //
    // calculate e
    //            i,j
    //
    // the first index, i, corresponds to the primes p  dividing a
    //                                                i
    // the second index, j, corresponds to the j-th root mod p
    //                                                        i
    //
    // e    = c                mod a
    //  0,j    d-1,(j,0,...,0)
    //
    // e    = 0  i > 0
    //  i,0
    //
    // e    = c                         - c              mod a
    //  i,j    d-1,(0,...,0,j,0,....,0)    d-1,(0,...,0)             i > 0, j > 0
    //                      ^
    //                      i-th set to j
    //

    double a_d(a_.get_double());
    double a_d_2(a_d);
    a_d_2 *= a_d;
    //
    // f  is dependent on m
    //  0                  0
    //
    VeryLong f0_vl(ppc_.N_);
    f0_vl -= ppc_.c_d_ * pow<VeryLong, long int>(m0_, ppc_.d_);
    f0_ = f0_vl.get_double();
    f0_ /= a_d_2;
    f0_ /= (pow<VeryLong, long int>(m0_, ppc_.d_-1)).get_double();

    if (ppc_.debug_)
    {
        std::cout << "f        = " << f0_ << std::endl;
        std::cout << " 0" << std::endl;
    }

    N_inv_c_d_ = ppc_.N_.inverse(a_) * ppc_.c_d_;

    // i = 0
    {
        f_[0].reserve(ppc_.d_);
        // j = 0
        {
            double f_i_j(ppc_.minus_c_d_d_d_);
            f_i_j *= x_[0][0];
            VeryLong e_i_j(ppc_.N_);

            //
            // dependency on m    here
            //                mu
            //                  0
            //
            e_i_j -= ppc_.c_d_ * pow<VeryLong, long int>(m_mu_0_, ppc_.d_);
            e_i_j *= m_mu_0_;
            e_i_j *= N_inv_c_d_;
            e_i_j /= a_;
            e_i_j %= a_;
            e_i_j *= a_;
            f_i_j -= e_i_j.get_double();
            f_i_j /= a_d_2;
            f_[0].push_back(f_i_j);
        }

        //for (long int j = 1; j < ppc_.d_; ++j)
        for (long int j = 1; j < (long int)ppc_.roots_mod_p_[combination_(0)].size(); ++j)
        {
            double f_i_j(ppc_.minus_c_d_d_d_);
            f_i_j *= x_[0][j];
            //
            // dependency on m    here
            //                mu
            //                  0
            //
            VeryLong m_mu = m_mu_0_ + m_[0][j] - m_[0][0];
            VeryLong e_i_j(ppc_.N_);
            e_i_j -= ppc_.c_d_ * pow<VeryLong, long int>(m_mu, ppc_.d_);
            e_i_j *= m_mu;
            e_i_j *= N_inv_c_d_;
            e_i_j /= a_;
            e_i_j %= a_;
            e_i_j *= a_;
            f_i_j -= e_i_j.get_double();
            f_i_j /= a_d_2;
            f_[0].push_back(f_i_j);
        }
    }
    for (int i = 1; i < primes_to_combine_; ++i)
    {
        f_[i].reserve(ppc_.d_);
        // j = 0
        double f_i_j(ppc_.minus_c_d_d_d_);
        f_i_j *= x_[i][0];
        f_i_j /= a_d_2;
        f_[i].push_back(f_i_j);

        //for (long int j = 1; j < ppc_.d_; ++j)
        for (long int j = 1; j < (long int)ppc_.roots_mod_p_[combination_(i)].size(); ++j)
        {
            double f_i_j(ppc_.minus_c_d_d_d_);
            f_i_j *= x_[i][j];
            VeryLong m_mu = m_mu_0_ + m_[i][j] - m_[i][0];
            VeryLong e_i_j(ppc_.N_);
            e_i_j -= ppc_.c_d_ * pow<VeryLong, long int>(m_mu, ppc_.d_);
            e_i_j *= m_mu;
            e_i_j -= c_d_1_0_;
            e_i_j *= N_inv_c_d_;
            e_i_j /= a_;
            e_i_j %= a_;
            e_i_j *= a_;
            f_i_j -= e_i_j.get_double();
            f_i_j /= a_d_2;
            f_[i].push_back(f_i_j);
        }
    }
}

void PolynomialPairCalculator::XYZ::make_flists()
{
    // Make two lists
    //
    //      u
    //     ---
    // f + \  f
    //  0  /   i,mu
    //     ---     i
    //     i=1
    //
    //      l
    //     ---
    //   - \  f
    //     /   i,mu
    //     ---     i
    //     i=u+1
    //
    // where u = [ l / 2 ]

    long int u = (primes_to_combine_ / 2) + 1;
    if ((primes_to_combine_ & 1) == 0)
    {
        --u;
    }
    flist1_mu_length_ = u;
    flist2_mu_length_ = primes_to_combine_ - u;

    long int d_u = std::pow(ppc_.d_, u);
    long int d_l_u = std::pow(ppc_.d_, primes_to_combine_ - u);
    flist1_.clear();
    flist1_.reserve(ppc_.d_powers_[u]);
    flist2_.clear();
    for (size_t n = 0; n < d_u; ++n)
    {
        double sum = f0_;
        for (size_t i = 0; i < u; ++i)
        {
            sum += f_[i][(n / ppc_.d_powers_[i]) % ppc_.d_];
        }
        flist1_.push_back(flist_item(modZ(sum), n));
    }
    for (size_t n = 0; n < d_l_u; ++n)
    {
        double sum = 0.0;
        for (size_t i = 0; i < primes_to_combine_ - u; ++i)
        {
            sum += f_[u + i][(n / ppc_.d_powers_[i]) % ppc_.d_];
        }
        flist2_.add(flist_item(modZ(sum), n));
    }
    if (ppc_.debug_)
    {
        std::cout << "flist1 : " << std::endl;
        for (auto& i1: flist1_)
        {
            i1.display(ppc_.d_, u);
        }
        std::cout << "flist2 : " << std::endl;
        flist2_.display();
    }
}

void PolynomialPairCalculator::XYZ::find_good_mu()
{
    const double fudge(10.0);
    double epsilon = ppc_.c_d_2_max_.get_double() / m0_.get_double();
    epsilon /= fudge;
    if (ppc_.debug_)
    {
        std::cout << "epsilon  = " << epsilon << std::endl;
    }

    for (auto& i1: flist1_)
    {
        flist2_.check_for_match(ppc_.d_, i1, epsilon, flist1_mu_length_, flist2_mu_length_, good_mu_, ppc_.debug_);
    }
}

void PolynomialPairCalculator::XYZ::process_good_mu(std::vector<Kleinjung_poly_info>& top_polys)
{
    for (size_t i = 0; i < good_mu_.size(); ++i)
    {
        VeryLong c_d_1 = PolynomialPairCalculator::calculate_c_d_1(ppc_.N_, ppc_.c_d_, ppc_.d_, primes_to_combine_, a_, N_inv_c_d_, m_mu_0_, c_d_1_0_, m_, good_mu_[i]);
        VeryLong b = calculate_m_mu(m_, good_mu_[i]);
        if (ppc_.debug_)
        {
            if (b < ppc_.m_)
            {
                std::cout << "b = " << b << " is less than m~ = " << ppc_.m_ << std::endl;
            }
            else
            {
                std::cout << "b = " << b << " is not less than m~ = " << ppc_.m_ << std::endl;
            }
        }
        if (ppc_.debug_)
        {
            check_c_d(ppc_.c_d_, ppc_.N_, a_, b, ppc_.d_);
            check_c_d_1(c_d_1, ppc_.N_, ppc_.c_d_, b, a_, ppc_.d_);
        }
        Polynomial<VeryLong> fm = base_m_polynomial(ppc_.N_, a_, b, ppc_.c_d_, c_d_1, ppc_.d_);
        if (ppc_.debug_)
        {
            display_mu(good_mu_[i]);
            std::cout << "fm = " << fm << std::endl;
        }
        if (ppc_.debug_)
        {
            VeryLong fm_b_a_mod_N = fm.evaluate_homogeneous(b, a_) % ppc_.N_;
            std::cout << "fm(b,a) mod N = " << fm_b_a_mod_N << std::endl;
        }
        VeryLong c_d_2 = fm.coefficient(ppc_.d_ - 2);
        if (ppc_.debug_)
        {
            std::cout << "c     = " << c_d_2 << std::endl;
            std::cout << " d-2 " << std::endl;
        }
        if (ppc_.debug_)
        {
            double f1 = c_d_2.get_double() / m0_.get_double();
            std::cout << "c    / m  = " << f1 << std::endl;
            std::cout << " d-2    0" << std::endl;
        }
        VeryLong k = c_d_2 / m0_;
        VeryLong r = c_d_2 % m0_;
        if (r > m0_ / 2L)
        {
            k += 1L;
            r -= m0_;
        }
        if (ppc_.debug_)
        {
            std::cout << "c    / m  = " << k << std::endl;
            std::cout << " d-2    0" << std::endl;
            std::cout << "c    % m  = " << r << std::endl;
            std::cout << " d-2    0" << std::endl;
        }
        fm.set_coefficient(ppc_.d_ - 1, c_d_1 + k * a_);
        fm.set_coefficient(ppc_.d_ - 2, c_d_2 - k * b);
        if (ppc_.debug_)
        {
            std::cout << "adjusted fm = " << fm << std::endl;
            std::cout << "a           = " << a_ << std::endl;
            std::cout << "b           = " << b << std::endl;
        }
        if (ppc_.debug_)
        {
            VeryLong fm_b_a_mod_N = fm.evaluate_homogeneous(b, a_) % ppc_.N_;
            std::cout << "fm(b,a) mod N = " << fm_b_a_mod_N << std::endl;
        }

#if 0
        {
            VeryLongModular::set_default_modulus(ppc_.N_);
            VeryLongModular tmp1 = VeryLongModular(b) / VeryLongModular(a_);
            VeryLong m = tmp1.get_very_long();
            long int s = (long int)PolynomialOptimizer::minimize_I_over_s(Polynomial<VeryLong>::convert_to_double<double>(fm), 1000.0);
            //std::cout << "1. s = " << s << std::endl;
            VeryLong s_vl = s;
            //MPFloat::set_precision(100);
            //MPFloat I_F_S = 0.0;
            double I_F_S = 0.0;
            VeryLong new_b;
            VeryLong new_m;
            //Polynomial<VeryLong> new_poly = PolynomialOptimizer::minimize_I<MPFloat>(fm, a_, b, m, s_vl, I_F_S, new_b, new_m);
            Polynomial<VeryLong> new_poly = PolynomialOptimizer::minimize_I<double>(fm, a_, b, m, s_vl, I_F_S, new_b, new_m);
            //std::cout << "2. new_poly = " << new_poly << std::endl;
            //std::cout << "3. new_b = " << new_b << std::endl;
            //std::cout << "4. new_m = " << new_m << std::endl;
            if (ppc_.debug_)
            {
                VeryLong fm_b_a_mod_N = new_poly.evaluate_homogeneous(new_b, a_) % ppc_.N_;
                std::cout << "fm(b,a) mod N = " << fm_b_a_mod_N << std::endl;
            }
            b = new_b;
            fm = new_poly;
        }
#endif

        Kleinjung_poly_info pi(a_, b, fm);
        const size_t max_top_polys = 20;
        if (top_polys.size() < max_top_polys || pi < top_polys[max_top_polys - 1])
        {
            if (ppc_.debug_)
            {
                std::cout << "Adding (" << a_ << ", " << b << ", " << fm << ") to top_polys" << std::endl;
            }
            top_polys.push_back(pi);
            std::sort(top_polys.begin(), top_polys.end());
            VeryLong top_c_d_1 = top_polys[0].fm_.coefficient(0);
            if (top_c_d_1 != ppc_.prev_top_c_d_1_)
            {
                ppc_.prev_top_c_d_1_ = top_c_d_1;
                std::cout << "count = " << ppc_.count_ << std::endl;
                for (size_t i = 0; i < 3 && i < top_polys.size(); i++)
                {
                    std::cout << "fm[" << i << "] = " << top_polys[i].fm_ << std::endl;
                }
            }
        }
        if (top_polys.size() >= max_top_polys)
        {
            top_polys.erase(top_polys.end() - 1);
        }
        ++ppc_.count_;
        if (ppc_.count_ % 10000 == 0)
        {
            std::cout << "count_ = " << ppc_.count_ << std::endl;
        }
    }
}

void PolynomialPairCalculator::XYZ::generate(std::vector<Kleinjung_poly_info>& top_polys)
{
    if (!make_a())
    {
        return;
    }

    make_x();

    make_f();

    make_flists();

    find_good_mu();

    process_good_mu(top_polys);

#if 0
    for (int iteration = 1; iteration < 30; ++iteration)
    {
        PolynomialPairCalculator::XYZ newer_xyz(*this, iteration * 30000L);
        newer_xyz.good_mu_.clear();
        newer_xyz.find_good_mu();

        newer_xyz.process_good_mu(top_polys);
    }
#endif
}

PolynomialPairCalculator::XYZ::XYZ(const PolynomialPairCalculator::XYZ& xyz, const VeryLong& q, const VeryLong& s)
    : ppc_(xyz.ppc_), combination_(xyz.combination_), primes_to_combine_(xyz.primes_to_combine_), p_(xyz.p_),
      a_(xyz.a_), m0_(xyz.m0_), x_(xyz.x_), m_(xyz.m_), m_mu_0_(xyz.m_mu_0_), c_d_1_0_(xyz.c_d_1_0_),
      f_(xyz.f_), f0_(xyz.f0_), N_inv_c_d_(xyz.N_inv_c_d_), flist1_(xyz.flist1_), flist1_mu_length_(xyz.flist1_mu_length_),
      flist2_(xyz.flist2_), flist2_mu_length_(xyz.flist2_mu_length_)
{
    const VeryLong zero(0L);
    a_ *= q;

    //
    // m  is set here
    //  0
    //
    m0_ = ppc_.m_;
    if (m0_ % a_ != zero)
    {
        m0_ -= m0_ % a_;
    }

    if (ppc_.debug_)
    {
        std::cout << "m        = " << m0_ << std::endl;
        std::cout << " 0" << std::endl;
    }

    // i = 0
    {
        VeryLong q_inv = q.inverse(p_[0]);
        VeryLong a_inv = xyz.a_.inverse(q);
        VeryLong extra_term = xyz.a_;
        extra_term *= a_inv;
        extra_term *= s;
        extra_term %= a_;

        for (size_t j = 0; j < ppc_.roots_mod_p_[combination_(0)].size(); ++j)
        {
            VeryLong& x_i_j = m_[0][j];
            x_i_j -= xyz.m0_;
            x_i_j *= q;
            x_i_j *= q_inv;
            // add to m_[0][j] the term that comes from q, s
            // so that the sum of m_[i][mu[j]] over i always includes this term
            //
            //                   d
            // and so solves a  X  = N mod a
            //                d
            //
            // where a is the new product including q
            //
            x_i_j += extra_term;

            x_i_j %= a_;
            x_[0][j] = x_i_j.get_double();
            //
            // m  is added into m    here
            //  0                0 j
            //
            // but x    is not changed for each adjustment of m
            //      0 j                                        0
            //
            x_i_j += m0_;
        }
    }

    for (int i = 1; i < primes_to_combine_; ++i)
    {
        VeryLong q_inv = q.inverse(p_[i]);
        for (size_t j = 0; j < ppc_.roots_mod_p_[combination_(i)].size(); ++j)
        {
            VeryLong& x_i_j = m_[i][j];
            x_i_j *= q;
            x_i_j *= q_inv;
            x_i_j %= a_;
            x_[i][j] = x_i_j.get_double();
        }
    }

    m_mu_0_ = zero;
    for (size_t i = 0; i < m_.size(); ++i)
    {
        //
        // m    includes m  from m
        //  mu            0       0 0
        //    0
        //
        m_mu_0_ += m_[i][0];
    }

    //
    // c        is also dependent on m
    //  d-1 mu                         0
    //        0
    //
    c_d_1_0_ = ppc_.N_;
    c_d_1_0_ -= ppc_.c_d_ * pow<VeryLong, long int>(m_mu_0_, ppc_.d_);
    c_d_1_0_ *= m_mu_0_;

    f_.clear();
    f_.resize(primes_to_combine_);
    make_f();

    make_flists();
}

PolynomialPairCalculator::XYZ::XYZ(const PolynomialPairCalculator::XYZ& xyz, const VeryLong& m_adjustment)
    : ppc_(xyz.ppc_), combination_(xyz.combination_), primes_to_combine_(xyz.primes_to_combine_), p_(xyz.p_),
      a_(xyz.a_), m0_(xyz.m0_), x_(xyz.x_), m_(xyz.m_), m_mu_0_(xyz.m_mu_0_), c_d_1_0_(xyz.c_d_1_0_),
      f_(xyz.f_), f0_(xyz.f0_), N_inv_c_d_(xyz.N_inv_c_d_), flist1_(xyz.flist1_), flist1_mu_length_(xyz.flist1_mu_length_),
      flist2_(xyz.flist2_), flist2_mu_length_(xyz.flist2_mu_length_)
{
    const VeryLong zero(0L);
    //
    // m  is set here
    //  0
    //
    m0_ -= a_ * m_adjustment;

    if (ppc_.debug_)
    {
        std::cout << "m        = " << m0_ << std::endl;
        std::cout << " 0" << std::endl;
    }

    // i = 0
    {
        for (size_t j = 0; j < ppc_.roots_mod_p_[combination_(0)].size(); ++j)
        {
            VeryLong& x_i_j = m_[0][j];
            x_i_j -= xyz.m0_;
            x_i_j += m0_;
        }
    }

    //
    // don't need to change x    for i > 0
    //                       i j

    m_mu_0_ = zero;
    for (size_t i = 0; i < m_.size(); ++i)
    {
        //
        // m    includes m  from m
        //  mu            0       0 0
        //    0
        //
        m_mu_0_ += m_[i][0];
    }

    //
    // c        is also dependent on m
    //  d-1 mu                         0
    //        0
    //
    c_d_1_0_ = ppc_.N_;
    c_d_1_0_ -= ppc_.c_d_ * pow<VeryLong, long int>(m_mu_0_, ppc_.d_);
    c_d_1_0_ *= m_mu_0_;

    double a_d(a_.get_double());
    double a_d_2(a_d);
    a_d_2 *= a_d;

    double old_f0(f0_);

    VeryLong f0_vl(ppc_.N_);
    f0_vl -= ppc_.c_d_ * pow<VeryLong, long int>(m0_, ppc_.d_);
    f0_ = f0_vl.get_double();
    f0_ /= a_d_2;
    f0_ /= (pow<VeryLong, long int>(m0_, ppc_.d_-1)).get_double();

    //for (size_t j = 0; j < (size_t)ppc_.d_; ++j)
    for (size_t j = 0; j < f_[0].size(); ++j)
    {
        double old_f = f_[0][j];
        f_[0][j] -= ((VeryLong(ppc_.d_) * ppc_.c_d_ * m_adjustment) % a_).get_double() / a_.get_double();
        double new_f = f_[0][j];
        // f_[0][j] contributes to each combination mu whose index k satisfies k % d = j
        size_t index = j;
        while (index < flist1_.size())
        {
            flist_item& fli = flist1_[index];
            fli.sum_ -= old_f0;
            fli.sum_ += f0_;
            fli.sum_ -= old_f;
            fli.sum_ += new_f;
            index += ppc_.d_;
        }
    }
}

void PolynomialPairCalculator::XYZ::generate(long int iterations, const VeryLong& q, const VeryLong& s, std::vector<Kleinjung_poly_info>& top_polys)
{
    long int iteration = 0;

    PolynomialPairCalculator::XYZ new_xyz(*this, q, s);

    new_xyz.good_mu_.clear();
    new_xyz.find_good_mu();

    new_xyz.process_good_mu(top_polys);

    for (iteration = 1; iteration < iterations; ++iteration)
    {
        PolynomialPairCalculator::XYZ newer_xyz(new_xyz, iteration * 10000L);
        newer_xyz.good_mu_.clear();
        newer_xyz.find_good_mu();

        newer_xyz.process_good_mu(top_polys);
    }
}

bool PolynomialPairCalculator::generate(long int degree)
{
    const VeryLong zero(0L);
    // Choose a, the product of several small primes and
    // find b, such that c_d b^d = N mod a
    double ALS_MAX = Skewed_config.MAX_ALS();
    double min_als_max = 1000.0;
    d_ = degree;
    // Some bounds
    long int small_prime_limit = Skewed_config.SMALL_PRIME_LIMIT();
    long int large_prime_min = Skewed_config.LARGE_PRIME_MIN();
    size_t max_primes_to_combine = Skewed_config.MAX_PRIMES_TO_COMBINE();
    long int d_power(1L);
    for (size_t i = 0; i < max_primes_to_combine; ++i)
    {
        d_powers_.push_back(d_power);
        d_power *= d_;
    }

    // m_ is approximately
    //     __________
    //    /
    // d /   N / c
    //  /         d
    //
    m_ = (N_ / c_d_).nth_root(d_);

    const VeryLong M = N_.nth_root(d_ + 1);
    //const VeryLong M = N_.nth_root(d_ + 1) / 10L;

    const VeryLong c_d_max = (pow<VeryLong, long int>(M, 2*d_ - 2) / N_).nth_root(d_ - 3);

    c_d_1_max_ = (M * M) / m_;

    c_d_2_max_ = (pow<VeryLong, long int>(M, 2 * d_ - 6) / pow<VeryLong, long int>(m_, d_ - 4)).nth_root(d_ - 2);

    prev_top_c_d_1_ = N_;

    std::ostringstream oss;
    oss << "kleinjung:" << std::endl;
    oss << "c        = " << c_d_ << std::endl;
    oss << " d" << std::endl;
    oss << "_" << std::endl;
    oss << "m        = " << m_ << std::endl;
    oss << "M        = " << M << std::endl;
    oss << "c        = " << c_d_max << std::endl;
    oss << " d,max" << std::endl;
    oss << "c        = " << c_d_1_max_ << std::endl;
    oss << " d-1,max" << std::endl;
    oss << "c        = " << c_d_2_max_ << std::endl;
    oss << " d-2,max" << std::endl;
    if (debug_)
    {
        std::cout << oss.str();
    }

    VeryLong minus_c_d = -c_d_;
    VeryLong minus_c_d_d = minus_c_d * VeryLong(d_);
    minus_c_d_d_d_ = minus_c_d_d.get_double();
    std::vector<VeryLong> c;
    c.resize(d_ + 1);
    c[0] = -N_;
    for (int i = 1; i < d_; i++) c[i] = zero;
    c[d_] = c_d_;
    poly_ = Polynomial<VeryLong>(c);
    if (debug_)
    {
        std::cout << "poly_ = " << poly_ << std::endl;
    }

    VeryLong::generate_prime_table();

    // primes is a vector of small primes satisfying p mod d = 1
    // and such that
    //                d
    //   f(X) =   c  X  - N
    //             d
    //
    // has d roots mod p. Also p does not divide c
    //                                            d

    // for each p in primes_, roots_mod_p_ is a vector containing the d roots of f(X) mod p

    long int p = VeryLong::firstPrime();
    //p = VeryLong::nextPrime();
    while (p < small_prime_limit)
    {
        if (p % d_ == 1L && c_d_ % p != 0L)
        //if (c_d_ % p != 0L)
        {
            std::vector<LongModular> roots;
            find_roots_mod_p<VeryLong, long int, LongModular>(poly_, p, roots);
            //std::cout << "p = " << p << ", roots.size() = " << roots.size() << std::endl;
            //if (roots.size() > 0)
            if (static_cast<long int>(roots.size()) == d_)
            {
                primes_.push_back(p);
                roots_mod_p_.push_back(roots);
            }
        }
        p = VeryLong::nextPrime();
    }

    p = zpnextb(large_prime_min);
    while (p < large_prime_min + 200)
    {
        if (p % d_ == 1L && c_d_ % p != 0L)
        //if (c_d_ % p != 0L)
        {
            std::vector<LongModular> roots;
            find_roots_mod_p<VeryLong, long int, LongModular>(poly_, p, roots);
            std::cout << "p = " << p << ", roots.size() = " << roots.size() << std::endl;
            //if (roots.size() > 0)
            if (static_cast<long int>(roots.size()) == d_)
            {
                primes_.push_back(p);
                roots_mod_p_.push_back(roots);
            }
            else if (static_cast<long int>(roots.size()) == 1L)
            {
                extra_primes_.push_back(p);
            }
        }
        p = zpnext();
    }

    // Now we want to find a and b, satisfying
    //         d
    //     c  b  = N mod a
    //      d
    // because in the end we want to find c  (0 <= i < d) such that
    //                                     i
    //
    //  d
    // ---      i  d-i
    // \    c  b  a    = N
    // /     i
    // ---
    // i=0
    //
    // and also we want b to be close to the d-th root of N / c
    //                                                         d
    // The linear polynomial will be
    // a X - b Y
    //
    //if (debug_)
    {
        std::cout << "primes : " << primes_.size() << std::endl;
        for (size_t i = 0; i < primes_.size(); i++) std::cout << primes_[i] << " ";
        std::cout << std::endl;
    }
    if (primes_.size() < 2)
    {
        return true;
    }

    count_ = 0;
    VeryLong best_b;
    VeryLong best_a;
    VeryLong best_c_d_1;
    double best_als = 1000000.0;
    Polynomial<VeryLong> kj;

    long int primes_to_combine = std::min(max_primes_to_combine, primes_.size());
    for (Combinations combination(primes_to_combine, primes_.size()); !combination.done(); combination.next())
    {
        if (debug_)
        {
            std::cout << "Combination : ";
            combination.display();
        }
        XYZ xyz(*this, combination, primes_to_combine);
        xyz.generate(top_polys_);

        const int extra_primes_to_use = 2;
        if ((int)extra_primes_.size() >= extra_primes_to_use)
        {
            Combinations combination_1(extra_primes_to_use, extra_primes_.size());
            do
                //for (Combinations combination_1(2, extra_primes_.size()); ; combination_1.next())
            {
                VeryLong q(1L);
                for (int i = 0; i < extra_primes_to_use; ++i)
                {
                    q *= extra_primes_[combination_1(i)];
                }
                if (q * xyz.a_ > c_d_1_max_)
                {
                    continue;
                }

                // find root mod q which is sum of
                //  r_i * (q / p_i) * ((q / p_i)^-1 mod p_i)
                VeryLong s(0L);
                for (int i = 0; i < extra_primes_to_use; ++i)
                {
                    VeryLong P_i = q / extra_primes_[combination_1(i)];
                    VeryLong P_i_inv = P_i.inverse(extra_primes_[combination_1(i)]);
                    s += VeryLong(extra_roots_mod_p_[combination_1(i)].get_long()) * P_i * P_i_inv;
                }
                s %= q;
#ifdef CHECK
                {
                    // check
                    //
                    //                 d
                    // should have a  s  = N mod q
                    //              d
                    //
                    VeryLong check = poly_.evaluate(s);
                    if (check % q != zero)
                    {
                        std::cout << "Problem : s doesn't work mod q" << std::endl;
                        std::cout << "s = [" << s << "], q = [" << q << "]" << std::endl;
                        std::cout << "check = " << check << std::endl;
                        std::cout << "check % q = " << check % q << std::endl;
                    }
                }
#endif

                const long int iterations = 30;
                xyz.generate(iterations, q, s, top_polys_);
            }
            while (combination_1.next());
        }
    }

    std::vector<Kleinjung_poly_info> minimized_polys;
    if (!top_polys_.empty())
    {
        std::cout << "Number of polynomials to minimize = " << top_polys_.size() << std::endl;
    }
    for (size_t i = 0; i < top_polys_.size(); i++)
    {
        VeryLong a = top_polys_[i].a_;
        VeryLong b = top_polys_[i].b_;
        Polynomial<VeryLong> fm = top_polys_[i].fm_;
#ifdef CHECK
        {
            VeryLong check = fm.evaluate_homogeneous(b, a);
            const VeryLong zero(0L);
            if (check % N_ != zero)
            {
                std::cout << "Problem: fm(b, a) != 0 % N" << std::endl;
            }
        }
#endif

        double tmp_als = PolynomialOptimizer::average_log_size(fm, 80000);
//      std::cout << "tmp_als = " << tmp_als << ", ALS_MAX = " << ALS_MAX << std::endl;
        if (tmp_als < min_als_max)
        {
            min_als_max = tmp_als;
            std::cout << "min_als_max = " << min_als_max << std::endl;
        }
        if (tmp_als > ALS_MAX) continue;
        std::cout << "tmp_als = " << tmp_als << std::endl;
        std::cout << "fm = " << fm << std::endl;

        VeryLong new_b;
        VeryLong new_m;
        Polynomial<VeryLong> try_poly = fm;
        Polynomial<VeryLong> min_poly;
        VeryLongModular::set_default_modulus(N_);
        VeryLongModular tmp1 = VeryLongModular(b) / VeryLongModular(a);
        VeryLong m = tmp1.get_very_long();
        VeryLong try_m = m;
        VeryLong try_b = b;
        double change = 1.0;
        double prev_als = 100.0;
        double the_average_log_size = 0.0;
        VeryLong s_vl;
        while (change > 0 && fabs(change) > 0.1)
        {
            min_poly = PolynomialOptimizer::minimize_I(try_poly, a, try_b, try_m, s_vl, the_average_log_size, new_b, new_m);
            change = prev_als - the_average_log_size;
            if (change > 0.0)
            {
                std::cout << "Average log size = " << the_average_log_size << std::endl;
                try_poly = min_poly;
                try_m = new_m;
                try_b = new_b;
                prev_als = the_average_log_size;
            }
            else
            {
                the_average_log_size = prev_als;
                min_poly = try_poly;
                new_m = try_m;
                new_b = try_b;
            }
        }

        VeryLong check = min_poly.evaluate_homogeneous(new_b, a);
        const VeryLong zero(0L);
        if (check % N_ != zero)
        {
            std::cout << "Problem: min_poly(b, a) != 0 % N" << std::endl;
        }

        Kleinjung_poly_info pi(a, new_b, min_poly);

        minimized_polys.push_back(pi);

        if (the_average_log_size < best_als)
        {
            best_als = the_average_log_size;
            kj = min_poly;
            best_a = a;
            best_b = new_b;
            std::cout << "best_a = " << best_a << std::endl;
            std::cout << "best_b = " << best_b << std::endl;
            std::cout << "best_als = " << best_als << std::endl;
            std::cout << "kj = " << kj << std::endl;
        }
    }

    if (!top_polys_.empty())
    {
        std::cout << "Number of minimized polynomials = " << minimized_polys.size() << std::endl;
    }
    if (!minimized_polys.empty())
    {
        std::cout << oss.str();
    }
    for (size_t i = 0; i < minimized_polys.size(); ++i)
    {
        VeryLong a = minimized_polys[i].a_;
        VeryLong b = minimized_polys[i].b_;
        Polynomial<VeryLong> kj = minimized_polys[i].fm_;

        std::cout << "a = " << a << std::endl;
        std::cout << "b = " << b << std::endl;

        VeryLongModular::set_default_modulus(N_);
        VeryLongModular tmp1 = VeryLongModular(b) / VeryLongModular(a);
        VeryLong m = tmp1.get_very_long();
        long int s = (long int)PolynomialOptimizer::minimize_I_over_s(Polynomial<VeryLong>::convert_to_double<double>(kj), 1000.0);

        VeryLong s_vl = s;
        double I_F_S = 0.0;
        VeryLong new_b;
        VeryLong new_m;
        Polynomial<VeryLong> new_poly = PolynomialOptimizer::minimize_I(kj, a, b, m, s_vl, I_F_S, new_b, new_m);
        VeryLong check = new_poly.evaluate_homogeneous(new_b, a);
        const VeryLong zero(0L);
        if (check % N_ != zero)
        {
            std::cout << "Problem: new_poly(b, a) != 0 % N" << std::endl;
        }
        double als = PolynomialOptimizer::average_log_size(new_poly, s_vl.get_long());
        std::cout << "kj = " << kj << std::endl;
        std::cout << "f = " << new_poly << std::endl;
        std::cout << "a = " << a << std::endl;
        std::cout << "b = " << new_b << std::endl;
        std::cout << "m = " << new_m << std::endl;
        std::cout << "s = " << s_vl << std::endl;
        std::cout << "als = " << als << std::endl;
        I_F_S = PolynomialOptimizer::average_log_size(new_poly, s_vl);
        double alpha = PolynomialOptimizer::alpha_F(new_poly, 2000, 200);
        double E_F = I_F_S + alpha;
        std::cout << "alpha = " << alpha << std::endl;
        std::cout << "E_F = " << E_F << std::endl;
        std::vector<PolynomialOptimizer::Poly_info> poly_list;
        if (E_F < Skewed_config.PRINTING_BOUND())
        {
            poly_list.push_back(PolynomialOptimizer::Poly_info(new_poly, s, I_F_S, alpha));
        }

        Polynomial<VeryLong> better_poly = PolynomialOptimizer::adjust_root_properties(Skewed_config, new_poly, a, new_b, s_vl, als, poly_list);
        {
            double I_F_S = PolynomialOptimizer::average_log_size(better_poly, s_vl);
            double alpha = PolynomialOptimizer::alpha_F(better_poly, 2000, 200);
            double E_F = I_F_S + alpha;
            std::cout << "alpha = " << alpha << std::endl;
            std::cout << "E_F = " << E_F << std::endl;
        }
        std::sort(poly_list.begin(), poly_list.end());
        std::cout << poly_list.size() << " polynomials to examine ..." << std::endl;
        const size_t max_polys_to_examine = 2000;
        if (poly_list.size() > max_polys_to_examine)
        {
            std::cout << "but only examining the best " << max_polys_to_examine << std::endl;
        }

        size_t examined = 0;
        for (auto iter = poly_list.begin();
                examined < max_polys_to_examine && iter != poly_list.end();
                ++iter)
        {
            VeryLong test_b = new_b;
            Polynomial<VeryLong> test_poly = iter->p;
            std::cout << "Examining " << test_poly << ", (I_F_S = " << iter->I_F_S << ", alpha = " << iter->alpha << ", E_F = " << iter->I_F_S + iter->alpha << ")" << std::endl;
            VeryLong test_s = iter->s;
            VeryLong translated_b;
            VeryLong translated_s;
            Polynomial<VeryLong> translated_poly = PolynomialOptimizer::translate(test_poly, a, test_b, test_s, translated_b, translated_s);
            double translated_I_F_S = PolynomialOptimizer::average_log_size(translated_poly, translated_s);
            double translated_alpha = PolynomialOptimizer::alpha_F(translated_poly, 2000, 200);
            double translated_E_F = translated_I_F_S + translated_alpha;
            if (translated_E_F < Skewed_config.PRINTING_BOUND())
            {
                std::cout << "test poly = " << test_poly << std::endl;
                std::cout << "E(F) = " << iter->I_F_S + iter->alpha << std::endl;
                std::cout << "translated poly = " << translated_poly << std::endl;
                std::cout << "E(F) = " << translated_E_F << std::endl;
                std::cout << "+++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
                output_file_ << translated_poly << std::endl;
                VeryLongModular::set_default_modulus(N_);
                VeryLongModular tmp1 = VeryLongModular(translated_b) / VeryLongModular(a);
                VeryLong m = tmp1.get_very_long();
                output_file_ << "m = " << m << std::endl;
                output_file_ << "a = " << a << std::endl;
                output_file_ << "b = " << translated_b << std::endl;
                output_file_ << "s = " << translated_s << std::endl;
                output_file_ << "alpha = " << translated_alpha << std::endl;
                output_file_ << "E(F) = " << translated_E_F << std::endl;
                output_file_ << translated_poly.evaluate_homogeneous(translated_b, a) << std::endl;
                output_file_ << std::endl;
                output_file_ << "==============================================================================" << std::endl;
                output_file_ << std::flush;
            }
            ++examined;
        }

        std::cout << "better_poly = " << better_poly << std::endl;
        double better_I_F_S = PolynomialOptimizer::average_log_size(better_poly, s_vl);
        alpha = PolynomialOptimizer::alpha_F(better_poly, 2000, 200);
        E_F = better_I_F_S + alpha;
        std::cout << "alpha = " << alpha << std::endl;
        std::cout << "E(F) = " << E_F << std::endl;
        if (E_F < Skewed_config.PRINTING_BOUND())
        {
            output_file_ << better_poly << std::endl;
            VeryLongModular::set_default_modulus(N_);
            VeryLongModular tmp1 = VeryLongModular(new_b) / VeryLongModular(a);
            VeryLong m = tmp1.get_very_long();
            output_file_ << "m = " << m << std::endl;
            output_file_ << "a = " << a << std::endl;
            output_file_ << "b = " << new_b << std::endl;
            output_file_ << "s = " << s_vl << std::endl;
            output_file_ << "alpha = " << alpha << std::endl;
            output_file_ << "E(F) = " << E_F << std::endl;
            output_file_ << better_poly.evaluate_homogeneous(new_b, a) << std::endl;
            output_file_ << std::endl;
            output_file_ << "==============================================================================" << std::endl;
            output_file_ << std::flush;
        }
    }
    if (debug_)
    {
        std::cout << "Exiting kleingjung ..." << std::endl;
    }
    return true;
}

// Skewed polynomial selection (Procedure 5.1.6 in Murphy's thesis)
void skewed_polynomial_selection()
{
    const VeryLong zero(0L);
    timing_file("kleinjung");
    VeryLong N = Skewed_config.N();
    int degree = Skewed_config.DEGREE();
    // start by finding a range for a_d and choosing a_d's that have a lot
    // of small factors
    const double MAX_FRACTION = Skewed_config.MAX_FRACTION();
    double ALS_MAX = Skewed_config.MAX_ALS();

    std::cout << "Skewed polynomial selection procedure starting" << std::endl;
    Skewed_config.display();

    std::fstream Output_file(Skewed_config.OUTPUT_FILE().c_str(), std::ios::out|std::ios::app);

    VeryLong min_ad = Skewed_config.MIN_AD();
    VeryLong max_ad = Skewed_config.MAX_AD();
    double min_log_ad = ln(min_ad);
    double max_log_ad = ln(max_ad);

    // Step 1 - fix a cofactor c which is a product of many small p^k
    //          with log(c) < log(a(d))
    const int NUMBER_OF_PRIMES = 21;

    long int primes[NUMBER_OF_PRIMES] =
    {
        2, 4, 8, 16, 32, 3, 9, 27, 5, 25, 7, 49, 11, 13, 17, 19,
        23, 29, 31, 37, 41
    };
    int finished = 0;
    //long int try_count = 0;
    //long int s; // best s for skewed region of integration
    VeryLong s_vl;
    Polynomial<VeryLong> min_poly;
    VeryLong ad;

    while (!finished)
    {
        //long int p = primes[genrand() % NUMBER_OF_PRIMES];
        VeryLong c = Skewed_config.C_START();

        if (Skewed_config.C_FACTOR() != 0L)
        {
            while (ln(c) < min_log_ad - ln(Skewed_config.C_FACTOR()))
            {
                long int p = primes[genrand() % NUMBER_OF_PRIMES];
                if (ln(c * VeryLong(p)) < max_log_ad) c = c * VeryLong(p);
            }
        }

        VeryLong q = min_ad  / c;
        VeryLong c_inc = c;

        if (q > zero) c *= q;
        // Generate multiples of c in the range

        ad = c;
        if (Skewed_config.C_RESTART() != zero)
        {
            ad = Skewed_config.C_RESTART();
        }

        while (ln(ad) < max_log_ad)
        {
            if (Skewed_config.NON_MONIC())
            {
                PolynomialPairCalculator ppc(N, ad, Output_file, Skewed_config.DEBUG());
                ppc.generate(degree);
            }
            else
            {
                VeryLong tmp = N / ad;
                //try_count++;
                // calculate x = (N - a_d m^d) / m^(d-1)
                // then integer part is a_(d-1) and fractional part ~ a_(d-2) / m
                VeryLong m;
                // calculate an approximation to m accurate to about 16 dec. places
                double md = pow((N.get_double() / ad.get_double()), 1.0 / (double)degree);
                // set up array of powers
                VeryLong m_powers[6];
                double md_powers[6];
                md_powers[0] = 1.0;
                md_powers[1] = md;
                for (int ii = 2; ii < degree + 1; ii++)
                {
                    md_powers[ii] = md * md_powers[ii - 1];
                }
                m_powers[0] = 1L;
                m = VeryLong(md);

                int good_m_found = 0;
                VeryLong quotient;
                VeryLong nleft;
                while (!good_m_found)
                {
                    double err = 0.0;
                    for (int ii = 1; ii < degree + 1; ii++)
                    {
                        m_powers[ii] = m * m_powers[ii - 1];
                    }

                    nleft = N - ad * m_powers[degree]; // N - a_d m^d
                    //cout << "nleft = " << nleft << std::endl;
                    err = nleft.get_double() / (ad.get_double() * (double)degree * md_powers[degree - 1]);
                    const VeryLong zero(0L);
                    if (fabs(err) < Skewed_config.GOOD_M_CUTOFF() ||
                            VeryLong(err) == zero) good_m_found = 1;
                    else
                    {
                        quotient = VeryLong(err);
                        m += quotient;
                        //             std::cout << "err = " << err << ", quotient = " << quotient << ", m = " << m << std::endl;
                    }
                }

                VeryLong remainder = nleft % m_powers[degree - 1];
                quotient = nleft / m_powers[degree - 1];
                VeryLong rr = m_powers[degree - 1] - remainder;
                if (remainder > rr)
                {
                    remainder = rr;
                    quotient += 1L;
                    remainder = -remainder;
                }
                double fraction = remainder.get_double() / md_powers[degree - 1];

                //std::cout << "fraction = " << fraction << std::endl;

                if (fabs(fraction) <= MAX_FRACTION)
                {
                    nleft = remainder;
                    Polynomial<VeryLong> poly = base_m_polynomial_1(N, degree, m);
                    Polynomial<VeryLong> fm = adjust_base_m_polynomial(poly, m);

                    //cout << "quotient = " << quotient << std::endl;
                    //std::cout << "m = " << m << std::endl;
                    //std::cout << "f = " << fm << std::endl;
                    VeryLong new_ad = fm.coefficient(degree);
                    if (new_ad < ad)
                    {
                        std::cerr << "Problem: bad ad " << new_ad << " is < previous ad = " << ad << std::endl;
                    }
                    else ad = new_ad;

                    //try_count = 0;
                    // try adjusting fm by adding a cubic adjustment iX^2(X - m) where i = -1, 0, 1
                    // but only for degree 5 and above.
                    std::vector<VeryLong> ca_coeff;
                    ca_coeff.resize(4);
                    ca_coeff[3] = 1L;
                    ca_coeff[2] = m * VeryLong(-1L);
                    ca_coeff[1] = 0L;
                    ca_coeff[0] = 0L;
                    Polynomial<VeryLong> cubic_adjustment(ca_coeff);
                    long int ca_start = -1;
                    long int ca_end = 2;
                    if (degree > 4) fm = fm - cubic_adjustment;
                    else
                    {
                        ca_start = 0;
                        ca_end = 1;
                    }

                    for (long int ca = ca_start; ca < ca_end; ca++)
                    {
                        double average_log_size = ALS_MAX + 1;

                        VeryLong new_m;
                        Polynomial<VeryLong> try_poly = fm;
                        VeryLong try_m = m;
                        double change = 1.0;
                        double prev_als = 100.0;
                        while (change > 0 && fabs(change) > 0.1)
                        {
                            VeryLong dummy;
                            min_poly = PolynomialOptimizer::minimize_I(try_poly, 1L, try_m, try_m, s_vl, average_log_size, dummy, new_m);
                            change = prev_als - average_log_size;
                            if (change > 0.0)
                            {
                                //std::cout << "Average log size = " << average_log_size << std::endl;
                                try_poly = min_poly;
                                try_m = new_m;
                                prev_als = average_log_size;
                            }
                            else
                            {
                                average_log_size = prev_als;
                                min_poly = try_poly;
                                new_m = try_m;
                            }
                        }
                        //std::cout << "Average log size = " << average_log_size << std::endl;
                        if (average_log_size < ALS_MAX)
                        {
                            // step 3
                            int step3_done = 0;
                            double best_E_F = 100.0;
                            Polynomial<VeryLong> best_step3_poly;
                            VeryLong best_step3_m;
                            VeryLong best_step3_s;
                            double best_step3_alpha = 0.0;
                            while (!step3_done)
                            {
                                Polynomial<VeryLong> better_poly;
                                std::vector<PolynomialOptimizer::Poly_info> poly_list;

                                better_poly = adjust_root_properties_orig(min_poly, new_m, s_vl, average_log_size, poly_list);
                                std::cout << poly_list.size() << " polynomials to examine ..." << std::endl;
                                std::cout << "best E(F) = " << best_E_F << std::endl;
                                double I_F_S = PolynomialOptimizer::average_log_size(better_poly, s_vl);
                                // can we improve size with a final translation?
                                VeryLong better_m;
                                VeryLong better_s;

                                Polynomial<VeryLong> even_better_poly = PolynomialOptimizer::translate(better_poly, 1L, new_m, s_vl, better_m, better_s);
                                double better_I_F_S = PolynomialOptimizer::average_log_size(even_better_poly, better_s);
                                double alpha = PolynomialOptimizer::alpha_F(better_poly, 2000, 200);
                                double E_F = I_F_S + alpha;
                                std::cout << "E(F) = " << E_F << std::endl;

                                double better_alpha = PolynomialOptimizer::alpha_F(even_better_poly, 2000, 200);
                                double better_E_F = better_I_F_S + better_alpha;
                                std::cout << "better E(F) = " << better_E_F << std::endl;
                                step3_done = 1;
                                if (E_F < best_E_F)
                                {
                                    min_poly = better_poly;
                                    best_E_F = E_F;
                                    best_step3_poly = better_poly;
                                    best_step3_m = new_m;
                                    best_step3_s = s_vl;
                                    best_step3_alpha = alpha;
                                    step3_done = 0;
                                }
                                if (better_E_F < best_E_F)
                                {
                                    min_poly = even_better_poly;
                                    best_E_F = better_E_F;
                                    best_step3_poly = even_better_poly;
                                    best_step3_m = better_m;
                                    best_step3_s = better_s;
                                    best_step3_alpha = better_alpha;
                                    new_m = better_m;
                                    s_vl = better_s;
                                    step3_done = 0;
                                }
                                if (best_E_F > Skewed_config.REPEAT_CUTOFF()) step3_done = 1;
                                if (best_E_F < Skewed_config.PRINTING_BOUND())
                                {
                                    Output_file << best_step3_poly << std::endl;
                                    Output_file << "m = " << best_step3_m << std::endl;
                                    Output_file << "s = " << best_step3_s << std::endl;
                                    Output_file << "alpha = " << best_step3_alpha << std::endl;
                                    Output_file << "E(F) = " << best_E_F << std::endl;
                                    Output_file << best_step3_poly.evaluate(new_m) << std::endl;
                                    Output_file << std::endl;
                                    Output_file << "==============================================================================" << std::endl;
                                    Output_file << std::flush;
                                }
                            }
                        }
                        fm = fm + cubic_adjustment;
                    }
                }
                //cout << "c_inc = " << c_inc << std::endl;
                //cout << "ad = " << ad << std::endl;
            }
            ad += c_inc;
        }
        if (Skewed_config.C_FACTOR() == 0L) finished = 1;
    }
}
