#ifndef POLYNOMIALOPTIMIZER_H
#define POLYNOMIALOPTIMIZER_H
#include "Polynomial.h"
#include "VeryLong.h"
#include "Config.h"

namespace PolynomialOptimizer
{
const double MAX_SAMPLE_RANGE = 1e15;
const int SAMPLE_SIZE = 100;

VeryLong b0_vl (const Polynomial<VeryLong>& f, const VeryLong& a, const VeryLong& b,const VeryLong& c0, const VeryLong& c1, const VeryLong& t);

// f(x_t) + (c1*(x - t) + c0)(a * (x - t) - b),
VeryLong b1_vl (const Polynomial<VeryLong>& f, const VeryLong& a, const VeryLong& b, const VeryLong& c0, const VeryLong& c1, const VeryLong& t);

VeryLong b2_vl (const Polynomial<VeryLong>& f, const VeryLong& a, const VeryLong& c1, const VeryLong& t);

VeryLong b3_vl (const Polynomial<VeryLong>& f, const VeryLong& t);

VeryLong b4_vl (const Polynomial<VeryLong>& f, const VeryLong& t);

VeryLong b5_vl (const Polynomial<VeryLong>& f);

VeryLong b_k_vl (const Polynomial<VeryLong>& f, int k, const VeryLong& t);

template <typename DOUBLE> DOUBLE b0 (const Polynomial<DOUBLE>& f, DOUBLE t);

template <typename DOUBLE> DOUBLE b0 (const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE b, DOUBLE c0, DOUBLE c1, DOUBLE t);

template <typename DOUBLE> DOUBLE b1 (const Polynomial<DOUBLE>& f, DOUBLE t);

template <typename DOUBLE> DOUBLE b1 (const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE b, DOUBLE c0, DOUBLE c1, DOUBLE t);

template <typename DOUBLE> DOUBLE b2 (const Polynomial<DOUBLE>& f, DOUBLE t);

template <typename DOUBLE> DOUBLE b2 (const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE c1, DOUBLE t);

template <typename DOUBLE> DOUBLE b3 (const Polynomial<DOUBLE>& f, DOUBLE t);

template <typename DOUBLE> DOUBLE b4 (const Polynomial<DOUBLE>& f);

template <typename DOUBLE> DOUBLE b4 (const Polynomial<DOUBLE>& f, DOUBLE t);

template <typename DOUBLE> DOUBLE b5 (const Polynomial<DOUBLE>& f);

template <typename DOUBLE> DOUBLE b_k (const Polynomial<DOUBLE>& f, int k, DOUBLE t);

template <typename DOUBLE> DOUBLE minimize_I_over_t(const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE b, DOUBLE s);

template <typename DOUBLE> DOUBLE minimize_I_over_t(const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE b, DOUBLE& c0, DOUBLE& c1, DOUBLE s);
template <typename DOUBLE> DOUBLE I(const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE b, DOUBLE& c0, DOUBLE& c1, DOUBLE t, DOUBLE s);

template <typename DOUBLE> DOUBLE J(const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE bb, DOUBLE c0, DOUBLE c1, DOUBLE t, DOUBLE s);

template <typename DOUBLE> DOUBLE minimize_I_over_t_and_s(const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE b, DOUBLE& c0, DOUBLE& c1, DOUBLE& t, DOUBLE& s);

// General J integral for arbitrary degree d polynomial.
// For degree d, J = sum_{k=0,2,...,2d} W(k) * s^(k-d) * S_k
// where W(k) = 4/((k+1)*(2d-k+1)) and S_k = sum_{i+j=k} b[i]*b[j]
template <typename DOUBLE>
DOUBLE J_general(const DOUBLE* b, int d, DOUBLE s)
{
    DOUBLE result = 0.0;
    // Precompute s powers: s_pow[k] = s^(k-d) for k = 0, 2, 4, ..., 2d
    // We build them incrementally: s^(k-d) = s^((k-2)-d) * s^2
    DOUBLE s2 = s * s;
    DOUBLE s_pow = 1.0;
    // Start with s^(-d)
    for (int i = 0; i < d; i++) s_pow /= s;
    // s_pow is now s^(-d) = s^(0 - d)

    for (int k = 0; k <= 2 * d; k += 2)
    {
        // Compute S_k = sum_{i+j=k, 0<=i<=d, 0<=j<=d} b[i]*b[j]
        int i_min = (k > d) ? k - d : 0;
        int i_max = (k < d) ? k : d;
        DOUBLE S_k = 0.0;
        for (int i = i_min; i <= i_max; i++)
        {
            int j = k - i;
            if (i == j)
                S_k += b[i] * b[j];
            else if (i < j)
                S_k += 2.0 * b[i] * b[j];
        }
        // W(k) = 4.0 / ((k+1) * (2*d - k + 1))
        DOUBLE W_k = 4.0 / ((DOUBLE)(k + 1) * (DOUBLE)(2 * d - k + 1));
        result += W_k * s_pow * S_k;
        s_pow *= s2;
    }
    return result;
}

template <typename DOUBLE>
DOUBLE J(const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE bb, DOUBLE c0, DOUBLE c1, DOUBLE t, DOUBLE s)
{
    int d = f.deg();
    std::vector<DOUBLE> b(d + 1);
    b[0] = b0(f,a,bb,c0,c1,t);
    b[1] = b1(f,a,bb,c0,c1,t);
    b[2] = b2(f,a,c1,t);
    for (int i = 3; i <= d; i++)
    {
        b[i] = b_k(f, i, t);
    }

    return J_general(b.data(), d, s);
}

template <typename DOUBLE>
DOUBLE J(const Polynomial<DOUBLE>& f, DOUBLE s)
{
    int d = f.deg();
    std::vector<DOUBLE> b(d + 1);
    for (int i = 0; i <= d; i++)
    {
        b[i] = f.coefficient(i);
    }

    return J_general(b.data(), d, s);
}

template <typename DOUBLE>
DOUBLE minimize_I_over_s_1(const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE b, DOUBLE& c0, DOUBLE& c1, DOUBLE t)
{
    DOUBLE s = 1000.0;
    DOUBLE last_s = 0.0;
    DOUBLE s_min = 1;
    DOUBLE s_max = MAX_SAMPLE_RANGE;
    DOUBLE s_sample_size = SAMPLE_SIZE;
    DOUBLE s_delta = (s_max - s_min) / s_sample_size;
    DOUBLE s_try = s_min;
    DOUBLE s_diff = 0.0;
    DOUBLE value;
    DOUBLE min_value = J(f, a, b, c0, c1, t, s);
    // min_value = I(f, m, c0, c1, t, s);
    // return s;
    int done = 0;
    //int iterations = 0;

    while (!done)
    {
        s_try = s_min;
        for (int l = 0; l < s_sample_size; l++)
        {
            value = I(f, a, b, c0, c1, t, s_try);
            if (value < min_value)
            {
                min_value = value;
                last_s = s;
                s = s_try;
            }
            s_try += s_delta;
        }
        s_min = s - s_sample_size*s_delta / 3.0;
        if (s_min <= (DOUBLE)0.0) s_min = 1;
        s_max = s + s_sample_size*s_delta / 3.0;
        s_delta = (s_max - s_min) / s_sample_size;
        s_diff = last_s - s;
        if (s_diff < (DOUBLE)0.0) s_diff = -s_diff;
        last_s = s;

        const DOUBLE epsilon = 1.0;
        //iterations++;
        if (s_diff < epsilon && s_delta < (DOUBLE)1.0) done = 1;
    }
    VeryLong s_vl(s);
    value = J<DOUBLE>(f, a, b, c0, c1, t, s_vl);
    DOUBLE value1 = J<DOUBLE>(f, a, b, c0, c1, t, s_vl + VeryLong(1.0));
    if (value < value1) return s_vl;
    else return (s_vl + VeryLong(1.0));
}

template <typename DOUBLE>
DOUBLE minimize_I_over_s(const Polynomial<DOUBLE>& f, DOUBLE s_start)
{
    DOUBLE s = s_start;
    DOUBLE last_s = 0.0;
    DOUBLE s_min = 1;
    DOUBLE s_max = MAX_SAMPLE_RANGE;
    DOUBLE s_sample_size = SAMPLE_SIZE;
    DOUBLE s_delta = (s_max - s_min) / s_sample_size;
    DOUBLE s_try = s_min;
    DOUBLE s_diff = 0.0;
    DOUBLE value;
    DOUBLE min_value = J(f, s);
    // return s;
    int done = 0;
    int iterations = 0;

    while (!done)
    {
        s_try = s_min;
        for (int l = 0; l < s_sample_size; l++)
        {
            value = J(f, s_try);
            if (value < min_value)
            {
                min_value = value;
                last_s = s;
                s = s_try;
            }
            s_try += s_delta;
        }
        s_min = s - s_sample_size*s_delta / 3.0;
        if (s_min <= (DOUBLE)0.0) s_min = 1;
        s_max = s + s_sample_size*s_delta / 3.0;
        s_delta = (s_max - s_min) / s_sample_size;
        s_diff = fabs(last_s - s);

        last_s = s;
        if (iterations > 1000)
        {
            std::cout << s_diff << std::endl;
        }
        const DOUBLE epsilon = 1.0;
        iterations++;
        if ((f.deg() < 5 && s_diff == (DOUBLE)0.0) || (s_diff < epsilon && s_delta < (DOUBLE)1.0)) done = 1;
    }
    VeryLong s_vl(s);
    value = J<DOUBLE>(f, s_vl);
    DOUBLE value1 = J<DOUBLE>(f, s_vl + VeryLong(1.0));
    if (value < value1) return s_vl;
    else return (s_vl + VeryLong(1.0));
}

template <typename DOUBLE>
DOUBLE minimize_I_over_s(const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE b,
                         DOUBLE c0, DOUBLE c1, DOUBLE t, DOUBLE s_start)
{
    DOUBLE s = s_start;
    DOUBLE last_s = 0.0;
    DOUBLE s_min = 1;
    DOUBLE s_max = MAX_SAMPLE_RANGE;
    DOUBLE s_sample_size = SAMPLE_SIZE;
    DOUBLE s_delta = (s_max - s_min) / s_sample_size;
    DOUBLE s_try = s_min;
    DOUBLE s_diff = 0.0;
    DOUBLE value;
    DOUBLE min_value = J(f, a, b, c0, c1, t, s);
    // return s;
    int done = 0;
    int iterations = 0;

    while (!done)
    {
        s_try = s_min;
        for (int l = 0; l < s_sample_size; l++)
        {
            value = J(f, a, b, c0, c1, t, s_try);
            if (value < min_value)
            {
                min_value = value;
                last_s = s;
                s = s_try;
            }
            s_try += s_delta;
        }
        s_min = s - s_sample_size*s_delta / 3.0;
        if (s_min <= (DOUBLE)0.0) s_min = 1;
        s_max = s + s_sample_size*s_delta / 3.0;
        s_delta = (s_max - s_min) / s_sample_size;
        s_diff = fabs(last_s - s);

        last_s = s;
        if (iterations > 1000)
        {
            std::cout << s_diff << std::endl;
        }
        const DOUBLE epsilon = 1.0;
        iterations++;
        if ((f.deg() < 5 && s_diff == (DOUBLE)0.0) || (s_diff < epsilon && s_delta < (DOUBLE)1.0)) done = 1;
    }
    VeryLong s_vl(s);
    value = J<DOUBLE>(f, a, b, c0, c1, t, s_vl);
    DOUBLE value1 = J<DOUBLE>(f, a, b, c0, c1, t, s_vl + VeryLong(1.0));
    if (value < value1) return s_vl;
    else return (s_vl + VeryLong(1.0));
}

template <typename DOUBLE>
Polynomial<VeryLong> minimize_I(const Polynomial<VeryLong>& fm,
                                const VeryLong& a1,
                                const VeryLong& b1,
                                const VeryLong& m1,
                                VeryLong& best_s,
                                DOUBLE& I_F_S,
                                VeryLong& new_b,
                                VeryLong& new_m)
{
    bool debug = false;
    if (debug)
    {
        std::cout << "1. fm = " << fm << std::endl;
        std::cout << "2. a1 = " << a1 << std::endl;
        std::cout << "3. b1 = " << b1 << std::endl;
        std::cout << "4. m1 = " << m1 << std::endl;
        std::cout << "5. best_s = " << best_s << std::endl;
    }
    Polynomial<DOUBLE> f = Polynomial<VeryLong>::convert_to_double<DOUBLE>(fm);
    if (debug)
        std::cout << "a. f = " << f << std::endl;
    DOUBLE a = a1;
    DOUBLE b = b1;

    DOUBLE c0 = 0.0;
    DOUBLE c1 = 0.0;
    DOUBLE t = 0.0;
    DOUBLE s = 1000.0;
    // s = 1;
    int finished = 0;
    DOUBLE min_value = J(f, a, b, c0, c1, t, s); // initial value
    if (debug)
        std::cout << "b. min_value = " << min_value << std::endl;
    DOUBLE local_minimum = min_value;
    DOUBLE min_c0;
    DOUBLE min_c1;
    DOUBLE min_t = t;
    DOUBLE min_s = s;
    int iterations = 0;
    min_t = minimize_I_over_t(f, a, b, min_c0, min_c1, s);
    if (debug)
        std::cout << "c. min_t = " << min_t << std::endl;
    min_s = minimize_I_over_s_1(f, a, b, min_c0, min_c1, min_t);
    if (debug)
        std::cout << "d. min_s = " << min_s << std::endl;
    min_value = I(f, a, b, min_c0, min_c1, min_t, min_s);
    if (debug)
        std::cout << "d. min_value = " << min_value << std::endl;
    if (min_value < local_minimum)
    {
        c0 = min_c0;
        c1 = min_c1;
        t = min_t;
        s = min_s;
        local_minimum = min_value;
    }
    else
    {
        min_t = t;
        min_s = s;
    }
    if (debug)
        std::cout << "(c0,c1,t,s) = (" << c0 << "," << c1 << "," << t << "," << s << "), local_minimum = " << local_minimum << ", iterations = " << iterations << std::endl;
    while (!finished && iterations < 20)
    {
        min_value = minimize_I_over_t_and_s(f, a, b, min_c0, min_c1, min_t, min_s);
        if (min_value < local_minimum)
        {
            local_minimum = min_value;
            c0 = min_c0;
            c1 = min_c1;
            t = min_t;
            s = min_s;
        }
        else finished = 1;
        if (debug)
            std::cout << "(c0,c1,t,s) = (" << c0 << "," << c1 << "," << t << "," << s << "), local_minimum = " << local_minimum << ", iterations = " << iterations << std::endl;
        iterations++;
    }
    if (debug)
        std::cout << "(c0,c1,t,s) = (" << c0 << "," << c1 << "," << t << "," << s << "), local_minimum = " << local_minimum << ", iterations = " << iterations << std::endl;

    // check ...
    // Now take nearest integers to c0, c1, t, and s
    // recalculate s for these values
    VeryLong t_vl(t);
    VeryLong c0_int = VeryLong(c0);
    VeryLong c1_int = VeryLong(c1);
    DOUBLE s_double = minimize_I_over_s<DOUBLE>(f, a, b, c0_int, c1_int, t_vl, s);
    VeryLong s_vl(s_double);
    //cout << "s_double = " << s_double << endl;

    //min_value = J(f, m, c0_int.get_double(), c1_int.get_double(), t_int, s_int);
    min_value = J<DOUBLE>(f, a, b, c0_int, c1_int, t_vl, s_vl);
    //cout << "(c0,c1,t,s) = (" << c0_int << "," << c1_int << "," << t_int << "," << s_int << "), min_value = " << min_value << endl;
    if (debug)
        std::cout << "(c0,c1,t,s) = (" << c0_int << "," << c1_int << "," << t_vl << "," << s_vl << "), min_value = " << min_value << std::endl;

    std::vector<VeryLong> bb;
    int deg = fm.deg();
    bb.resize(deg + 1);
    bb[0] = b0_vl(fm,a1,b1,c0_int,c1_int,t_vl);
    bb[1] = b1_vl(fm,a1,b1,c0_int,c1_int,t_vl);
    bb[2] = b2_vl(fm,a1,c1_int,t_vl);
    for (int i = 3; i <= deg; i++)
    {
        bb[i] = b_k_vl(fm, i, t_vl);
    }

    Polynomial<VeryLong> min_poly(bb);

    //std::cout << "min poly = " << min_poly << std::endl;
    //std::cout << "new value for b = " << b1 + a1 * t_vl << std::endl;
    //std::cout << "new value for m = " << m1 + t_vl << std::endl;
    //std::cout << "s = " << s_vl << std::endl;
    //cout << min_poly.evaluate(m1 + t_int) << std::endl;
    //cout << "alpha = " << alpha_F(min_poly, 100,50) << std::endl;
    best_s = s_vl;
    I_F_S = log(sqrt(min_value/4.0));
    new_m = m1 + t_vl;
    new_b = b1 + a1 * t_vl;
    if (debug)
    {
        std::cout << "6. min_poly = " << min_poly << std::endl;
        std::cout << "7. a1 = " << a1 << std::endl;
        std::cout << "8. new_b = " << new_b << std::endl;
        std::cout << "9. new_m = " << new_m << std::endl;
        std::cout << "10. best_s = " << best_s << std::endl;
    }
    return min_poly;
}
// b0 - b5 are coefficients of adjusted f(x_t) + (c1*x_t + c0)(x_t - m)
// where x_t = x - t
// Note that at this = N at x_t = m, since f(m) = N, i.e. x = m + t
//
// We also have versions for f(x_t) + (c1*x_t + c0)(a * x_t - b),
// so at x_t = b / a (mod N), this = f(b/a) = 0, so x - t = b / a mod N
// i.e. x = b/a + t, or a * x = b + a * t
// so if g(x) = f(x_t) + (c1*x_t + c0)(a*x_t - b)
// g(b/a + t) = 0 mod N
// or
// g(b + a * t) = 0 mod N,
// i.e. adjusting m by t <=> adjusting b by a * t

template <typename DOUBLE>
DOUBLE b0 (const Polynomial<DOUBLE>& f, DOUBLE t)
{
    DOUBLE tt = 1.0;
    DOUBLE res = 0.0;
    res = f.coefficient(0);

    for (int i = 1; i <= f.deg(); i++)
    {
        tt *= -t;
        res += f.coefficient(i) * tt;
    }
    return res;
}

template <typename DOUBLE>
DOUBLE b0 (const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE b, DOUBLE c0, DOUBLE c1, DOUBLE t)
{
    DOUBLE res = b0(f, t);
    res -= b * c0;
    res -= t * c0 * a;
    res += b * c1 * t;
    res += t * c1 * t * a;
    return res;
}

template <typename DOUBLE>
DOUBLE b1 (const Polynomial<DOUBLE>& f, DOUBLE t)
{
    DOUBLE tt = 1.0;
    DOUBLE res = 0.0;
    long int multiplier = 1;
    res = f.coefficient(1);

    for (int i = 2; i <= f.deg(); i++)
    {
        multiplier++;
        tt *= -t;
        res += f.coefficient(i) * tt * multiplier;
    }
    return res;
}

template <typename DOUBLE>
DOUBLE b1 (const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE b, DOUBLE c0, DOUBLE c1, DOUBLE t)
{
    DOUBLE res = b1(f, t);
    res -= b * c1;
    res += c0 * a;
    res -= 2.0 * c1 * t * a;
    return res;
}

template <typename DOUBLE>
DOUBLE b2 (const Polynomial<DOUBLE>& f, DOUBLE t)
{
    DOUBLE tt = 1.0;
    DOUBLE res = 0.0;
    long int multiplier = 1;
    int increment = 2;
    res = f.coefficient(2);

    for (int i = 3; i <= f.deg(); i++)
    {
        //cout << sign << " : " << multiplier << endl;
        multiplier += increment;
        increment++;
        tt *= -t;
        res += f.coefficient(i) * tt * multiplier;
    }
    return res;
}

template <typename DOUBLE>
DOUBLE b2 (const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE c1, DOUBLE t)
{
    DOUBLE res = b2(f, t);
    res += c1 * a;
    return res;
}

template <typename DOUBLE>
DOUBLE b3 (const Polynomial<DOUBLE>& f, DOUBLE t)
{
    DOUBLE tt = 1.0;
    DOUBLE res = 0.0;
    long int multiplier = 1;
    int increment = 3;
    int increment1 = 3;
    res = f.coefficient(3);

    for (int i = 4; i <= f.deg(); i++)
    {
        //cout << sign << " : " << multiplier << endl;
        multiplier += increment;
        increment += increment1;
        increment1++;
        tt *= -t;
        res += f.coefficient(i) * tt * multiplier;
    }
    return res;
}

template <typename DOUBLE>
DOUBLE b4 (const Polynomial<DOUBLE>& f)
{
    if (f.deg() < 4) return 0.0;
    return f.coefficient(4);
}

template <typename DOUBLE>
DOUBLE b4 (const Polynomial<DOUBLE>& f, DOUBLE t)
{
    if (f.deg() < 4) return 0.0;
    return b_k(f, 4, t);
}

template <typename DOUBLE>
DOUBLE b5 (const Polynomial<DOUBLE>& f)
{
    if (f.deg() < 5) return 0.0;
    if (f.deg() == 5) return f.coefficient(5);
    // For degree > 5, b5 with t=0 is just the coefficient
    return f.coefficient(5);
}

// General b_k: k-th coefficient of f(x-t) when expanded in powers of x.
// i.e. b_k = sum_{i=k}^{d} C(i,k) * f_i * (-t)^(i-k)
// where f(x) = sum_{i=0}^{d} f_i * x^i and f(x-t) = sum_{k=0}^{d} b_k * x^k
template <typename DOUBLE>
DOUBLE b_k (const Polynomial<DOUBLE>& f, int k, DOUBLE t)
{
    int d = f.deg();
    if (d < k) return 0.0;
    if (d == k) return f.coefficient(k);
    // Compute sum_{i=k}^{d} C(i,k) * f_i * (-t)^(i-k)
    DOUBLE res = f.coefficient(k);
    DOUBLE neg_t_power = 1.0;
    DOUBLE binom = 1.0; // C(k,k) = 1, will build up C(i,k)
    for (int i = k + 1; i <= d; i++)
    {
        neg_t_power *= -t;
        binom *= (DOUBLE)i / (DOUBLE)(i - k);
        res += (DOUBLE)binom * f.coefficient(i) * neg_t_power;
    }
    return res;
}


// if partial derivatives of I wrt c0 and c1 are zero, then we get
// two simultaneous equations for c0 and c1 of the form
// a c0 + b c1 = d
// c c0 + d c1 = f
// these functions define the coefficients a - f
template <typename DOUBLE>
DOUBLE a(DOUBLE aa, DOUBLE bb, DOUBLE t, DOUBLE s)
{
    DOUBLE res = bb +  aa * t;
    res *= res;
    res *= 8.0/11.0;
    res /= s*s*s*s*s;
    res += (8.0 * aa * aa)/ (27.0 * s * s * s);
    //std::cout << "a() : (aa,bb,t,s) = (" << std::setprecision(20) << aa << "," << bb << "," << t << "," << s << "), res = " << res << std::endl;
    return res;
}

template <typename DOUBLE>
DOUBLE b(DOUBLE aa, DOUBLE bb, DOUBLE t, DOUBLE s)
{
    DOUBLE res = bb + aa * t;
    res *= res;
    res *= - t * 8.0 / 11.0;
    res /= s*s*s*s*s;
    DOUBLE res1 = 2.0*bb + 3.0*aa*t;
    res1 *= -(8.0 * aa)/ (27.0 * s * s * s);
    res += res1;
    //std::cout << "b() : (aa,bb,t,s) = (" << aa << "," << bb << "," << t << "," << s << "), res = " << res << std::endl;
    return res;
}

template <typename DOUBLE>
DOUBLE c(DOUBLE aa, DOUBLE bb, DOUBLE t, DOUBLE s)
{
    DOUBLE res = bb + aa * t;
    res *= res;
    res *= - t * 8.0 / 11.0;
    res /= s*s*s*s*s;
    DOUBLE res1 = 2.0*bb + 3.0*aa*t;
    res1 *= -(8.0 * aa)/ (27.0 * s * s * s);
    res += res1;
    //std::cout << "c() : (aa,bb,t,s) = (" << aa << "," << bb << "," << t << "," << s << "), res = " << res << std::endl;
    return res;
}

template <typename DOUBLE>
DOUBLE d(DOUBLE aa, DOUBLE bb, DOUBLE t, DOUBLE s)
{
    DOUBLE res = bb + aa * t;
    res *= res;
    res *= 8.0 * t * t / (11.0 * s * s * s * s * s);
    DOUBLE s3 = s * s * s;
    DOUBLE res1 = aa * t * (bb + aa * t) * 16.0 / (27.0 * s3);
    res += res1;
    res1 = bb + 2.0*aa*t;
    res1 *= res1;
    res1 *= 8.0 / (27.0 * s3);
    res += res1;
    res += (8.0 * aa * aa)/(35.0 * s);
    //std::cout << "d() : (aa,bb,t,s) = (" << aa << "," << bb << "," << t << "," << s << "), res = " << res << std::endl;
    return res;
}

// functions to calculate optimum values of c0 and c1 for given s and t
// we have b0 = A + c1 t(m+t) - c0 (m+t)
//         b1 = B + c0 - (m+2t)c1
//         b2 = C + c1

template <typename DOUBLE>
DOUBLE A (const Polynomial<DOUBLE>& f, DOUBLE t)
{
    DOUBLE tt = 1.0;
    DOUBLE res = 0.0;
    res = f.coefficient(0);

    for (int i = 1; i <= f.deg(); i++)
    {
        tt *= -t;
        res += f.coefficient(i) * tt;
    }
    return res;
}

template <typename DOUBLE>
DOUBLE B (const Polynomial<DOUBLE>& f, DOUBLE t)
{
    DOUBLE tt = 1.0;
    DOUBLE res = 0.0;
    long int multiplier = 1;
    res = f.coefficient(1);

    for (int i = 2; i <= f.deg(); i++)
    {
        multiplier++;
        tt *= -t;
        res += f.coefficient(i) * tt * multiplier;
    }
    return res;
}

template <typename DOUBLE>
DOUBLE C (const Polynomial<DOUBLE>& f, DOUBLE t)
{
    DOUBLE tt = 1.0;
    DOUBLE res = 0.0;
    long int multiplier = 1;
    int increment = 2;
    res = f.coefficient(2);

    for (int i = 3; i <= f.deg(); i++)
    {
        //cout << sign << " : " << multiplier << endl;
        multiplier += increment;
        increment++;
        tt *= -t;
        res += f.coefficient(i) * tt * multiplier;
    }
    return res;
}

template <typename DOUBLE>
DOUBLE e(const Polynomial<DOUBLE>&f, DOUBLE aa, DOUBLE bb, DOUBLE t, DOUBLE s)
{
    DOUBLE res = bb + aa * t;
    res *= A(f, t) * 8.0 / (11.0 * s * s * s * s * s);
    DOUBLE s3 = s * s * s;
    DOUBLE res1 = -8.0 * aa * B(f, t) / (27.0 * s3);
    res += res1;
    res1 = 8.0 * (bb + aa * t) * C(f,t) / (27.0 * s3);
    res += res1;
    res1 = 8.0 * (bb + aa * t) * b4(f, t) / (35.0 * s);
    res += res1;
    res1 = -8.0 * aa * b3(f, t) / (35.0 * s);
    res += res1;
    res1 = -8.0 * aa * b5(f) * s / 35.0;
    res += res1;
    return res;
}

template <typename DOUBLE>
DOUBLE ff(const Polynomial<DOUBLE>&f, DOUBLE aa, DOUBLE bb, DOUBLE t, DOUBLE s)
{
    DOUBLE res = -8.0 * t * (bb + aa * t) * A(f, t) / (11.0 * s * s * s * s * s);
    DOUBLE s3 = s * s * s;
    res += -8.0 * t * (bb + aa * t) * C(f, t) / (27.0 * s3);
    res += -8.0 * aa * A(f, t) / (27.0 * s3);
    res += 8.0 * (bb + 2.0*aa*t) * B(f, t) / (27.0 * s3);
    res += -8.0 * aa * C(f, t) / (35.0 * s);
    res += -8.0 * t * (bb + aa * t) * b4(f, t) / (35.0 * s);
    res += 8.0 * (bb + 2.0*aa*t) * b3(f, t) / (35.0 * s);
    res += -8.0 * aa * b4(f, t) * s / 35.0;
    res += 8.0 * (bb + 2.0*aa*t) * b5(f) * s / 35.0;
    return res;
}

template <typename DOUBLE>
void best_c0_and_c1(const Polynomial<DOUBLE>& f, DOUBLE aa, DOUBLE bb, DOUBLE t, DOUBLE s, DOUBLE& c0, DOUBLE& c1)
{
    // Compute optimal c0, c1 by solving the 2x2 system arising from
    // dJ/dc0 = 0, dJ/dc1 = 0.
    //
    // The linear polynomial adjustment (c1*(x-t) + c0)(a*(x-t) - b) modifies:
    //   b0 -> b0 - (b+at)*c0 + (b+at)*c1*t
    //   b1 -> b1 + a*c0 - (b+2at)*c1
    //   b2 -> b2 + a*c1
    //
    // J is quadratic in (c0,c1) so we solve the linear system.
    // Let p0 = -(b+at), q0 = (b+at)*t   [db0/dc0, db0/dc1]
    // Let p1 = a,        q1 = -(b+2at)   [db1/dc0, db1/dc1]
    // Let p2 = 0,        q2 = a           [db2/dc0, db2/dc1]

    int d = f.deg();
    DOUBLE bat = bb + aa * t;
    DOUBLE bat2 = bb + 2.0 * aa * t;
    DOUBLE p0 = -bat, q0 = bat * t;
    DOUBLE p1 = aa,   q1 = -bat2;
    DOUBLE p2 = 0.0,  q2 = aa;

    // Compute base coefficients (without c0, c1 contribution)
    std::vector<DOUBLE> b_base(d + 1);
    b_base[0] = b0(f, t);
    b_base[1] = b1(f, t);
    b_base[2] = b2(f, t);
    for (int i = 3; i <= d; i++)
    {
        b_base[i] = b_k(f, i, t);
    }

    // J = sum_{k=0,2,...,2d} W(k) * s^(k-d) * S_k
    // where S_k = sum_{i+j=k} b[i]*b[j]
    // b[0] = b_base[0] + p0*c0 + q0*c1
    // b[1] = b_base[1] + p1*c0 + q1*c1
    // b[2] = b_base[2] + p2*c0 + q2*c1
    // b[i] = b_base[i] for i >= 3
    //
    // dJ/dc0 = sum_k W(k)*s^(k-d) * sum_{i+j=k} 2*(db[i]/dc0)*b[j]
    //        = 2 * sum_k W(k)*s^(k-d) * sum_{i+j=k, i<=2} p_i * b[j]
    // At optimum with c0=c1=0, we get a linear system:
    //   A00*c0 + A01*c1 = -E0
    //   A10*c0 + A11*c1 = -E1
    // where:
    //   A00 = 2*sum_k W(k)*s^(k-d) * sum_{i+j=k,i<=2} p_i*p_j (with p_j=0 for j>2)
    //         + terms involving p_i*b_base[j] cross-contributions...
    //
    // Actually, since J is quadratic in c0,c1, we can write:
    //   J = J0 + (gradient)^T * [c0,c1] + 0.5 * [c0,c1]^T * H * [c0,c1]
    // and the minimum is at H*[c0,c1] = -gradient (at c0=c1=0)
    //
    // Compute the Hessian and gradient numerically for correctness:
    // H[i][j] = d²J/dc_i dc_j, gradient[i] = dJ/dc_i at c0=c1=0

    // Precompute s-power weights: w[k] = W(k) * s^(k-d) for even k from 0 to 2d
    std::vector<DOUBLE> w(2 * d + 1, 0.0);
    DOUBLE s_pow = 1.0;
    for (int i = 0; i < d; i++) s_pow /= s;
    DOUBLE s2 = s * s;
    for (int k = 0; k <= 2 * d; k += 2)
    {
        w[k] = 4.0 / ((DOUBLE)(k + 1) * (DOUBLE)(2 * d - k + 1)) * s_pow;
        s_pow *= s2;
    }

    // Gradient: dJ/dc0 and dJ/dc1 evaluated at c0=c1=0
    // dJ/dc0 = 2 * sum_{even k} w[k] * sum_{i+j=k, 0<=i<=2, 0<=j<=d} p_i * b_base[j]
    // dJ/dc1 = 2 * sum_{even k} w[k] * sum_{i+j=k, 0<=i<=2, 0<=j<=d} q_i * b_base[j]
    DOUBLE grad0 = 0.0, grad1 = 0.0;
    DOUBLE p[3] = {p0, p1, p2};
    DOUBLE q[3] = {q0, q1, q2};
    for (int k = 0; k <= 2 * d; k += 2)
    {
        DOUBLE g0_k = 0.0, g1_k = 0.0;
        for (int i = 0; i <= 2 && i <= k; i++)
        {
            int j = k - i;
            if (j >= 0 && j <= d)
            {
                g0_k += p[i] * b_base[j];
                g1_k += q[i] * b_base[j];
            }
        }
        // No second symmetric sum needed here; grad*2 below already accounts for i/j symmetry.
        grad0 += w[k] * g0_k;
        grad1 += w[k] * g1_k;
    }
    grad0 *= 2.0;
    grad1 *= 2.0;
    // Hessian: H[a][b] = 2 * sum_{even k} w[k] * sum_{i+j=k, i<=2, j<=2} (p_i if a=0, q_i if a=1) * (p_j if b=0, q_j if b=1)
    DOUBLE H00 = 0.0, H01 = 0.0, H11 = 0.0;
    for (int k = 0; k <= 2 * d; k += 2)
    {
        DOUBLE h00_k = 0.0, h01_k = 0.0, h11_k = 0.0;
        for (int i = 0; i <= 2 && i <= k; i++)
        {
            int j = k - i;
            if (j >= 0 && j <= 2)
            {
                h00_k += p[i] * p[j];
                h01_k += p[i] * q[j];
                h11_k += q[i] * q[j];
            }
        }
        H00 += w[k] * h00_k * 2.0;
        H01 += w[k] * h01_k * 2.0;
        H11 += w[k] * h11_k * 2.0;
    }

    DOUBLE det = H00 * H11 - H01 * H01;
    if (fabs(det) < (DOUBLE)1e-10) return;

    DOUBLE c0_d = -(H11 * grad0 - H01 * grad1) / det;
    DOUBLE c1_d = -(H00 * grad1 - H01 * grad0) / det;
    VeryLong c0_vl(c0_d);
    VeryLong c1_vl(c1_d);
    c0 = c0_vl;
    c1 = c1_vl;
}

template <typename DOUBLE>
DOUBLE I(const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE b, DOUBLE& c0, DOUBLE& c1, DOUBLE t, DOUBLE s)
{
    best_c0_and_c1(f, a, b, t, s, c0, c1);
    return J(f, a, b, c0, c1, t, s);
}

template <typename DOUBLE>
DOUBLE minimize_I_over_t(const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE b, DOUBLE s)
{
    DOUBLE t = 0.0;
    DOUBLE last_t = 0.0;
    DOUBLE t_min = -MAX_SAMPLE_RANGE;//-1000000000;
    DOUBLE t_max = +MAX_SAMPLE_RANGE;//+1000000000;
    DOUBLE t_sample_size = SAMPLE_SIZE;
    DOUBLE t_delta = (t_max - t_min) / t_sample_size;
    DOUBLE t_try = t_min;
    DOUBLE t_diff = 0.0;
    DOUBLE value;
    DOUBLE c0 = 0.0;
    DOUBLE c1 = 0.0;
    DOUBLE min_value = J(f, a, b, c0, c1, t, s);
    int done = 0;
    //int iterations = 0;

    while (!done)
    {
        t_try = t_min;
        for (int l = 0; l < t_sample_size; l++)
        {
            value = J(f, a, b, c0, c1, t_try, s);
            if (value < min_value)
            {
                min_value = value;
                last_t = t;
                t = t_try;
            }
            t_try += t_delta;
        }
        t_min = t - t_sample_size*t_delta / 3;
        t_max = t + t_sample_size*t_delta / 3;
        t_delta = (t_max - t_min) / t_sample_size;
        t_diff = last_t - t;
        if (t_diff < 0) t_diff = -t_diff;
        last_t = t;

        const DOUBLE epsilon = 1.0;
        //iterations++;
        if (t_diff < epsilon && t_delta < 1.0) done = 1;
    }
    DOUBLE value1;
    VeryLong t_vl(t);
    if (t > 0)
    {
        value = J<DOUBLE>(f, a, b, c0, c1, t_vl, s);
        value1 = J<DOUBLE>(f, a, b, c0, c1, t_vl + VeryLong(1.0), s);
        if (value < value1) return t_vl;
        else return (t_vl + VeryLong(1.0));
    }
    else
    {
        value = J<DOUBLE>(f, a, b, c0, c1, t_vl, s);
        value1 = J<DOUBLE>(f, a, b, c0, c1, -(-t_vl + VeryLong(1.0)), s);
        if (value < value1) return t_vl;
        else return -(-t_vl + VeryLong(1.0));
    }
}

template <typename DOUBLE>
DOUBLE minimize_I_over_t(const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE b, DOUBLE& c0, DOUBLE& c1, DOUBLE s)
{
    DOUBLE t = 0.0;
    DOUBLE last_t = 0.0;
    DOUBLE t_min = -MAX_SAMPLE_RANGE;//-1000000000;
    DOUBLE t_max = +MAX_SAMPLE_RANGE;//+1000000000;
    DOUBLE t_sample_size = SAMPLE_SIZE;
    DOUBLE t_delta = (t_max - t_min) / t_sample_size;
    DOUBLE t_try = t_min;
    DOUBLE t_diff = 0.0;
    DOUBLE value;
    DOUBLE min_value = I(f, a, b, c0, c1, t, s);
    //std::cout << "A. min_value = " << min_value << std::endl;
    int done = 0;
    //int iterations = 0;

    while (!done)
    {
        t_try = t_min;
        for (int l = 0; l < t_sample_size; l++)
        {
            //std::cout << "(f, a, b, c0, c1, t_try, s) = (" << f << ", " << a << ", " << b << ", " << c0 << ", " << c1 << ", " << t_try << ", " << s << ")" << std::endl;
            value = I(f, a, b, c0, c1, t_try, s);
            //std::cout << "B. value = " << value << std::endl;
            if (value < min_value)
            {
                min_value = value;
                last_t = t;
                t = t_try;
            }
            t_try += t_delta;
        }
        t_min = t - t_sample_size*t_delta / 3.0;
        t_max = t + t_sample_size*t_delta / 3.0;
        t_delta = (t_max - t_min) / t_sample_size;
        t_diff = last_t - t;
        if (t_diff < 0) t_diff = -t_diff;
        last_t = t;

        const DOUBLE epsilon = 1.0;
        //iterations++;
        if (t_diff < epsilon && t_delta < (DOUBLE)1.0) done = 1;
    }
    DOUBLE value1;
    VeryLong t_vl(t);
    if (t > (DOUBLE)0.0)
    {
        value = J<DOUBLE>(f, a, b, c0, c1, t_vl, s);
        //std::cout << "C. t_vl = " << t_vl << std::endl;
        //std::cout << "D. value = " << value << std::endl;
        value1 = J<DOUBLE>(f, a, b, c0, c1, t_vl + VeryLong(1.0), s);
        //std::cout << "E. value1 = " << value1 << std::endl;
        if (value < value1) return t_vl;
        else return (t_vl + VeryLong(1.0));
    }
    else
    {
        value = J<DOUBLE>(f, a, b, c0, c1, t_vl, s);
        //std::cout << "F. t_vl = " << t_vl << std::endl;
        //std::cout << "G. value = " << value << std::endl;
        value1 = J<DOUBLE>(f, a, b, c0, c1, -(-t_vl + VeryLong(1.0)), s);
        //std::cout << "H. value1 = " << value1 << std::endl;
        if (value < value1) return t_vl;
        else return -(-t_vl + VeryLong(1.0));
    }
}

template <typename DOUBLE>
DOUBLE minimize_I_over_t_and_s(const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE b, DOUBLE& c0, DOUBLE& c1, DOUBLE& t, DOUBLE& s)
{
    //std::cout << "minimize_I_over_t_and_s(" << a << "," << b << ")" << std::endl;
    const int sample_size = SAMPLE_SIZE;
    const DOUBLE r = 10.0; // length of vector
    static DOUBLE s_vec[sample_size];
    static DOUBLE t_vec[sample_size];
    static int first_time = 1;

    if (first_time)
    {
        first_time = 0;
        DOUBLE theta_inc = 2.0 * M_PI / sample_size;
        DOUBLE theta = 0.0;
        for (int i = 0; i < sample_size; i++)
        {
            s_vec[i] = cos (theta);
            t_vec[i] = sin (theta);
            theta += theta_inc;
        }
    }

    DOUBLE min_value = J(f, a, b, c0, c1, t, s);
    //DOUBLE initial_min_value = min_value;
    DOUBLE min_c0 = c0;
    DOUBLE min_c1 = c1;
    DOUBLE min_s = s;
    DOUBLE min_t = t;
    //int min_i = -1;
    int min_i = 0;
    DOUBLE s_try;
    DOUBLE t_try;
    for (int i = 0; i < sample_size; i++)
    {
        s_try = s + r * s_vec[i];
        if (s_try < 1) s_try = 1;
        // s_try = 1;
        t_try = t + r * t_vec[i];
        DOUBLE c0_try;
        DOUBLE c1_try;
        DOUBLE value = I(f, a, b, c0_try, c1_try, t_try, s_try);
        if (value < min_value)
        {
            min_value = value;
            min_i = i;
        }
    }

    //if (min_value < initial_min_value)
    {
        // use sampling to get minimum point in this direction
        int done = 0;
        DOUBLE r_min = 0;
        DOUBLE r_max = 1e15;//MAX_SAMPLE_RANGE;
        DOUBLE r_sample_size = SAMPLE_SIZE;
        DOUBLE r_delta = (r_max - r_min) / r_sample_size;
        DOUBLE r_try = r_min;
        DOUBLE r_diff = 0.0;
        DOUBLE rr = 0;
        DOUBLE value;
        DOUBLE c0_try;
        DOUBLE c1_try;
        DOUBLE last_r = 0.0;
        while (!done)
        {
            r_try = r_min;
            t_try = t + r_try * t_vec[min_i];
            s_try = s + r_try * s_vec[min_i];
            if (s_try < 1) s_try = 1;
            // s_try = 1;
            for (int l = 0; l < r_sample_size; l++)
            {
                value = I(f, a, b, c0_try, c1_try, t_try, s_try);
                if (value < min_value)
                {
                    min_value = value;
                    last_r = rr;
                    rr = r_try;
                    min_s = s_try;
                    min_t = t_try;
                    min_c0 = c0_try;
                    min_c1 = c1_try;
                }
                r_try += r_delta;
                t_try += r_delta * t_vec[min_i];
                s_try += r_delta * s_vec[min_i];
                if (s_try < 1) s_try = 1;
                // s_try = 1;
            }
            r_min = rr - r_sample_size*r_delta / 3.0;
            if (r_min < 0) r_min = 0;
            r_max = rr + r_sample_size*r_delta / 3.0;
            r_delta = (r_max - r_min) / r_sample_size;
            r_diff = last_r - rr;
            if (r_diff < 0) r_diff = -r_diff;
            last_r = rr;

            const DOUBLE epsilon = 1.0;
            if (r_diff < epsilon && r_delta < (DOUBLE)1.0) done = 1;
        }

        // min point found
        c0 = min_c0;
        c1 = min_c1;
        s = min_s;
        t = min_t;
        //std::cout << "New minimum found at (c0,c1,t,s) = (" << c0 << "," << c1 << "," << t << "," << s << "), value = " << min_value << std::endl;
    }
    return min_value;
}

Polynomial<VeryLong> translate(const Polynomial<VeryLong>& p, const VeryLong& a, const VeryLong& b, const VeryLong& s, VeryLong& better_b, VeryLong& better_s);
double average_log_size(const Polynomial<VeryLong>& p, long int s);
double alpha_F(const Polynomial<VeryLong>& poly, long int B, long int BOUND);
double alpha_F_murphy(const Polynomial<VeryLong>& poly, long int B, long int BOUND);
struct Poly_info
{
    Poly_info(const Polynomial<VeryLong>& pp, const VeryLong& ss,
              double I, double alp) : p(pp), s(ss), I_F_S(I), alpha(alp)
    {}
    Polynomial<VeryLong> p;
    VeryLong s;
    double I_F_S;
    double alpha;
    bool operator<(const Poly_info& pi) const
    {
        return (I_F_S + alpha < pi.I_F_S + pi.alpha);
    }
};
Polynomial<VeryLong> adjust_root_properties(const Skewed_selection_config& Skewed_config,
        const Polynomial<VeryLong>& min_poly,
        const VeryLong& a,
        const VeryLong& b,
        VeryLong& s,
        double average_log_size,
        std::vector<Poly_info>& poly_list);
};
#endif
