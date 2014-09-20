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

    template <typename DOUBLE> DOUBLE minimize_I_over_t(const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE b, DOUBLE s);

    template <typename DOUBLE> DOUBLE minimize_I_over_t(const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE b, DOUBLE& c0, DOUBLE& c1, DOUBLE s);
    template <typename DOUBLE> DOUBLE I(const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE b, DOUBLE& c0, DOUBLE& c1, DOUBLE t, DOUBLE s);

    template <typename DOUBLE> DOUBLE J(const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE bb, DOUBLE c0, DOUBLE c1, DOUBLE t, DOUBLE s);

    template <typename DOUBLE> DOUBLE minimize_I_over_t_and_s(const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE b, DOUBLE& c0, DOUBLE& c1, DOUBLE& t, DOUBLE& s);

    template <typename DOUBLE>
    DOUBLE J(const Polynomial<DOUBLE>& f, DOUBLE a, DOUBLE bb, DOUBLE c0, DOUBLE c1, DOUBLE t, DOUBLE s)
    {
        DOUBLE b[6];
        b[0] = b0(f,a,bb,c0,c1,t);
        b[1] = b1(f,a,bb,c0,c1,t);
        b[2] = b2(f,a,c1,t);
        b[3] = b3(f,t);
        b[4] = b4(f,t);
        b[5] = b5(f);
    
        DOUBLE s3 = s * s * s;
        DOUBLE s5 = s3 * s * s;
        DOUBLE rs = 1.0 / s;
        DOUBLE rs3 = 1.0 / s3;
        DOUBLE rs5 = 1.0 / s5;
    
        DOUBLE result = 0.0;
    
        result = (b[0] * b[0] * rs5 + b[5] * b[5] * s5) * 4.0 / 11.0;
        result += ((2.0 * b[0] * b[2] + b[1] * b[1]) * rs3 +
                   (2.0 * b[3] * b[5] + b[4] * b[4]) * s3) * 4.0 / 27.0;
        result += ((2.0 * b[0] * b[4] + 2.0 * b[1] * b[3] + b[2] * b[2]) * rs +
                   (2.0 * b[1] * b[5] + 2.0 * b[2] * b[4] + b[3] * b[3]) * s) * 4.0 / 35.0;
        return result;
    }

    template <typename DOUBLE>
    DOUBLE J(const Polynomial<DOUBLE>& f, DOUBLE s)
    {
        DOUBLE b[6];
        b[0] = f.coefficient(0);
        b[1] = f.coefficient(1);
        b[2] = f.coefficient(2);
        b[3] = f.coefficient(3);
        b[4] = b4(f);
        b[5] = b5(f);
    
        DOUBLE s3 = s * s * s;
        DOUBLE s5 = s3 * s * s;
        DOUBLE rs = 1.0 / s;
        DOUBLE rs3 = 1.0 / s3;
        DOUBLE rs5 = 1.0 / s5;
    
        DOUBLE result = 0.0;
    
        result = (b[0] * b[0] * rs5 + b[5] * b[5] * s5) * 4.0 / 11.0;
        result += ((2.0 * b[0] * b[2] + b[1] * b[1]) * rs3 +
                   (2.0 * b[3] * b[5] + b[4] * b[4]) * s3) * 4.0 / 27.0;
        result += ((2.0 * b[0] * b[4] + 2.0 * b[1] * b[3] + b[2] * b[2]) * rs +
                   (2.0 * b[1] * b[5] + 2.0 * b[2] * b[4] + b[3] * b[3]) * s) * 4.0 / 35.0;
        return result;
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
        int iterations = 0;
    
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
            iterations++;
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
        bb.resize(6);
        bb[0] = b0_vl(fm,a1,b1,c0_int,c1_int,t_vl);
        bb[1] = b1_vl(fm,a1,b1,c0_int,c1_int,t_vl);
        bb[2] = b2_vl(fm,a1,c1_int,t_vl);
        bb[3] = b3_vl(fm,t_vl);
        bb[4] = b4_vl(fm,t_vl);
        bb[5] = b5_vl(fm);
    
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
        if (f.deg() < 5) return f.coefficient(4);
        DOUBLE res = f.coefficient(4) - 5.0 * f.coefficient(5) * t;
        return res;
    }

    template <typename DOUBLE>
    DOUBLE b5 (const Polynomial<DOUBLE>& f)
    {
        if (f.deg() < 5) return 0.0;
        return f.coefficient(5);
    }

    //------------------------------------------
    // Functions for evaluating integral
    //------------------------------------------

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
    void best_c0_and_c1(const Polynomial<DOUBLE>& f, DOUBLE aa, DOUBLE bb, DOUBLE t, DOUBLE s, DOUBLE& c0, DOUBLE& c1, long int sample_size = 1)
    {
        DOUBLE det = a(aa, bb, t, s) * d(aa, bb, t, s);
        det -= b(aa, bb, t, s) * b(aa, bb, t, s);
        //std::cout << "best_c0_and_c1 : aa = " << aa << ", bb = " << bb << ", t = " << t << ", s = " << s << std::endl;
        //std::cout << "best_c0_and_c1 : det = " << det << std::endl;
        if (fabs(det) < (DOUBLE)1e-10) return;
        DOUBLE c0_d = (d(aa, bb, t, s) * e(f, aa, bb, t, s) - b(aa, bb, t, s) * ff(f, aa, bb, t, s)) / det;
        //std::cout << "best_c0_and_c1 : c0_d = " << c0_d << std::endl;
    
        DOUBLE c1_d = (a(aa, bb, t, s) * ff(f, aa, bb, t, s) - c(aa, bb, t, s) * e(f, aa, bb, t, s)) / det;
        //std::cout << "best_c0_and_c1 : c1_d = " << c1_d << std::endl;
        VeryLong c0_vl(c0_d);
        VeryLong c1_vl(c1_d);
        c0 = c0_vl;
        c1 = c1_vl;
        //std::cout << "best_c0_and_c1 : c0 = " << c0 << ", c1 = " << c1 << std::endl;
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
        int iterations = 0;
    
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
            iterations++;
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
        int iterations = 0;
    
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
            iterations++;
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
