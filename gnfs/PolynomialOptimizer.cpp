#include "PolynomialOptimizer.h"
#include "LongModular.h"
#include "discriminant.h"
#include <time.h>

namespace
{
    VeryLong b0_vl (const Polynomial<VeryLong>& f, const VeryLong& a, const VeryLong& b, const VeryLong& t)
    {
        VeryLong tt(1L);
        VeryLong res(0L);
        res = f.coefficient(0);
    
        for (int i = 1; i <= f.deg(); i++)
        {
            tt *= -1L * t;
            res += f.coefficient(i) * tt;
        }
        return res;
    }
    
    VeryLong b1_vl (const Polynomial<VeryLong>& f, const VeryLong& a, const VeryLong& b, const VeryLong& t)
    {
        VeryLong tt(1L);
        VeryLong res(0L);
        long int multiplier = 1;
        res = f.coefficient(1);
    
        for (int i = 2; i <= f.deg(); i++)
        {
            multiplier++;
            tt *= -1L * t;
            res += f.coefficient(i) * tt * multiplier;
        }
        return res;
    }
    
    VeryLong b2_vl (const Polynomial<VeryLong>& f, const VeryLong& a, const VeryLong& t)
    {
        VeryLong tt(1L);
        VeryLong res(0L);
        long int multiplier = 1;
        long int increment = 2;
        res = f.coefficient(2);
    
        for (int i = 3; i <= f.deg(); i++)
        {
            //cout << sign << " : " << multiplier << endl;
            multiplier += increment;
            increment++;
            tt *= -1L * t;
            res += f.coefficient(i) * tt * multiplier;
        }
        return res;
    }

    std::vector<VeryLong> F_sample;

    void sample_F_values(const Polynomial<VeryLong>& poly, long int BOUND)
    {
        F_sample.clear();
        const unsigned long int UPPER_LIMIT = 1000000L;
        unsigned long int a = 2;
        unsigned long int b = 2;
    
        for (long int i = 0; i < BOUND; i++)
        {
            for (long int j = 0; j < BOUND; j++)
            {
                a = genrand() % UPPER_LIMIT + 1;
                b = genrand() % UPPER_LIMIT + 1;
                while (gcd(a,b) != 1L)
                {
                    a = genrand() % UPPER_LIMIT + 1;
                    b = genrand() % UPPER_LIMIT + 1;
                }
                VeryLong value = poly.evaluate_homogeneous((long int)a,(long int)b);
                F_sample.push_back(value);
            }
        }
    }
    
    double cont_F(const Polynomial<VeryLong>& poly, long int p, long int BOUND, bool well_behaved, bool estimate)
    {
        double cont = 0.0;
        if (!estimate && well_behaved)
        {
            // use heuristic formula for cont_p(F)
            long int q_p = count_roots_mod_p<VeryLong, long int, LongModular>(poly, p);
            // check for projective root
            VeryLong c_d = poly.coefficient(poly.deg());
            long int p_k = 1L;
            while (c_d % p == 0L)
            {
                q_p++;
                c_d /= p;
                p_k *= p;
            }
            if (p_k > 1L)
            {
                // we know p^k divides c_d
                int i = poly.deg() - 1;
                while (i > 0)
                {
                    VeryLong c = poly.coefficient(i);
                    p_k /= p;
                    while (p_k > 1L && c % p_k != 0L)
                    {
                        p_k /= p;
                        q_p--;
                    }
                    i--;
                }
            }
    
            cont = q_p * p * 1.0 / (p * p - 1.0);
        }
        else
        {
            for (auto& i1: F_sample)
            {
                VeryLong value = i1;
                while (value % p == 0L)
                {
                    value /= p;
                    cont++;
                }
            }
            cont = cont / ((double)BOUND * (double)BOUND);
        }
        return cont;
    }
}

namespace PolynomialOptimizer
{
    VeryLong b0_vl (const Polynomial<VeryLong>& f, const VeryLong& a, const VeryLong& b,const VeryLong& c0, const VeryLong& c1, const VeryLong& t)
    {
        VeryLong res = ::b0_vl(f, a, b, t);
        res -= (b + a * t) * c0;
        res += (b + a * t) * c1 * t;
        return res;
    }
    
    // f(x_t) + (c1*(x - t) + c0)(a * (x - t) - b),
    VeryLong b1_vl (const Polynomial<VeryLong>& f, const VeryLong& a, const VeryLong& b, const VeryLong& c0, const VeryLong& c1, const VeryLong& t)
    {
        VeryLong res = ::b1_vl(f, a, b, t);
        res -= b * c1;
        res += c0 * a;
        res -= 2L * c1 * t * a;
        return res;
    }
    
    VeryLong b2_vl (const Polynomial<VeryLong>& f, const VeryLong& a, const VeryLong& c1, const VeryLong& t)
    {
        VeryLong res = ::b2_vl(f, a, t);
        res += c1 * a;
        return res;
    }
    
    VeryLong b3_vl (const Polynomial<VeryLong>& f, const VeryLong& t)
    {
        VeryLong tt(1L);
        VeryLong res(0L);
        long int multiplier = 1;
        long int increment = 3;
        long int increment1 = 3;
        res = f.coefficient(3);
    
        for (int i = 4; i <= f.deg(); i++)
        {
            //cout << sign << " : " << multiplier << endl;
            multiplier += increment;
            increment += increment1;
            increment1++;
            tt *= -1L * t;
            res += f.coefficient(i) * tt * multiplier;
        }
        return res;
    }
    
    VeryLong b4_vl (const Polynomial<VeryLong>& f, const VeryLong& t)
    {
        if (f.deg() < 4) return 0L;
        if (f.deg() < 5) return f.coefficient(4);
        VeryLong res = f.coefficient(4) - 5L * f.coefficient(5) * t;
        return res;
    }
    
    VeryLong b5_vl (const Polynomial<VeryLong>& f)
    {
        if (f.deg() < 5) return 0L;
        return f.coefficient(5);
    }

    Polynomial<VeryLong> translate(const Polynomial<VeryLong>& p, const VeryLong& a, const VeryLong& b, const VeryLong& s, VeryLong& better_b, VeryLong& better_s)
    {
        // find a better poly just by translation
        VeryLong better_t = minimize_I_over_t(Polynomial<VeryLong>::convert_to_double<double>(p), a.get_double(), b.get_double(), s.get_double());
        std::vector<VeryLong> bb;
        bb.resize(6);
        bb[0] = ::b0_vl(p,a,b,better_t);
        bb[1] = ::b1_vl(p,a,b,better_t);
        bb[2] = ::b2_vl(p,a,better_t);
        bb[3] = b3_vl(p,better_t);
        bb[4] = b4_vl(p,better_t);
        bb[5] = b5_vl(p);
        Polynomial<VeryLong> translated_poly(bb);
        better_b = b + better_t * a;
        better_s = minimize_I_over_s(Polynomial<VeryLong>::convert_to_double<double>(translated_poly), a.get_double(), better_b.get_double(), 0.0, 0.0, 0.0, 100000.0);
        return translated_poly;
    }

    double average_log_size(const Polynomial<VeryLong>& p, long int s)
    {
        return log(sqrt(PolynomialOptimizer::J(Polynomial<VeryLong>::convert_to_double<double>(p),(double)s)/4.0));
    }

    double alpha_F(const Polynomial<VeryLong>& poly, long int B, long int BOUND)
    {
        long int p = zpnextb(2);
        double alpha = 0.0;
        VeryLong disc = discriminant(poly);
        //std::cout << "disc(f) = " << disc << std::endl;
        double cont_p = 0.0;
        sample_F_values(poly, BOUND);
    
        while (p < B)
        {
            cont_p = cont_F(poly, p, BOUND, disc % p != 0, false);
            alpha += (1.0/(p - 1.0) - cont_p) * log((double)p);
            p = zpnext();
        }
        return alpha;
    }

    double alpha_F_murphy(const Polynomial<VeryLong>& poly, long int B, long int BOUND)
    {
        long int p = zpnextb(2);
        double alpha = 0.0;
        VeryLong disc = discriminant(poly);
        std::cout << "disc(f) = " << disc << std::endl;
        double cont_p = 0.0;
        sample_F_values(poly, BOUND);
    
        while (p < B)
        {
            // for small primes, use F-sample
            if (p < 100)
                cont_p = cont_F(poly, p, BOUND, disc % p != 0, true);
            else
                cont_p = cont_F(poly, p, BOUND, disc % p != 0, false);
            alpha += (1.0/(p - 1.0) - cont_p) * log((double)p);
            p = zpnext();
        }
        return alpha;
    }

    Polynomial<VeryLong> adjust_root_properties(const Skewed_selection_config& Skewed_config, 
                                                const Polynomial<VeryLong>& min_poly,
                                                const VeryLong& a,
                                                const VeryLong& b,
                                                VeryLong& s,
                                                double average_log_size,
                                                std::vector<Poly_info>& poly_list)
    {
        const VeryLong one(1L);
        std::cout << "Sieving for best root properties ..." << std::endl;
        std::cout << "min_poly = " << min_poly << std::endl;
        std::cout << "a = " << a << std::endl;
        std::cout << "b = " << b << std::endl;
        std::cout << "s = " << s << std::endl;
        std::cout << "average_log_size = " << average_log_size << std::endl;
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
        std::cout << "als = " << average_log_size << std::endl;
        std::cout << "PRINTING_BOUND = " << target_E_F << std::endl;
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
        long int k = 0;
    
        memset((char*)cont_array_data, 0, cont_array_data_size * sizeof(short));
        long int p = zpnextb(2);
        VeryLong leading_coefficient = min_poly.coefficient(min_poly.deg());
        time_t start = time(0);
        while (p < MAX_SMALL_PRIME)
        {
            p_k = p;
            k = 1;
            int projective_root = (leading_coefficient % p == 0L);
            double logp = log((double)p);
            while (p_k < MAX_P_K)
            {
                p_k *= p;
                k++;
            }
            p_k /= p;
            k--;
            //cout << "p^k = " << p_k << endl;
    
            double prob_here = (double)p / (p_k * (p + 1.0));
    
            {
                LongModular::set_default_modulus(p_k);
                //VeryLong mm = m % p_k;
                VeryLong aa = a % (long)p_k;
                VeryLong bb = b % (long)p_k;
                //VeryLong m_mod_p = m % p;
                //LongModular m_lm(mm.get_long());
                LongModular a_lm(aa.get_long());
                LongModular b_lm(bb.get_long());
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
                        //long int l_minus_m = (p_k + l - m_lm.get_long()) % p_k; // l - m (mod p^k)
                        long int al_minus_b = (p_k + a_lm.get_long() * l - b_lm.get_long()) % p_k; // a * l - b (mod p^k)
                        if (l < p_k)
                        {
                            // possible non-projective
                            //con0 = (f_value.get_long() + (p_k + j1) * l * l_minus_m) % p_k;
                            con0 = (f_value.get_long() + (p_k + j1) * l * al_minus_b) % p_k;
                            //con1 = l_minus_m;
                            con1 = al_minus_b;
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
                            //LongModular y = LongModular(p_k + p * (l - p_k));
                            LongModular y(static_cast<long int>(p_k + p * (l - p_k)));
                            LongModular con0_lm = f.evaluate_homogeneous(LongModular(1L), y);
                            //LongModular con1_lm = LongModular(1L) - m_lm * y;
                            LongModular con1_lm = LongModular(1L) - b_lm * y / a_lm;
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
                //p_k *= p;
                //k++;
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
        double best_s = s.get_double();
        double printing_bound = Skewed_config.PRINTING_BOUND() + 1.5;
        // calculate scaling for alpha estimate
    
        for (long int j1 = -MAX_J1; j1 <= MAX_J1; j1++)
        {
            for (long int j0 = -MAX_J0; j0 <= MAX_J0; j0++)
            {
                double cont = ((double)(cont_array [MAX_J1 + j1] [MAX_J0 + j0])*0.001);
                if (cont < best)
                {
                    best = cont;
                    j0_best = j0;
                    j1_best = j1;
                }
                if (j0 == 0 && j1 == 0)
                {
                    std::cout << "(j0,j1) = (" << j0 << "," << j1 << "), cont = " << cont << ", alpha_cutoff = " << alpha_cutoff << std::endl;
                }
                if (cont < alpha_cutoff)
                {
                    // worth a look
    
                    // std::cout << "(j0,j1) = (" << j0 << "," << j1 << "), cont = " << cont << ", alpha_cutoff = " << alpha_cutoff << std::endl;
                    std::vector<VeryLong> c;
                    c.resize(3);
                    c[0] = VeryLong(j0) * b;
                    c[1] = -one * (a * VeryLong(j0) + b * VeryLong(j1));
                    c[2] = VeryLong(j1) * a;
                    Polynomial<VeryLong> F = min_poly + Polynomial<VeryLong>(c);
    
                    //double I_F_S_0 = PolynomialOptimizer::average_log_size(F, 1000L);
                    s = PolynomialOptimizer::minimize_I_over_s(Polynomial<VeryLong>::convert_to_double<double>(F), s.get_double());
                    double I_F_S = PolynomialOptimizer::average_log_size(F, s);
                    if (I_F_S + cont < printing_bound)
                    {
                        poly_list.push_back(Poly_info(F, s, I_F_S, cont));
                    }
    
                    // since calculating alpha is fairly slow, only
                    // do it if we have a hope of beating the best E so far
                    if (I_F_S + cont < best_E)
                    {
                        double alpha = cont;
                        double E = I_F_S + alpha;
                        if (E < best_E)
                        {
                            best_alpha = alpha;
                            best_F = F;
                            best_E = E;
                            best_s = s.get_double();
                            std::cout << "j0 = " << j0 << " and j1 = " << j1 << ", cont = " << cont << std::endl;
                            //cout << F.evaluate(m) << endl;
                            std::cout << F.evaluate_homogeneous(b, a) << std::endl;
                            std::cout << "F = " << F << std::endl;
                            std::cout << "a = " << a << std::endl;
                            std::cout << "b = " << b << std::endl;
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
            //c[0] = j0_best * m;
            //c[1] = -1L * (j0_best + m * j1_best);
            c[0] = VeryLong(j0_best) * b;
            c[1] = -one * (a * VeryLong(j0_best) + b * VeryLong(j1_best));
            c[2] = VeryLong(j1_best) * a;
            Polynomial<VeryLong> F = min_poly + Polynomial<VeryLong>(c);
            //cout << F.evaluate(m) << endl;
            std::cout << F.evaluate_homogeneous(b, a) << std::endl;
            std::cout << "F = " << F << std::endl;
            double alpha = PolynomialOptimizer::alpha_F(F, 2000, 200);
            if (alpha < best_alpha)
            {
                best_alpha = alpha;
                best_F = F;
            }
            std::cout << "alpha = " << alpha << std::endl;
        }
    
        s = PolynomialOptimizer::minimize_I_over_s(Polynomial<VeryLong>::convert_to_double<double>(best_F), best_s);
        std::cout << "Best F = " << best_F << std::endl;
        std::cout << "a = " << a << std::endl;
        std::cout << "b = " << b << std::endl;
        std::cout << "alpha = " << best_alpha << std::endl;
        std::cout << "New s = " << s << std::endl;
    
        return best_F;
    }

}
