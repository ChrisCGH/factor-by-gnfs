#ifndef _MPFLOAT_H
#define _MPFLOAT_H

#include <gmp.h>
#include <iostream>
#include <string.h>

class MPFloat
{
public:
    friend class VeryLong;
    MPFloat()
    {
        mpf_init(d_);
    }

    MPFloat(const MPFloat& mpf)
    {
        mpf_init(d_);
        mpf_set(d_, mpf.d_);
    }

    MPFloat(double d)
    {
        mpf_init(d_);
        mpf_set_d(d_, d);
    }

    MPFloat(long int li)
    {
        mpf_init(d_);
        mpf_set_si(d_, li);
    }

    MPFloat(int i)
    {
        mpf_init(d_);
        mpf_set_si(d_, i);
    }

    ~MPFloat()
    {
        mpf_clear(d_);
    }

    MPFloat& operator=(const MPFloat& mpf)
    {
        if (this == &mpf) return *this;
        mpf_set(d_, mpf.d_);
        return *this;
    }

    int operator==(const MPFloat& mpf) const
    {
        return (mpf_cmp(d_, mpf.d_) == 0);
    }

    int operator==(long int li) const
    {
        return (mpf_cmp_si(d_, li) == 0);
    }

    int operator==(int i) const
    {
        return (mpf_cmp_si(d_, i) == 0);
    }

    int operator!=(const MPFloat& mpf) const
    {
        return (mpf_cmp(d_, mpf.d_) != 0);
    }

    int operator!=(long int li) const
    {
        return (mpf_cmp_si(d_, li) != 0);
    }

    int operator<(const MPFloat& mpf) const
    {
        return (mpf_cmp(d_, mpf.d_) < 0);
    }

    int operator<(long int li) const
    {
        return (mpf_cmp_si(d_, li) < 0);
    }

    int operator<(int i) const
    {
        return (mpf_cmp_si(d_, i) < 0);
    }

    int operator<=(const MPFloat& mpf) const
    {
        return (mpf_cmp(d_, mpf.d_) <= 0);
    }

    int operator<=(long int li) const
    {
        return (mpf_cmp_si(d_, li) <= 0);
    }

    int operator>(const MPFloat& mpf) const
    {
        return (mpf_cmp(d_, mpf.d_) > 0);
    }

    int operator>(long int li) const
    {
        return (mpf_cmp_si(d_, li) > 0);
    }

    int operator>=(const MPFloat& mpf) const
    {
        return (mpf_cmp(d_, mpf.d_) >= 0);
    }

    int operator>=(long int li) const
    {
        return (mpf_cmp_si(d_, li) >= 0);
    }

    friend MPFloat operator+(const MPFloat& mpf1, const MPFloat& mpf2)
    {
        MPFloat mpf;
        mpf_add(mpf.d_, mpf1.d_, mpf2.d_);
        return mpf;
    }

    friend MPFloat operator+(const MPFloat& mpf1, const double& d2)
    {
        MPFloat mpf;
        mpf_add(mpf.d_, mpf1.d_, MPFloat(d2).d_);
        return mpf;
    }

    friend MPFloat operator+(const double& d1, const MPFloat& mpf2)
    {
        return (mpf2 + d1);
    }

    friend MPFloat operator-(const MPFloat& mpf1, const MPFloat& mpf2)
    {
        MPFloat mpf;
        mpf_sub(mpf.d_, mpf1.d_, mpf2.d_);
        return mpf;
    }

    friend MPFloat operator-(const MPFloat& mpf1, const double& d2)
    {
        MPFloat mpf;
        mpf_sub(mpf.d_, mpf1.d_, MPFloat(d2).d_);
        return mpf;
    }

    friend MPFloat operator-(const double& d1, const MPFloat& mpf2)
    {
        MPFloat mpf;
        mpf_sub(mpf.d_, MPFloat(d1).d_, mpf2.d_);
        return mpf;
    }

    friend MPFloat operator*(const MPFloat& mpf1, const MPFloat& mpf2)
    {
        MPFloat mpf;
        mpf_mul(mpf.d_, mpf1.d_, mpf2.d_);
        return mpf;
    }

    friend MPFloat operator*(const MPFloat& mpf1, const double& d2)
    {
        MPFloat mpf;
        mpf_mul(mpf.d_, mpf1.d_, MPFloat(d2).d_);
        return mpf;
    }

    friend MPFloat operator*(const MPFloat& mpf1, const long int& li2)
    {
        MPFloat mpf;
        mpf_mul_ui(mpf.d_, mpf1.d_, li2);
        return mpf;
    }

    friend MPFloat operator*(const double& d1, const MPFloat& mpf2)
    {
        MPFloat mpf;
        mpf_mul(mpf.d_, MPFloat(d1).d_, mpf2.d_);
        return mpf;
    }

    friend MPFloat operator/(const MPFloat& mpf1, const MPFloat& mpf2)
    {
        MPFloat mpf;
        mpf_div(mpf.d_, mpf1.d_, mpf2.d_);
        return mpf;
    }

    friend MPFloat operator/(const MPFloat& mpf1, const double& d2)
    {
        MPFloat mpf;
        mpf_div(mpf.d_, mpf1.d_, MPFloat(d2).d_);
        return mpf;
    }

    friend MPFloat operator/(const double& d1, const MPFloat& mpf2)
    {
        MPFloat mpf;
        mpf_div(mpf.d_, MPFloat(d1).d_, mpf2.d_);
        return mpf;
    }

    MPFloat& operator+=(const MPFloat& mpf)
    {
        mpf_add(d_, d_, mpf.d_);
        return *this;
    }

    MPFloat& operator-=(const MPFloat& mpf)
    {
        mpf_sub(d_, d_, mpf.d_);
        return *this;
    }

    MPFloat& operator*=(const MPFloat& mpf)
    {
        mpf_mul(d_, d_, mpf.d_);
        return *this;
    }

    MPFloat& operator/=(const MPFloat& mpf)
    {
        mpf_div(d_, d_, mpf.d_);
        return *this;
    }

    MPFloat& operator+=(const double& d)
    {
        mpf_add(d_, d_, MPFloat(d).d_);
        return *this;
    }

    MPFloat& operator-=(const double& d)
    {
        mpf_sub(d_, d_, MPFloat(d).d_);
        return *this;
    }

    MPFloat& operator*=(const double& d)
    {
        mpf_mul(d_, d_, MPFloat(d).d_);
        return *this;
    }

    MPFloat& operator/=(const double& d)
    {
        mpf_div(d_, d_, MPFloat(d).d_);
        return *this;
    }

    operator double() const
    {
        return mpf_get_d(d_);
    }

    MPFloat operator-() const
    {
        MPFloat mpf(0.0);
        MPFloat mpf1;
        mpf_sub(mpf1.d_, mpf.d_, d_);
        return mpf1;
    }

    friend MPFloat fabs(const MPFloat& mpf)
    {
        if (mpf < MPFloat(0.0))
        {
            return -mpf;
        }
        return mpf;
    }

    static void set_precision(unsigned long int PREC)
    {
        mpf_set_default_prec(PREC);
    }

    friend std::ostream& operator<<(std::ostream& os, MPFloat& mpf)
    {
#if 0
        streambuf* sb = os.rdbuf();
        _IO_FILE* iof = static_cast<_IO_FILE*>(sb);
        FILE* f = (FILE*)iof;
        mpf_out_str(f, 10, 0, mpf.d_);
#endif
#if 1
        char tmp[1024];
        char tmp1[1024];
        mp_exp_t EXP;
        mpf_get_str(tmp, &EXP, 10, 0, mpf.d_);
        if (tmp[0] == '-')
        {
            tmp1[0] = tmp[0];
            tmp1[1] = tmp[1];
            tmp1[2] = '.';
            tmp1[3] = '\0';
            strcat(tmp1, tmp + 2);
        }
        else
        {
            tmp1[0] = tmp[0];
            tmp1[1] = '.';
            tmp1[2] = '\0';
            strcat(tmp1, tmp + 1);
        }
        os << tmp1 << "e" << EXP;
#endif
        return os;
    }

    static MPFloat random()
    {
        static int first_time = 1;
        static gmp_randstate_t rand_state;
        if (first_time)
        {
            first_time = 0;
            gmp_randinit (rand_state, GMP_RAND_ALG_LC, 128);
        }
        MPFloat mpf;
        mpf_urandomb(mpf.d_, rand_state, 100);
#if 1
        mpf *= MPFloat(3.14159e5);
        mpf -= MPFloat(3.14159e5 / 2.0);
#endif
        return mpf;
    }

private:
    mpf_t d_;
};
#endif
