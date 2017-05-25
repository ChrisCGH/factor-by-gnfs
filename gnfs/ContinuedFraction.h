#ifndef _CONTINUED_FRACTION_H
#define _CONTINUED_FRACTION_H
#include "Quotient.h"
#include <vector>

template <class Z> class Continued_Fraction
{
public:
    Continued_Fraction(const double number)
    {
        double x = number;
        long int a = (long int)x;
        x = x - a;
        x = 1 / x;
        cf.push_back(Z(a));
        for (int i = 0; i < 1000; i++)
        {
            a = (long int)x;
            if (a == 0) break;
            x = x - a;
            x = 1/x;
            cf.push_back(Z(a));
        }
    }
    Continued_Fraction(const Z& numerator, const Z& denominator)
    {
        Quotient<Z> x(numerator, denominator);
        std::cout << x << std::endl;
        Z a = numerator / denominator;
        const Quotient<Z> one(1L);
        x = x - Quotient<Z>(a);
        x = one / x;
        cf.push_back(a);
        for (int i = 0; i < 1000; i++)
        {
            std::cout << "x = " << x << std::endl;
            std::cout << x.numerator() << std::endl;
            std::cout << x.denominator() << std::endl;
            a = x.numerator() / x.denominator();
            if (a == 0L) break;
            x = x - Quotient<Z>(a);
            if (x.numerator() == 0L) break;
            x = one / x;
            cf.push_back(a);
        }

    }
    ~Continued_Fraction()
    {}
    Quotient<Z> nearest_rational_approximation(int depth)
    {
        int i = depth;
        if (i >= (int)cf.size()) i = cf.size() - 1;
        if (i < 1) i = 1;
        const Z one(1L);
        const Quotient<Z> one_q(1L, 1L);
        Quotient<Z> q(one, cf[i]);
        i--;
        while (i > 0)
        {
            q = Quotient<Z>(cf[i]) + q;
            q = one_q / q;
            i--;
        }
        q = Quotient<Z>(cf[i]) + q;
        return q;
    }
    void display()
    {
        for (int i = 0; i < cf.size(); i++)
        {
            std::cout << i << ":" << cf[i] << std::endl;
        }
    }


private:
    std::vector<Z> cf;
};

#endif

