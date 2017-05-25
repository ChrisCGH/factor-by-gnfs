#ifndef __QUOTIENT_H
#define __QUOTIENT_H

#include <iostream>
#include <string>
#include "gcd.h"
#include "pow.h"

template <class I> class Quotient
{
public:
    Quotient() : a(0L), b(1L)
    {}
    Quotient(I i) : a(i), b(1L)
    {}
    Quotient(I i, I j) : a(i), b(j)
    {
        if (b.is_zero())
        {
            throw std::string("Quotient::Quotient() : divide by zero, quotient undefined");
        }
        else simplify();
    }
    Quotient(const Quotient& q) : a(q.a), b(q.b)
    {}
    int operator== (const Quotient& q) const
    {
        return (a * q.b == b * q.a);
    }
    bool is_zero() const
    {
        return (a.is_zero());
    }
    bool is_one() const
    {
        return (a == b);
    }
    bool operator!= (const Quotient& q) const
    {
        return !(*this == q);
    }

    Quotient& operator= (const Quotient& q)
    {
        if (this == &q) return *this;

        a = q.a;
        b = q.b;
        return *this;
    }

    I numerator() const
    {
        return a;
    }
    I denominator() const
    {
        const I zero(0L);
        const I one(1L);
        if (a != zero) return b;
        return one;
    }
    bool operator< (const Quotient& q) const
    {
        const I zero(0L);
        Quotient result = *this;
        result -= q;
        if ((result.a < zero && result.b > zero) ||
                (result.a > zero && result.b < zero))
        {
            return true;
        }
        return false;
    }

    bool operator<= (const Quotient& q) const
    {
        if (*this == q) return true;
        if (*this < q) return true;
        return false;
    }

    int operator> (const Quotient& q) const
    {
        return (!(*this <= q));
    }

    int operator>= (const Quotient& q) const
    {
        return (!(*this < q));
    }

    Quotient operator- () const
    {
        return Quotient(-a, b);
    }

    friend Quotient operator+ (const Quotient& q1, const Quotient& q2)
    {
        if (q1.is_zero()) return q2;
        if (q2.is_zero()) return q1;
        Quotient result;
        result.a = q1.a * q2.b + q1.b * q2.a;
        result.b = q1.b * q2.b;
        result.simplify();
        return result;
    }
    friend Quotient operator- (const Quotient& q1, const Quotient& q2)
    {
        if (q2.is_zero()) return q1;
        if (q1.is_zero()) return -q2;
        Quotient result;
        result.a = q1.a * q2.b - q1.b * q2.a;
        result.b = q1.b * q2.b;
        result.simplify();
        return result;
    }
    Quotient& operator-= (const Quotient& q2)
    {
        if (q2.is_zero()) return *this;
        if (this->is_zero())
        {
            *this = -q2;
            return *this;
        }

        a *= q2.b;
        a -= b * q2.a;
        b *= q2.b;
        simplify();
        return *this;
    }
    Quotient& operator+= (const Quotient& q2)
    {
        if (q2.is_zero()) return *this;
        if (is_zero())
        {
            *this = q2;
            return *this;
        }

        a *= q2.b;
        a += b * q2.a;
        b *= q2.b;
        simplify();
        return *this;
    }
#if 0
    Quotient& operator*=(long int li)
    {
        if (is_zero() || li == 1L) return *this;
        if (li == 0L || a == b)
        {
            a = li;
            b = 1L;
            return *this;
        }

        if (b % li == 0L)
        {
            b /= li;
        }
        else
        {
            a *= li;
        }
        return *this;
    }
#endif

    Quotient& operator*=(const Quotient& q)
    {
        if (is_zero() || q.a == q.b) return *this;
        if (q.is_zero() || a == b)
        {
            *this = q;
            return *this;
        }

        a *= q.a;
        b *= q.b;
        if (this != &q) simplify();
        return *this;
    }

    Quotient& operator*=(const I& i)
    {
        if (is_zero() || i == I(1L)) return *this;
        if (i == I(0L) || a == b)
        {
            *this = i;
            return *this;
        }

        a *= i;
        simplify();
        return *this;
    }

    friend Quotient operator* (const Quotient& q1, const Quotient& q2)
    {
        if (q1.is_zero()) return q1;
        if (q2.is_zero()) return q2;
        if (q1.a == q1.b) return q2;
        if (q2.a == q2.b) return q1;

        static Quotient result;
        result = q1;
        result.a *= q2.a;
        result.b *= q2.b;
        result.simplify();
        return result;
    }

    friend Quotient operator/ (const Quotient& q1, const Quotient& q2)
    {
        if (q2.is_zero())
        {
            throw std::string("Quotient::operator/() : divide by zero, quotient undefined");
        }

        if (q2.a == q2.b) return q1;
        Quotient result;
        result.a = q1.a;
        result.a *= q2.b;
        result.b = q1.b;
        result.b *= q2.a;
        result.simplify();
        return result;
    }

    friend std::ostream& operator<< (std::ostream& os, const Quotient& q)
    {
        os << q.a;
        if (q.b != 1L && q.a != 0L) os << "/" << q.b;

        return os;
    }

    /*
          friend std::istream& operator>> (std::istream& is, Quotient& q)
          {
             I a;
             I b;
             char temp;
             is >> a;
             is >> temp;
             is >> b;

             q = Quotient<I>(a,b);
             return is;
          }
    */

    void simplify()
    {
        const I zero(0L);
        const I one(1L);
        if (b < zero)
        {
            a = -a;
            b = -b;
        }
        if (b.is_one()) return;
        if (a.is_zero()) return;
        static I c;
        c = ::gcd(a,b);
        if (c.is_one()) return;
        a /= c;
        b /= c;
    }

    friend Quotient exp(const Quotient& q1, long int li)
    {
        Quotient res;
        res.a = pow(I(q1.a), I(li));
        res.b = pow(I(q1.b), I(li));
        return res;
    }

    Quotient& operator/=(const Quotient& q)
    {
        const I one(1L);
        if (this == &q)
        {
            if (is_zero())
            {
                throw std::string("Quotient::operator/=() : divide by zero, quotient undefined");
            }
            a = one;
            b = one;
        }
        else
        {
            if (q.is_zero())
            {
                throw std::string("Quotient::operator/=() : divide by zero, quotient undefined");
            }
            else
            {
                if (q.a != q.b)
                {
                    a *= q.b;
                    b *= q.a;
                    simplify();
                }
            }
        }
        return *this;
    }

private:
    // represents a/b
    I a;
    I b;

};

template <>
inline void Quotient<long int>::simplify()
{
    if (b < 0L)
    {
        a = -a;
        b = -b;
    }
    if (b == 1L) return;
    if (a == 0L) return;
    long int c = ::gcd(a,b);
    if (c == 1L) return;
    a /= c;
    b /= c;
}

template <>
inline Quotient<long int>::Quotient(long int i, long int j) : a(i), b(j)
{
    if (!b)
    {
        throw std::string("Quotient::Quotient() : divide by zero, quotient undefined");
    }
    else simplify();
}

template <>
inline Quotient<long int>& Quotient<long int>::operator+= (const Quotient<long int>& q2)
{
    if (q2.a == 0L) return *this;
    if (a == 0L)
    {
        *this = q2;
        return *this;
    }

    a *= q2.b;
    a += b * q2.a;
    b *= q2.b;
    simplify();
    return *this;
}

template <>
inline bool Quotient<long int>::is_zero() const
{
    return (a == 0L);
}
#if 0
template <>
Quotient<long int>& Quotient<long int>::operator*=(long int li)
{
}
#endif
#endif



