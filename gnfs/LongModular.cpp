#include "LongModular.h"
#include <sstream>
#include <cstdlib>
#include "legendre.h"

#if 1
inline unsigned long int myzmulmods(unsigned long int a, unsigned long int b, unsigned long int n)
{
    unsigned long int i =  (unsigned long int)( ((unsigned long long)a * (unsigned long long)b) % (unsigned long long)n);
//	long int i =  (long int)( ((__int64)a * (__int64)b) % (__int64)n);
    //if (i < 0) i += n;
    return i;
}
#endif
//#define myzmulmods zmulmods

unsigned long int LongModular::modulus_ = 0;
bool LongModular::staticInitDone = false;

LongModular::LongModular() : _li(0)
{
    staticInit();
}

LongModular::~LongModular()
{}

LongModular::LongModular(const long int li) : _li(li)
{
    staticInit();
    if (li < 0)
    {
        long int tmp = li;
        while (tmp < 0) tmp += modulus_;
        _li = tmp;
    }
    _li = _li % modulus_;
}

LongModular::LongModular(const long long int lli) : _li(0)
{
    staticInit();
    if (lli < 0)
    {
        long long int tmp = lli;
        while (tmp < 0) tmp += modulus_;
        _li = tmp % modulus_;
    }
    else
    {
        _li = lli % modulus_;
    }
}

LongModular::LongModular(const LongModular& vl) : _li(vl._li)
{}

LongModular::LongModular(const long int modulus, const long int l) : _li(l)
{
    modulus_ = modulus;
    if (l < 0)
    {
        long int tmp = l;
        while (tmp < 0) tmp += modulus_;
        _li = tmp;
    }
    //while (_li < 0) _li += modulus_;
    _li = _li % modulus_;
}

LongModular::LongModular(const std::string& s) : _li(0L)
{
    staticInit();
    long int li = std::atoi(s.c_str());
    if (li < 0)
    {
        long int tmp = li;
        while (tmp < 0) tmp += modulus_;
        _li = tmp;
    }
    else
    {
        _li = li;
    }
    _li = _li % modulus_;
}

int LongModular::operator== (const LongModular& vl) const
{
    return (_li == vl._li);
}

int LongModular::operator!= (const LongModular& vl) const
{
    return (_li != vl._li);
}

LongModular& LongModular::operator= (const LongModular& vl)
{
    if (this == &vl) return *this;
    _li = vl._li;
    return *this;
}

LongModular operator+ (const LongModular& vl1, const LongModular& vl2)
{
    LongModular vl;
    vl._li = (vl1._li + vl2._li) % LongModular::modulus_;
    return vl;
}

LongModular operator- (const LongModular& vl1, const LongModular& vl2)
{
    LongModular vl;
    if (vl1._li < vl2._li)
    {
        vl._li = LongModular::modulus_;
        vl._li -= vl2._li;
        vl._li += vl1._li;
    }
    else vl._li = vl1._li - vl2._li;
//   if (vl._li < 0) vl._li = LongModular::modulus_ + vl._li;
    vl._li %= LongModular::modulus_;
    return vl;
}

LongModular operator* (const LongModular& vl1, const LongModular& vl2)
{
    LongModular vl;
    vl._li = myzmulmods(vl1._li, vl2._li, LongModular::modulus_);
    return vl;
}

LongModular operator/ (const LongModular& vl1, const LongModular& vl2)
{
    if (vl2._li == 0L)
    {
        throw std::string("operator/(LongModular) : divide by zero, quotient undefined");
    }
    LongModular vl;
    vl._li = myzmulmods(vl1._li, zinvs(vl2._li, LongModular::modulus_), LongModular::modulus_);
    return vl;
}

LongModular& LongModular::operator+= (const LongModular& vl)
{
    _li = (_li + vl._li) % LongModular::modulus_;
    return *this;
}

LongModular& LongModular::operator-= (const LongModular& vl)
{
    if (_li < vl._li) _li = modulus_ + _li - vl._li;
    else _li = (_li - vl._li);
    //if (_li < 0) _li += LongModular::modulus_;
    _li %= LongModular::modulus_;
    return *this;
}

LongModular& LongModular::operator*= (const LongModular& vl)
{
#if 0
    _li = myzmulmods(_li, vl._li, LongModular::modulus_);
#else
    _li = (unsigned long int)( ((unsigned long long)_li * (unsigned long long)vl._li) % LongModular::modulus_);

#endif
    return *this;
}

LongModular& LongModular::add_product(const LongModular& x, const LongModular& y, const LongModular& z)
{
    unsigned long long int r = (static_cast<unsigned long long int>(x._li) * y._li) % modulus_;
    _li = (_li + r * z._li) % modulus_;
    return *this;
}

LongModular LongModular::operator- () const
{
    LongModular tmp(*this);
    tmp._li = LongModular::modulus_ - _li;
    return tmp;
}

LongModular& LongModular::operator/= (const LongModular& vl)
{
    if (vl._li == 0L)
    {
        throw std::string("LongModular::operator/= : divide by zero, quotient undefined");
    }
    _li = myzmulmods(_li, zinvs(vl._li, LongModular::modulus_), LongModular::modulus_);
    return *this;
}

ostream& operator<< (ostream& os, const LongModular& vl)
{
    os << vl._li;
    return os;
}

LongModular exp(const LongModular& vl1, const LongModular& vl2)
{
    LongModular vl;
    vl._li = zexpmods(vl1._li, vl2._li, LongModular::modulus_);
    return vl;
}

LongModular& LongModular::exp(const long int& li2)
{
    _li = zexpmods(_li, li2, LongModular::modulus_);
    return *this;
}

// Algorithm 1.5.1 (Square Root Mod p)
LongModular LongModular::square_root() const
{
    long int p = modulus_;
    // write p - 1 as 2^e . q with q odd
    long int q = p - 1;
    long int e = 0;
    while ((q & 1) == 0)
    {
        q >>= 1;
        ++e;
    }
    // 1. [Find generator]
    long int n = 1;
    while (legendre<long int>(n, p) != -1)
    {
        n++;
    }

    long int z = zexpmods(n, q, p);

    // 2. [Initialize]
    long int y = z;
    long int r = e;
    long int x = zexpmods(_li, (q - 1L) / 2L, p);
    long int b = myzmulmods(_li, x, p);
    b = myzmulmods(b, x, p);
    x = myzmulmods(_li, x, p);

    while (true)
    {
        // 3. [Find exponent]
        if (b == 1L)
        {
            return LongModular(x);
        }

        long int m = 1;
        long int b2 = myzmulmods(b, b, p);
        long int b2m = b2;
        while (b2m != 1L)
        {
            b2m = myzmulmods(b2m, b2, p);
            ++m;
        }

        if (m == r)
        {
            std::ostringstream oss;
            oss << _li << " is not a quadratic residue mod " << p;
            throw oss.str();
        }

        // 4. [Reduce exponent]
        if (r - m - 1 <= 0)
        {
            std::ostringstream oss;
            oss << "r - m - 1 = " << r - m - 1 << " <= 0";
            throw oss.str();
        }
        long int rm1 = r - m - 1;
        long int power = 1 << rm1;
        long int t = zexpmods(y, power, p);
        y = myzmulmods(t, t, p);
        r = m;
        x = myzmulmods(x, t, p);
        b = myzmulmods(b, y, p);
    }
}
