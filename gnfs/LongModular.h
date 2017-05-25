#ifndef __LONGMODULAR_H
#define __LONGMODULAR_H
#include <iostream>
using std::ostream;
#include "VeryLong.h"
extern "C"
{
#include "lip.h"
}

class LongModular
{
public:
    friend class VeryLong;
    LongModular();
    LongModular(const long int li);
    LongModular(const long long int li);
    LongModular(const LongModular& vl);
    LongModular(const long int modulus, const long int li);
    LongModular(const std::string& s);

    ~LongModular();
    int operator== (const LongModular& vl) const;
    int operator!= (const LongModular& vl) const;
    int operator< (const LongModular& vl) const
    {
        return 0;
    }
    int operator> (const LongModular& vl) const
    {
        return 0;
    }
    int operator< (const long int li) const
    {
        return 0;
    }
    LongModular& operator= (const LongModular& vl);
    LongModular& operator+= (const LongModular& vl);
    LongModular& operator-= (const LongModular& vl);
    LongModular& operator*= (const LongModular& vl);
    LongModular& operator/= (const LongModular& vl);
    LongModular operator- () const;
    //LongModular& inverse ();
    LongModular& exp(const long int& li);

    friend LongModular exp(const LongModular& vl1, const LongModular& vl2);
    //friend LongModular exp(const LongModular& vl1, const long int& li2);
    friend LongModular operator+ (const LongModular& vl1, const LongModular& vl2);
    friend LongModular operator- (const LongModular& vl1, const LongModular& vl2);
    friend LongModular operator* (const LongModular& vl1, const LongModular& vl2);
    friend LongModular operator/ (const LongModular& vl1, const LongModular& vl2);

    friend ostream& operator<< (ostream& os, const LongModular& vl);

    LongModular& add_product(const LongModular& x, const LongModular& y, const LongModular& z);

    unsigned long int get_long() const
    {
        return _li;
    }
    int is_zero() const
    {
        return (_li == 0L);
    }

    static void set_default_modulus(const unsigned long int modulus)
    {
        modulus_ = modulus;
        staticInitDone = true;
    }
    static unsigned long int get_default_modulus()
    {
        return modulus_;
    }
    LongModular square_root() const;
private:
    static unsigned long int modulus_;
    static void staticInit()
    {
        if (!staticInitDone)
        {
            modulus_ = 1L;
            staticInitDone = true;
        }
    }
    static bool staticInitDone;
    unsigned long int _li;
};
#endif
