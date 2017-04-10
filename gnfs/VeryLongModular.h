#ifndef __VERYLONGMODULAR_H
#define __VERYLONGMODULAR_H
#include <iostream>
#include "VeryLong.h"
extern "C"
{
#include <gmp.h>
}

class VeryLongModular
{
   public:
      friend class VeryLong;
      VeryLongModular();
      VeryLongModular(const VeryLong& modulus);
      VeryLongModular(const VeryLongModular& vl);
      VeryLongModular(const VeryLong& modulus, char* s);
      VeryLongModular(const VeryLong& modulus, long l);
      VeryLongModular(const VeryLong& modulus, const VeryLong& vl);
      explicit VeryLongModular(const long int li);
      VeryLongModular(const long long int lli);
      ~VeryLongModular();
      int operator== (const VeryLongModular& vl) const;
      int operator== (const VeryLong& vl) const;
      int operator== (const long int li) const
      {
         return (*this == VeryLong(li));
      }
      int is_zero() const;
      int is_one() const;
      int operator!= (const VeryLongModular& vl) const;
      /*  int operator!= (const VeryLong& vl) const
          {
             return !(*this == vl);
          }*/
      int operator< (const VeryLongModular& vl) const
      {
         return (mpz_cmp(vl_, vl.vl_) < 0);
      }
      int operator> (const VeryLongModular& vl) const
      {
         return (mpz_cmp(vl_, vl.vl_) > 0);
      }
      int operator< (const long int li) const
      {
         return 0;
      }
      VeryLongModular& operator= (const VeryLongModular& vl);
      VeryLongModular& operator= (const VeryLong& vl);
      VeryLongModular& operator+= (const VeryLongModular& vl);
      VeryLongModular& operator-= (const VeryLongModular& vl);
      VeryLongModular& operator*= (const VeryLongModular& vl);
      VeryLongModular& operator/= (const VeryLongModular& vl);
      VeryLongModular operator- () const;
      friend VeryLongModular inverse (const VeryLongModular& vl);
      VeryLongModular& exp(const long int& li);

      friend VeryLongModular exp(const VeryLongModular& vl1, const VeryLongModular& vl2);
      friend VeryLongModular exp(const VeryLongModular& vl1, const long int& li2);
      friend VeryLongModular operator+ (const VeryLongModular& vl1, const VeryLongModular& vl2);
      friend VeryLongModular operator- (const VeryLongModular& vl1, const VeryLongModular& vl2);
      friend VeryLongModular operator* (const VeryLongModular& vl1, const VeryLongModular& vl2);
      friend VeryLongModular operator/ (const VeryLongModular& vl1, const VeryLongModular& vl2);
      friend VeryLongModular operator% (const VeryLongModular& vl1, const VeryLongModular& vl2);

      friend ostream& operator<< (ostream& os, const VeryLongModular& vl);

      VeryLongModular& add_product(const VeryLongModular& x, const VeryLongModular& y, const VeryLongModular& z);

      VeryLong get_very_long() const
      {
         VeryLong vl;
         mpz_set(vl.vl_, vl_);
         return vl;
      }

      static void set_default_modulus(const VeryLong& modulus)
      {
         staticInit();
         mpz_set(modulus_, modulus.vl_);
      }
      static VeryLong get_default_modulus()
      {
         VeryLong vl;
         mpz_set(vl.vl_, modulus_);
         return vl;
      }
   private:
      static mpz_t modulus_;
      mpz_t vl_;
      static int staticInitDone;
      static void staticInit()
      {
         if (!staticInitDone)
         {
            mpz_init(modulus_);
            mpz_set_ui(modulus_, 1L);
            staticInitDone = 1;
         }
      }
};
#endif
