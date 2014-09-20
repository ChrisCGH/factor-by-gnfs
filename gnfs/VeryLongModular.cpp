#include "VeryLongModular.h"

mpz_t VeryLongModular::modulus_;
int VeryLongModular::staticInitDone = 0;

VeryLongModular::VeryLongModular()
{
   staticInit();
   mpz_init(vl_);
}

VeryLongModular::~VeryLongModular()
{
   mpz_clear(vl_);
}

VeryLongModular::VeryLongModular(const VeryLong& vl)
{
   staticInit();
   mpz_init_set(vl_, vl.vl_);
   mpz_mod(vl_, vl_, modulus_);
}

VeryLongModular::VeryLongModular(const VeryLongModular& vl)
{
   mpz_init_set(vl_, vl.vl_);
}

VeryLongModular::VeryLongModular(const VeryLong& modulus, char* s)
{
   mpz_set(modulus_, modulus.vl_);
   mpz_init_set_str(vl_, s, 10);
   mpz_mod(vl_, vl_, modulus_);
}

VeryLongModular::VeryLongModular(const VeryLong& modulus, long l)
{
   mpz_set(modulus_, modulus.vl_);
   mpz_init_set_si(vl_, l);
   mpz_mod(vl_, vl_, modulus_);
}

VeryLongModular::VeryLongModular(const VeryLong& modulus, const VeryLong& vl)
{
   mpz_set(modulus_, modulus.vl_);
   mpz_init_set(vl_, vl.vl_);
   mpz_mod(vl_, vl_, modulus_);
}

VeryLongModular::VeryLongModular(const long int li)
{
   staticInit();
   mpz_init_set_si(vl_, li);
   mpz_mod(vl_, vl_, modulus_);
}

VeryLongModular::VeryLongModular(const long long int lli)
{
   staticInit();
   long long int num = lli;
   if (lli < 0) num = -num;
   unsigned long int q = (unsigned long int)(num >> 32);
   unsigned long int r = (unsigned long int)(num & 0x00000000ffffffff);
   mpz_init_set_ui(vl_, q);
   mpz_mul_2exp(vl_, vl_, 32);
   mpz_add_ui(vl_, vl_, r);
   mpz_mod(vl_, vl_, modulus_);
   if (lli < 0)
   {
      mpz_sub(vl_, modulus_, vl_);
   }
}

int VeryLongModular::operator== (const VeryLongModular& vl) const
{
   return (mpz_cmp(vl_, vl.vl_) == 0);
}

int VeryLongModular::is_zero() const
{
   return (mpz_sgn(vl_) == 0);
}

int VeryLongModular::is_one() const
{
   const VeryLongModular one(1L);
   return (*this == one);
}

int VeryLongModular::operator== (const VeryLong& vl) const
{
   static VeryLong tmp;
   mpz_mod(tmp.vl_, vl.vl_, modulus_);
   if (vl.is_zero())
   {
      return tmp.is_zero();
   }

   return (mpz_cmp(vl_, tmp.vl_) == 0);
}

int VeryLongModular::operator!= (const VeryLongModular& vl) const
{
   return (mpz_cmp(vl_, vl.vl_) != 0);
}

VeryLongModular& VeryLongModular::operator= (const VeryLongModular& vl)
{
   if (this == &vl) return *this;
   mpz_set(vl_, vl.vl_);
   return *this;
}

VeryLongModular& VeryLongModular::operator= (const VeryLong& vl)
{
   mpz_set(vl_, vl.vl_);
   mpz_mod(vl_, vl_, modulus_);
   return *this;
}

VeryLongModular operator+ (const VeryLongModular& vl1, const VeryLongModular& vl2)
{
   static VeryLongModular vl;
   mpz_add(vl.vl_, vl1.vl_, vl2.vl_);
   mpz_mod(vl.vl_, vl.vl_, VeryLongModular::modulus_);
   return vl;
}

VeryLongModular operator- (const VeryLongModular& vl1, const VeryLongModular& vl2)
{
   static VeryLongModular vl;
   mpz_sub(vl.vl_, vl1.vl_, vl2.vl_);
   mpz_mod(vl.vl_, vl.vl_, VeryLongModular::modulus_);
   return vl;
}

VeryLongModular operator* (const VeryLongModular& vl1, const VeryLongModular& vl2)
{
   static VeryLongModular vl;
   mpz_mul(vl.vl_, vl1.vl_, vl2.vl_);
   mpz_mod(vl.vl_, vl.vl_, VeryLongModular::modulus_);
   return vl;
}

VeryLongModular operator/ (const VeryLongModular& vl1, const VeryLongModular& vl2)
{
   if (vl2 == VeryLongModular(0L))
   {
      throw std::string("operator/(VeryLongModular) : divide by zero, quotient undefined");
   }
   static VeryLongModular vl;
   mpz_invert(vl.vl_, vl2.vl_, VeryLongModular::modulus_);
   mpz_mul(vl.vl_, vl.vl_, vl1.vl_);
   mpz_mod(vl.vl_, vl.vl_, VeryLongModular::modulus_);
   return vl;
}

VeryLongModular operator% (const VeryLongModular& vl1, const VeryLongModular& vl2)
{
   return VeryLongModular(0L);
}

VeryLongModular& VeryLongModular::operator+= (const VeryLongModular& vl)
{
   mpz_add(vl_, vl_, vl.vl_);
   mpz_mod(vl_, vl_, modulus_);
   return *this;
}

VeryLongModular& VeryLongModular::operator-= (const VeryLongModular& vl)
{
   mpz_sub(vl_, vl_, vl.vl_);
   mpz_mod(vl_, vl_, modulus_);
   return *this;
}

VeryLongModular& VeryLongModular::operator*= (const VeryLongModular& vl)
{
   mpz_mul(vl_, vl_, vl.vl_);
   mpz_mod(vl_, vl_, modulus_);
   return *this;
}

VeryLongModular& VeryLongModular::add_product(const VeryLongModular& x, const VeryLongModular& y, const VeryLongModular& z)
{
    VeryLongModular r = x;
    r *= y;
    r *= z;
    mpz_add(vl_, vl_, r.vl_);
    mpz_mod(vl_, vl_, modulus_);
    return *this;
}

VeryLongModular VeryLongModular::operator- () const
{
   VeryLongModular tmp(*this);
   mpz_sub(tmp.vl_, modulus_, vl_);
   mpz_mod(tmp.vl_, tmp.vl_, modulus_);
   return tmp;
}

VeryLongModular inverse(const VeryLongModular& vl)
{
   VeryLongModular tmp;
   mpz_invert(tmp.vl_, vl.vl_, VeryLongModular::modulus_);
   return tmp;
}

VeryLongModular& VeryLongModular::operator/= (const VeryLongModular& vl)
{
   if (vl == VeryLongModular(0L))
   {
      throw std::string("VeryLongModular::operator/= : divide by zero, quotient undefined");
   }
   VeryLong tmp;
   mpz_invert(tmp.vl_, vl.vl_, modulus_);
   mpz_mul(vl_, vl_, tmp.vl_);
   mpz_mod(vl_, vl_, modulus_);
   return *this;
}

ostream& operator<< (ostream& os, const VeryLongModular& vl)
{
   static char tmp[10240];
   mpz_get_str(tmp, 10, vl.vl_);
   os << tmp;
//   os << "(";
//   zswrite(str, vl._modulus);
//   os << str << ")";
   return os;
}

VeryLongModular exp(const VeryLongModular& vl1, const VeryLongModular& vl2)
{
   VeryLongModular vl;
   mpz_powm(vl.vl_, vl1.vl_, vl2.vl_, VeryLongModular::modulus_);
   return vl;
}

VeryLongModular exp(const VeryLongModular& vl1, const long int& li2)
{
   VeryLongModular vl;
   mpz_powm_ui(vl.vl_, vl1.vl_, li2, VeryLongModular::modulus_);
   return vl;
}

VeryLongModular& VeryLongModular::exp(const long int& li2)
{
   mpz_powm_ui(vl_, vl_, li2, modulus_);
   return *this;
}

