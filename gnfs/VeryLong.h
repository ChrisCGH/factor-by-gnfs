#ifndef __VERYLONG_H
#define __VERYLONG_H
#include <iostream>
#include <stdexcept>
#include <stdint.h>
using std::ostream;
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#ifndef WIN32
#include <netinet/in.h>
#endif
#include <vector>
#include <gmp.h>
#include "mod.h"
#include "MPFloat.h"
#define USE_QS 1
class VeryLong
{
   public:
      friend class VeryLongModular;
      VeryLong()
      {
         mpz_init(vl_);
      }
      VeryLong(const VeryLong& vl)
      {
         mpz_init_set(vl_, vl.vl_);
      }
      VeryLong(const char* s)
      {
         mpz_init_set_str(vl_, s, 0);
      }
      VeryLong(const std::string& s)
      {
         std::string s1(s);
         mpz_init_set_str(vl_, s1.c_str(), 0);
      }
      VeryLong(long l)
      {
         mpz_init_set_si(vl_, l);
      }
      VeryLong(unsigned long l)
      {
         mpz_init_set_ui(vl_, l);
      }
      VeryLong(unsigned long long l)
      {
         uint32_t q = static_cast<uint32_t>(l >> 32);
         uint32_t r = static_cast<uint32_t>(l & 0x00000000ffffffff);
         mpz_init_set_ui(vl_, q);
         mpz_mul_2exp(vl_, vl_, 32);
         mpz_add_ui(vl_, vl_, r);
      }
      VeryLong(long long l)
      {
         int negative = 0;
         if (l < 0)
         {
            l = -l;
            negative = 1;
         }
         uint32_t q = static_cast<uint32_t>(l >> 32);
         uint32_t r = static_cast<uint32_t>(l & 0x00000000ffffffff);
         mpz_init_set_ui(vl_, q);
         mpz_mul_2exp(vl_, vl_, 32);
         mpz_add_ui(vl_, vl_, r);
         if (negative) mpz_neg(vl_, vl_);
      }
      VeryLong(double d)
      {
         mpz_init_set_d(vl_, d);
      }
      explicit VeryLong(long double d)
      {
         mpz_init_set_d(vl_, d);
      }
      ~VeryLong()
      {
         mpz_clear(vl_);
      }
      bool operator== (const VeryLong& vl) const
      {
         return (mpz_cmp(vl_, vl.vl_) == 0);
      }
      bool operator== (long int li) const
      {
         return (mpz_cmp_si(vl_, li) == 0);
      }
      bool operator!= (long int li) const
      {
         return (mpz_cmp_si(vl_, li) != 0);
      }
      bool is_zero() const
      {
         return (mpz_sgn(vl_) == 0);
      }
      bool is_one() const
      {
         const VeryLong one(1L);
         return (*this == one);
      }
      operator MPFloat() const
      {
          MPFloat mpf;
          mpf_set_z(mpf.d_, vl_);
          
          return mpf;
      }
      operator double() const
      {
          return mpz_get_d(vl_);
      }
      bool operator!= (const VeryLong& vl) const
      {
         return (mpz_cmp(vl_, vl.vl_) != 0);
      }
      bool operator< (const VeryLong& vl) const
      {
         return (mpz_cmp(vl_, vl.vl_) < 0);
      }
      bool operator< (long int li) const
      {
         return (mpz_cmp_si(vl_, li) < 0);
      }
      bool operator<= (const VeryLong& vl) const
      {
         return (mpz_cmp(vl_, vl.vl_) <= 0);
      }
      bool operator<= (long int li) const
      {
         return (mpz_cmp_ui(vl_, li) <= 0);
      }
      bool operator> (const VeryLong& vl) const
      {
         return (mpz_cmp(vl_, vl.vl_) > 0);
      }
      bool operator> (long int li) const
      {
         return (mpz_cmp_si(vl_, li) > 0);
      }
      bool operator>= (const VeryLong& vl) const
      {
         return (mpz_cmp(vl_, vl.vl_) >= 0);
      }
      bool operator>= (long int li) const
      {
         return (mpz_cmp_si(vl_, li) >= 0);
      }
      VeryLong& operator= (const VeryLong& vl)
      {
         if (this == &vl) return *this;
         mpz_set(vl_, vl.vl_);
         return *this;
      }
      VeryLong& operator+= (const VeryLong& vl)
      {
         mpz_add(vl_, vl_, vl.vl_);
         return *this;
      }
      VeryLong& operator-= (const VeryLong& vl)
      {
         mpz_sub(vl_, vl_, vl.vl_);
         return *this;
      }
      VeryLong& operator*= (const VeryLong& vl)
      {
         mpz_mul(vl_, vl_, vl.vl_);
         return *this;
      }
      VeryLong& operator/= (const VeryLong& vl)
      {
         mpz_fdiv_q(vl_, vl_, vl.vl_);
         return *this;
      }
      VeryLong& operator/= (long int li)
      {
         mpz_fdiv_q_ui(vl_, vl_, li);
         return *this;
      }
      VeryLong operator- () const
      {
         VeryLong vl;
         mpz_neg(vl.vl_, vl_);
         return vl;
      }
      void negate()
      {
         mpz_neg(vl_, vl_);
      }
      double get_double() const
      {
         return mpz_get_d(vl_);
      }
      long double get_long_double() const
      {
         return mpz_get_d(vl_);
      }
      friend VeryLong exp(const VeryLong& vl1, const VeryLong& vl2);
      friend double ln(const VeryLong& vl);
      friend double log10(const VeryLong& vl);
      friend VeryLong operator+ (const VeryLong& vl1, const VeryLong& vl2)
      {
         VeryLong res;
         mpz_add(res.vl_, vl1.vl_, vl2.vl_);
         return res;
      }
      friend VeryLong operator- (const VeryLong& vl1, const VeryLong& vl2)
      {
         VeryLong res;
         mpz_sub(res.vl_, vl1.vl_, vl2.vl_);
         return res;
      }
      friend VeryLong operator* (const VeryLong& vl1, const VeryLong& vl2)
      {
         VeryLong res;
         mpz_mul(res.vl_, vl1.vl_, vl2.vl_);
         return res;
      }
      friend VeryLong operator* (const VeryLong& vl1, long int li)
      {
         VeryLong res;
         mpz_mul_si(res.vl_, vl1.vl_, li);
         return res;
      }
      friend VeryLong operator* (long int li, const VeryLong& vl1)
      {
         VeryLong res;
         mpz_mul_si(res.vl_, vl1.vl_, li);
         return res;
      }
      friend VeryLong operator/ (const VeryLong& vl1, const VeryLong& vl2)
      {
         VeryLong res;
         mpz_fdiv_q(res.vl_, vl1.vl_, vl2.vl_);
         return res;
      }
      friend VeryLong divide(const VeryLong& vl1, const VeryLong& vl2)
      {
         VeryLong res;
         mpz_fdiv_q(res.vl_, vl1.vl_, vl2.vl_);
         return res;
      }
      friend VeryLong operator/ (const VeryLong& vl1, const long int& li2)
      {
         VeryLong res;
         mpz_fdiv_q_ui(res.vl_, vl1.vl_, li2);
         return res;
      }
      VeryLong& operator%= (const VeryLong& vl)
      {
         mpz_mod(vl_, vl_, vl.vl_);
         return *this;
      }
      friend VeryLong operator% (const VeryLong& vl1, const VeryLong& vl2)
      {
         VeryLong res;
         mpz_mod(res.vl_, vl1.vl_, vl2.vl_);
         return res;
      }
      friend long int operator% (const VeryLong& vl1, const long int& li2)
      {
         return mpz_fdiv_ui(vl1.vl_, li2);
      }
      friend int jacobi_symbol(const VeryLong& a, const VeryLong& n)
      {
         return mpz_jacobi(a.vl_, n.vl_);
      }
      VeryLong nth_root(const long int n) const
      {
         VeryLong res;
         mpz_root(res.vl_, vl_, n);
         return res;
      }
      friend VeryLong gcd(const VeryLong& a, const VeryLong& b)
      {
         VeryLong res;
         mpz_gcd(res.vl_, a.vl_, b.vl_);
         return res;
      }
      friend VeryLong extended_gcd(const VeryLong& a, const VeryLong& b, VeryLong& xa, VeryLong& xb)
      {
         VeryLong res;
         mpz_gcdext(res.vl_, xa.vl_, xb.vl_, a.vl_, b.vl_);
         return res;
      }
      VeryLong inverse(VeryLong p) const
      {
         // inverse mod p 
         VeryLong res;
         mpz_invert(res.vl_, vl_, p.vl_);
         return res;
      }
      friend void combine(VeryLong& x, const VeryLong& v, const VeryLong& q)
      {
         // replace x by x - q * v
         VeryLong tmp;
         mpz_mul(tmp.vl_, v.vl_, q.vl_);
         mpz_sub(x.vl_, x.vl_, tmp.vl_);
      }
      bool is_probable_prime(long int number_of_tests = 10) const
      {
         return !(mpz_probab_prime_p(vl_, number_of_tests) == 0);
      }
      bool is_prime_power() const
      {
         return !(mpz_perfect_power_p(vl_) == 0);
      }
      bool is_square() const
      {
         return (mpz_perfect_square_p(vl_) != 0);
      }
      bool factorise_pollardrho(VeryLong* factor, VeryLong* new_n, long int tries = 0);
      bool factorise_qs(VeryLong* factor, VeryLong* new_n);
      bool factorise_squfof(VeryLong* factor, VeryLong* new_n);
      bool factorise_p_minus_1(VeryLong* factor, VeryLong* new_N);
      bool factorise_trial_division(std::vector<VeryLong> * factors, VeryLong* new_N) const;
      bool factorise_trial_division(std::vector<long int> * factors, VeryLong* new_N) const;
      bool factorise_trial_division(std::vector<long int> * factors, VeryLong* new_N, long int B) const;
      bool factorise_trial_division(VeryLong* factor, VeryLong* new_N) const;
      bool factorise_ecm(VeryLong* factor, VeryLong* new_N);
      bool factorise_fermat(VeryLong* factor, VeryLong* new_N);
      bool factorise(std::vector<VeryLong> * factors) const;
      bool factorise_no_trial(std::vector<VeryLong> * factors) const;
      bool factorise_no_trial(std::vector<VeryLong> * factors, VeryLong& unfactored_part) const;
      bool factorise(std::vector<VeryLong>* factors, VeryLong& unfactored_part) const;
      static bool factorise_no_trial_ll(long long int lli, std::vector<long int>* factors);
      static bool factorise_no_trial_l(long int li, std::vector<long int> * factors);

      friend ostream& operator<< (ostream& os, const VeryLong& vl)
      {
         static char tmp[10240];
         mpz_get_str(tmp, 10, vl.vl_);
         os << tmp;
         return os;
      }
      long int get_long() const
      {
         //return mpz_get_ui(vl_);
         return mpz_get_si(vl_);
      }
      long long int get_long_long() const
      {
         mpz_t v;
         mpz_init(v);
         mpz_div_2exp(v, vl_, 32);
         long long ll = mpz_get_si(v);
         ll <<= 32;
         mpz_mod_2exp(v, vl_, 32);
         ll += mpz_get_ui(v);
         mpz_clear(v);
         return ll;
      }
      static VeryLong random(const VeryLong& a, const VeryLong& b);
      friend VeryLong abs(const VeryLong& vl)
      {
         VeryLong vl1;
         mpz_abs(vl1.vl_, vl.vl_);
         return vl1;
      }
      VeryLong abs() const
      {
         VeryLong vl;
         mpz_abs(vl.vl_, vl_);
         return vl;
      }
      static void addPrime(const VeryLong& p)
      {
         extraPrimes_.push_back(p);
      }
      static long int firstPrime()
      {
         PrimeTablePtr_ = PrimeTable_;
         CurrentPrime_ = 2L;
         return CurrentPrime_;
      }
      static long int nextPrime()
      {
         // 3 marks end of table (differences of primes is never = 3)
         if (*PrimeTablePtr_ != 3)
         {
            CurrentPrime_ += *PrimeTablePtr_;
            PrimeTablePtr_++;
         }
         return CurrentPrime_;
      }
      static void setDebug();
      static void clearDebug();

      static void generate_prime_table();
      static void set_max_prime(long int B);
      static void clear_prime_table();
      static long int get_max_prime()
      {
         return Max_prime_;
      }
      friend bool kleinjung(const VeryLong&, const VeryLong&, VeryLong&, VeryLong&, std::fstream&);
   protected:
      mpz_t vl_;
      static std::vector<long int> Prime_power;
      static std::vector<long int> Prime_power2;
      static long int Max_prime_;
      static int PrimeTableSize_;
      static unsigned char* PrimeTable_;
      static unsigned char* PrimeTablePtr_;
      static long int CurrentPrime_;
      static long int Max_diff;
      static void generate_prime_powers();
      static std::vector<VeryLong> extraPrimes_;
};
#define FASTVERYLONG 1
#ifdef FASTVERYLONG
class FastVeryLong : public VeryLong
{
   public:
      FastVeryLong() : VeryLong(), ll_(0), l_(0), fast_(false), faster_(false)
      {}
      FastVeryLong(const FastVeryLong& fvl) : VeryLong(fvl), ll_(fvl.ll_), l_(fvl.l_), fast_(fvl.fast_), faster_(fvl.faster_)
      {}
      FastVeryLong(const VeryLong& vl) : VeryLong(vl), ll_(0), l_(0), fast_(false), faster_(false)
      {
         check();
         check1();
      }
      FastVeryLong(char* s) : VeryLong(s), ll_(0), l_(0), fast_(false), faster_(false)
      {
         check();
         check1();
      }
      FastVeryLong(long l) : VeryLong(l), ll_(0), l_(0), fast_(false), faster_(false)
      {}
      FastVeryLong(unsigned long l) : VeryLong(l), ll_(0), l_(0), fast_(false), faster_(false)
      {}
      FastVeryLong(long long l) : VeryLong(l), ll_(0), l_(0), fast_(false), faster_(false)
      {}
      FastVeryLong(double d) : VeryLong(d), ll_(0), l_(0), fast_(false), faster_(false)
      {}
      explicit FastVeryLong(long double d) : VeryLong(d), ll_(0), l_(0), fast_(false), faster_(false)
      {}
      FastVeryLong& operator=(const FastVeryLong& fvl)
      {
         if (this == &fvl) return *this;
         ll_ = fvl.ll_;
         l_ = fvl.l_;
         fast_ = fvl.fast_;
         faster_ = fvl.faster_;
         mpz_set(vl_, fvl.vl_);
         return *this;
      }
      ~FastVeryLong()
      {}
      void display(std::ostream& os) const
      {
         if (faster_) os << l_;
         else if (fast_) os << ll_;
         else os << static_cast<const VeryLong&>(*this);
      }
      FastVeryLong& operator/=(long int li)
      {
         if (faster_)
         {
            l_ /= li;
         }
         else if (fast_)
         {
            ll_ /= li;
            check1();
         }
         else
         {
            static_cast<VeryLong&>(*this) /= li;
            check();
            check1();
         }
         return *this;
      }
      friend FastVeryLong operator/ (const FastVeryLong& fvl1, const long int& li2)
      {
         FastVeryLong res(fvl1);
         if (res.faster_)
         {
            res.l_ /= li2;
         }
         else if (res.fast_)
         {
            res.ll_ /= li2;
         }
         else
         {
            mpz_fdiv_q_ui(res.vl_, fvl1.vl_, li2);
         }
         return res;
      }

      bool divisible_by(const long int& li)
      {
         if (faster_)
         {
            if (l_ % li == 0)
            {
               l_ /= li;
               return true;
            }
            return false;
         }

         if (fast_)
         {
            modasm2(ll_, li, r);
            if (r == 0)
            {
               ll_ /= li;
               check1();
               return true;
            }
            return false;
         }

         if (static_cast<const VeryLong&>(*this) % li == 0L)
         {
            static_cast<VeryLong&>(*this) /= li;
            check();
            check1();
            return true;
         }
         return false;
      }

      friend long int operator% (const FastVeryLong& fvl, const long int& li)
      {
         if (fvl.faster_)
         {
            return static_cast<long int>(fvl.l_ % li);
         }
         if (fvl.fast_)
         {
            modasm2(fvl.ll_, li, r);
            return r;
         }
         return static_cast<const VeryLong&>(fvl) % li;
      }
      bool is_probable_prime(long int number_of_tests = 1) const
      {
         if (faster_)
         {
            VeryLong tmp(l_);
            return tmp.is_probable_prime(number_of_tests);
         }
         if (fast_)
         {
            VeryLong tmp(ll_);
            return tmp.is_probable_prime(number_of_tests);
         }
         return !(mpz_probab_prime_p(vl_, number_of_tests) == 0);
      }
      bool operator<= (long int li) const
      {
         if (faster_) return (l_ <= li);
         if (fast_) return (ll_ <= li);
         return (mpz_cmp_ui(vl_, li) <= 0);
      }
      bool operator== (long int li) const
      {
         if (faster_) return (li == l_);
         if (fast_) return (li == ll_);
         return (static_cast<const VeryLong&>(*this) == VeryLong(li));
      }
      bool operator!= (const VeryLong& vl) const
      {
         if (faster_) return (vl != VeryLong(l_));
         if (fast_) return (vl != VeryLong(ll_));
         return (vl != *this);
      }

      bool operator==(const FastVeryLong& fvl) const
      {
         if (faster_)
         {
            return fvl.faster_ && l_ == fvl.l_;
         }
         else if (fast_)
         {
            return fvl.fast_ && ll_ == fvl.ll_;
         }
         else
         {
            return !fvl.faster_ && !fvl.fast_ && vl_ == fvl.vl_;
         }
      }

      long int get_long() const
      {
         if (faster_) return l_;
         if (fast_) return static_cast<long int>(ll_);
         return VeryLong::get_long();
      }
      long long int get_long_long() const
      {
         if (faster_) return l_;
         if (fast_) return ll_;
         return VeryLong::get_long_long();
      }
      bool factorise_no_trial(std::vector<long int> * factors) const
      {
#if 0
         // 1185264170860630907
         const VeryLong x("1185264170860630907");
         if (!(x != *this))
         {
            int y = 0;
            ++y;
         }
         check();
         check1();
#endif
         if (faster_) return factorise_no_trial_l(l_, factors);
         else if (fast_) return factorise_no_trial_ll(ll_, factors);
         else
         {
            return false;
            std::vector<VeryLong> f;
            if (VeryLong::factorise_no_trial(&f))
            {
               for (std::vector<VeryLong>::const_iterator it = f.begin();
                     it != f.end();
                     ++it)
               {
                  if (*it > VeryLong(max_long_)) return false;
                  factors->push_back(it->get_long());
               }
               return true;
            }
         }
         return false;
      }
      bool factorise_no_trial(std::vector<VeryLong> * factors) const
      {
         if (faster_)
         {
            VeryLong tmp(l_);
            return tmp.factorise_no_trial(factors);
         }
         if (fast_)
         {
            VeryLong tmp(ll_);
            return tmp.factorise_no_trial(factors);
         }
         return VeryLong::factorise_no_trial(factors);
      }
      friend ostream& operator<< (ostream& os, const FastVeryLong& vl)
      {
         vl.display(os);
         return os;
      }

      static const VeryLong max_long_long_;
      static const long long int max_long_;
      bool is_very_long() const
      {
         return (!fast_ && !faster_);
         //return static_cast<const VeryLong&>(*this) >= max_long_long_;
      }
      bool is_fast() const
      {
         return fast_;
      }
      bool is_faster() const
      {
         return faster_;
      }

   private:
      mutable long long int ll_;
      mutable long int l_;
      mutable bool fast_;
      mutable bool faster_;
      void check1() const
      {
         if (fast_ && !faster_ && ll_ <= max_long_)
         {
            faster_ = true;
            l_ = static_cast<long int>(ll_);
         }
      }
      void check() const
      {
         if (!fast_ && static_cast<const VeryLong&>(*this) <= max_long_long_)
         {
            fast_ = true;
            mpz_t v;
            mpz_init(v);
            mpz_div_2exp(v, vl_, 32);
            ll_ = mpz_get_si(v);
            ll_ <<= 32;
            mpz_mod_2exp(v, vl_, 32);
            ll_ += mpz_get_ui(v);
            mpz_clear(v);
         }
      }
};
#else
typedef VeryLong FastVeryLong;
#endif
double ln(const VeryLong& vl);

#endif
