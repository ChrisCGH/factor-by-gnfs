#pragma GCC diagnostic ignored "-Wundefined-var-template"
#include "NumberField.h"
#include "AlgebraicNumber.h"
#include "AlgebraicNumber_in_O_pO.h"
#include "Ideal.h"
#include "pow.h"
#include "discriminant.h"
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include "MPFloat.h"
#include <fstream>
#include <float.h>
#include "timings.h"
#ifdef WIN32
#define isnan(x) _isnan(x)
#define isinf(x) !_finite(x)
#endif
#ifdef __SUNPRO_CC
template <>
void complex<MPFloat>::_div(const MPFloat& __z1_r, const MPFloat& __z1_i,
                        const MPFloat& __z2_r, const MPFloat& __z2_i,
                        MPFloat& __res_r, MPFloat& __res_i) {
  MPFloat __ar = __z2_r >= 0L ? __z2_r : -__z2_r;
  MPFloat __ai = __z2_i >= 0L ? __z2_i : -__z2_i;

  if (__ar <= __ai) {
    MPFloat __ratio = __z2_r / __z2_i;
    MPFloat __denom = __z2_i * (1.0 + __ratio * __ratio);
    __res_r = (__z1_r * __ratio + __z1_i) / __denom;
    __res_i = (__z1_i * __ratio - __z1_r) / __denom;
  }
  else {
    MPFloat __ratio = __z2_i / __z2_r;
    MPFloat __denom = __z2_r * (1.0 + __ratio * __ratio);
    __res_r = (__z1_r + __z1_i * __ratio) / __denom;
    __res_i = (__z1_i - __z1_r * __ratio) / __denom;
  }
}
#endif
#define DEBUG_OUTPUT 1
namespace
{
    bool verbose()
    {
        if (std::getenv("NUMBER_FIELD_VERBOSE"))
        {
            return true;
        }
        return false;
    }
int isSquareFree(const std::vector<std::pair<VeryLong, int > >& factors, VeryLong ignorePrimeSquare = 0L)
{
   for (size_t i = 0; i < factors.size(); i++)
   {
      if (factors[i].second > 1)
      {
         if (ignorePrimeSquare == 0L) return 0;
         if (factors[i].first == ignorePrimeSquare && factors[i].second > 3) return 0;
      }
   }
   return 1;
}

int isFundamentalDiscriminant(const std::vector<std::pair<VeryLong, int > >& factors, const VeryLong& product)
{
   const VeryLong four(4L);
   if (product % four == 1L && isSquareFree(factors)) return 1;
   if (product % four == 0L)
   {
      VeryLong x = (product / four) % four;
      if (x == 2L || x == 3L)
      {
         if (isSquareFree(factors, 2L)) return 1;
      }
   }

   return 0;
}

void factorDiscriminant(const VeryLong& D,
                        VeryLong& D0,
                        VeryLong& F,
                        std::vector<std::pair<VeryLong, int > >& FF_factors)
{
   std::vector<VeryLong> factors;
   std::vector<std::pair<VeryLong, int > > F_factors;
   VeryLong D1 = D;
   if (D1 < 0L) D1 = -D1;
   D1.factorise(&factors);
   VeryLong check(1L);
   for (size_t j = 0; j < factors.size(); j++)
   {
      check *= factors[j];
   }

   VeryLong q = D1 / check;
   std::sort(factors.begin(), factors.end());
   // Now decompose D as D0 * F^2 where D0 is 1 or a fundamental discriminant
   // To do this collect factors together and remove all even powers of primes > 2
   std::vector<std::pair<VeryLong, int > > collectedFactors;
   VeryLong prev(0L);
   D0 = 1L;
   F = 1L;
   int index = -1;
   for (auto& f: factors)
   {
      if (prev != f)
      {
         index++;
         collectedFactors.push_back(std::pair<VeryLong, int>(f, 1));
         F_factors.push_back(std::pair<VeryLong, int>(f, 0));
         prev = f;
      }
      else
      {
         collectedFactors[index].second = 1 - collectedFactors[index].second;
         if (collectedFactors[index].second == 0)
         {
            F *= f;
            F_factors[index].second += 1;
         }
      }
   }
   for (auto& cf: collectedFactors)
   {
      if (cf.second > 0)
      {
         D0 *= cf.first;
      }
   }
   if (!isFundamentalDiscriminant(collectedFactors, D0))
   {
      if (collectedFactors[0].first == 2L && F % 2L == 0L)
      {
         collectedFactors[0].second += 2;
         F /= 2L;
         F_factors[0].second -= 1;
         D0 *= 4L;
         if (!isFundamentalDiscriminant(collectedFactors, D0))
         {
            std::cout << "Problem, can't decompose D = " << D << std::endl;
         }
      }
   }
   for (auto& ff: F_factors)
   {
      if (ff.second > 0)
      {
         FF_factors.push_back(ff);
      }
   }
}

};

NumberField::NumberField()
      : structureMatrix_(1,1),
      integralBasisAlpha_(1,1),
      integralBasisAlphaInv_(1,1),
      integralBasisTheta_(1,1),
      integralBasisThetaInv_(1,1)
{}
NumberField::NumberField(const Polynomial<VeryLong>& poly, const char* fbFile)
      : structureMatrix_(2*poly.deg() - 1, poly.deg()),
      min_poly_(poly),
      integralBasisAlpha_(poly.deg(), poly.deg()),
      integralBasisAlphaInv_(poly.deg(), poly.deg()),
      integralBasisTheta_(poly.deg(), poly.deg()),
      integralBasisThetaInv_(poly.deg(), poly.deg())
{
   if (poly.deg() >= MAX_DEGREE)
   {
      std::ostringstream oss;
      oss << "degree of number field (" << poly.deg() << ") too big, must be < " << MAX_DEGREE;
      throw std::string(oss.str());
   }
   // Convert poly to Polynomial<complex<double > > to find (complex) roots
   MPFloat::set_precision(100);
   std::vector<complex<MPFloat > > co(poly.deg() + 1, complex<MPFloat>(0.0, 0.0));
   int i = 0;
   for (i = 0; i <= poly.deg(); i++)
   {
      co[i] = complex<MPFloat>((MPFloat)poly.coefficient(i).get_double(), MPFloat(0.0));
   }
   Polynomial<complex<MPFloat > > cpoly(co);
   std::vector<complex<MPFloat > > roots(0, complex<MPFloat>(0.0, 0.0));
   find_roots_over_C_q1<MPFloat>(cpoly, roots);
#ifdef DEBUG_OUTPUT
   if (verbose())
   {
       std::cout << "NumberField::NumberField() : min_poly_ = " << min_poly_ << std::endl;
       std::cout << "NumberField::NumberField() : roots are " << std::endl;
   }
#endif
   for (size_t j = 0; j < roots.size(); j++)
   {
#ifdef DEBUG_OUTPUT
      if (verbose())
      {
          std::cout << roots[j] << std::endl;
      }
#endif
      double a = roots[j].real();
      double b = roots[j].imag();
#ifdef DEBUG_OUTPUT
      if (verbose())
      {
          std::cout << a << " + " << b << "i" << std::endl;
      }
#endif
      roots_.push_back(complex<long double>((long double)a, (long double)b));
#ifdef DEBUG_OUTPUT
      if (verbose())
      {
          std::cout << std::setprecision(20) << roots_[j] << std::endl;
          complex<MPFloat> res = cpoly.evaluate(roots[j]);
          std::cout << "min_poly_(z) = " << res << std::endl;
          std::cout << "|min_poly_(z)|^2 = " << std::norm(res) << std::endl;
      }
#endif
   }

   // assume that 1, alpha, alpha^2, ..., alpha^(degree - 1)
   // is a basis for the NumberField.
   // Then we need to know how to express alpha^k in terms of
   // these, for k >= degree, i.e.
   // alpha^(k+degree) = sum (i = 0,degree-1) A(k,i) * alpha^i
   // for k = 0,...,degree-2
   // So set up a matrix of rational numbers structureMatrix_ s.t.
   // alpha^i = sum structureMatrix_(i,j) * alpha^j
   // i = 0,...,2*degree - 2
   // j = 0,...,degree - 1
   int d = degree();
   for (i = 0; i < d; i++)
   {
      for (int j = 0; j < d; j++)
      {
         structureMatrix_(i,j) = Quotient<VeryLong>(0L);
      }
      structureMatrix_(i,i) = Quotient<VeryLong>(1L);
   }
   for (int j = 0; j < d; j++)
   {
      structureMatrix_(d,j) = Quotient<VeryLong>(-1L * min_poly_.coefficient(j), c_d());
   }
   for (i = d+1; i < 2*d - 1; i++)
   {
      structureMatrix_(i,0) = structureMatrix_(i-1,d-1) * Quotient<VeryLong>(-1L * min_poly_.coefficient(0), c_d());
      for (int j = 1; j < d; j++)
      {
         structureMatrix_(i,j) = structureMatrix_(i-1,j-1) - structureMatrix_(i-1,d-1) * Quotient<VeryLong>(min_poly_.coefficient(j), c_d());
      }
   }

   discriminant_ = ::discriminant(poly);
#ifdef DEBUG_OUTPUT
   if (verbose())
   {
       std::cout << "Discriminant of " << poly << " is " << discriminant_ << std::endl;
   }
#endif

   // calculate factor base
   if (!fbFile)
   {
      factorBase_ = new FactorBase(min_poly_, 10000L, "fb.dat");
   }
   else
   {
      try
      {
         factorBase_ = new FactorBase(fbFile);
      }
      catch (...)
      {
         std::cerr << "NumberField::NumberField : caught exception" << std::endl;
         //delete factorBase_;
         factorBase_ = new FactorBase(min_poly_, 1600L, fbFile);
         //factorBase_->write(fbFile);
      }
   }

   AlgebraicNumber::setNumberField(*this);
   timing_file("round2.tim");
   timing_start("Round2");
   Round2();
   timing_stop();
#ifdef DEBUG_OUTPUT
   if (verbose())
   {
       std::cout << "integralBasisAlpha_ = " << std::endl << integralBasisAlpha_;
       std::cout << "integralBasisAlphaInv_ = " << std::endl << integralBasisAlphaInv_;
   }
#endif
}

NumberField::~NumberField()
{}

int NumberField::conjugates() const
{
   return roots_.size();
}

complex<long double > NumberField::conjugate(int r) const
{
   if (r < 0 || r >= static_cast<int>(roots_.size()))
   {
       throw std::string("NumberField::conjugate(): index out of range");
   }
   return complex<long double>((long double)(double)roots_[r].real(),
                               (long double)(double)roots_[r].imag());
}

long double NumberField::ln_sigma(int j, const VeryLong& a, const VeryLong& b) const
{
   // approximate a - alpha_j * b where alpha_j is jth conjugate
   const long double OVERFLOW_LIMIT = 1e100;
   long double aa = a.get_long_double();
   long double bb = b.get_long_double();
   //std::cerr << "NumberField::ln_sigma(), aa = " << aa << ", bb = " << bb << std::endl;
   long double delta = 1.0;
   if (fabs(aa) > OVERFLOW_LIMIT || fabs(bb) > OVERFLOW_LIMIT)
   {
      delta = OVERFLOW_LIMIT;
      aa /= OVERFLOW_LIMIT;
      bb /= OVERFLOW_LIMIT;
   }
   //std::cerr << "NumberField::ln_sigma(), j = " << j << ", aa = " << aa << ", bb = " << bb << ", delta = " << delta << std::endl;
   complex<long double> alpha_j = conjugate(j);
   //complex<long double> alpha_j_ = std::conj(alpha_j);
   complex<long double> x = complex<long double>(aa) - alpha_j * bb;
   complex<long double> x_ = std::conj(x);
   long double mod2 = (x * x_).real();
   //long double mod2 = aa * aa + bb * bb * std::norm(alpha_j) - 2.0 * aa * bb * alpha_j.real();
   //std::cerr << "NumberField::ln_sigma(), mod2 = " << mod2 << std::endl;
   long double answer = 0.0;
   if (mod2 > 0.0)
   {
      answer = (long double)log((double)mod2) / 2.0 + (long double)log(delta);
   }
   else
   {
	   complex<long double> r = complex<long double>(aa) - alpha_j * complex<long double>(bb);
	   std::cerr << "Warning! : a - b alpha_" << j << "  very close to zero :" << std::endl;
	   std::cerr << "(a,b) = (" << a << "," << b << ")" << std::endl;
	   std::cerr << "alpha_" << j << " = " << alpha_j << std::endl;
	   std::cerr << "a - b alpha_" << j << " = " << r << std::endl;
   }
   //std::cerr << "NumberField::ln_sigma(), answer = " << answer << std::endl;
   return answer;
}

// Algorithm 6.1.8 (Zassenhaus's Round 2)
void NumberField::Round2()
{
   //timing_file("round2.tim");
   // First find an algebraic integer and its minimal monic polynomial T
   // c_d * alpha is an algebraic integer, so use that:
#ifdef DEBUG_OUTPUT
   if (verbose())
   {
       std::cout << "Round 2 :" << std::endl;
   }
#endif
   int d = degree();
   VeryLong c_d = NumberField::c_d();
   AlgebraicNumber theta = c_d * AlgebraicNumber::alpha();
   const Quotient<VeryLong> qone = Quotient<VeryLong>(1L);
   Quotient<VeryLong> cf = qone;
   // minimal poly of theta is T(X) = F(X,c_d) / c_d
   std::vector<VeryLong> c;
   c.resize(d + 1);
   VeryLong factor = 1L;
   c[d] = factor;
   for (int i = d - 1; i >= 0; --i)
   {
      c[i] = min_poly_.coefficient(i) * factor;
      factor *= c_d;
   }
   Polynomial<VeryLong> T(c);
   monic_min_poly_ = T;
#ifdef DEBUG_OUTPUT
   if (verbose())
   {
       std::cout << "Monic minimal polynomial = " << monic_min_poly_ << std::endl;
   }
#endif

   // Step 1. [Factor discriminant of polynomial]
   VeryLong D = ::discriminant(T);
#ifdef DEBUG_OUTPUT
   if (verbose())
   {
       std::cout << "Discriminant of " << T << " is " << D << std::endl;
   }
#endif
   VeryLong D1 = D;
   if (D1 < 0L) D1 *= -1L;
   VeryLong F(1L);
   VeryLong D0(1L);
   std::vector<std::pair<VeryLong, int > > F_factors;
   factorDiscriminant(D, D0, F, F_factors);

   // Step 2. [Initialize]
   std::vector<AlgebraicNumber> omega;
   omega.resize(d);
   const AlgebraicNumber one(VeryLong(1L));
   AlgebraicNumber x;
   x = one;
   omega[0] = one;
   for (auto omegaIter = omega.begin() + 1;
         omegaIter != omega.end();
         ++omegaIter)
   {
      x *= theta;
      *omegaIter = x;
   }

   int done3 = 0;
   while (!done3)
   {
      // Step 3. [Loop on factors of F]
//      cout << "Step 3 : F = " << F << endl;
      if (F == 1L)
      {
         // The omega[i] form an integral basis
         // The way we have done the computations,
         // the coefficients of omega[i] are in
         // terms of alpha, but since omega[1] = theta = c_d * alpha
         // it's easy enough to convert.
         integralBasisAlpha_.set_size(d, d);
         integralBasisTheta_.set_size(d, d);
         Quotient<VeryLong> G = Quotient<VeryLong>(1L);

         cf = qone;
         for (int r = 0; r < d; r++)
         {
            for (int c = 0; c < d; c++)
            {
               integralBasisAlpha_(r,c) = omega[c].coefficients()[r];
               integralBasisTheta_(r,c) = integralBasisAlpha_(r,c) / cf;
               if (r == c)
               {
                  G *= integralBasisTheta_(r,c);
               }
            }
            cf *= c_d;
         }
         invert(integralBasisAlpha_, integralBasisAlphaInv_);
         invert(integralBasisTheta_, integralBasisThetaInv_);
         AlgebraicNumber_in_O_pO::set_basis(integralBasisAlpha_, omega);
         AlgebraicNumber_in_O_pO_1::set_basis(integralBasisAlpha_, omega);

#ifdef DEBUG_OUTPUT
         if (verbose())
         {
             for (int i = 0; i < d; i++)
             {
                std::cout << "omega[" << i << "] = " << omega[i] << std::endl;
                std::cout << "N(omega[" << i << "]) = " << omega[i].norm() << std::endl;
             }
         }
#endif

         Quotient<VeryLong> qd = D * G * G;
         fieldDiscriminant_ = qd.numerator();
#ifdef DEBUG_OUTPUT
         if (verbose())
         {
             std::cout << "G = " << G << ", qd = " << qd << ", d = " << fieldDiscriminant_ << std::endl;
         }
#endif
         index_ = G.denominator();

#ifdef DEBUG_OUTPUT
         if (verbose())
         {
             std::cout << "integralBasisTheta_ = " << std::endl << integralBasisTheta_;
             std::cout << "integralBasisThetaInv_ = " << std::endl << integralBasisThetaInv_;
         }
#endif

         return;
      }
//timing_start("Point 1");
      // get p, the smallest prime factor of F
      size_t ii = 0;
      while (ii < F_factors.size() && F_factors[ii].second == 0) ii++;
      if (ii >= F_factors.size())
      {
         std::cout << "Problem, F has no factors!" << std::endl;
         return;
      }
      VeryLong p = F_factors[ii].first;
      //std::cerr << "p = " << p << std::endl;
      VeryLong pd = 1L;
      for (int i = 0; i < d; i++) pd *= p;


      // Step 4. [Factor modulo p]
      // Factor T mod p into distinct irreducible polynomials
//timing_stop();
//timing_start("Point 2: Factor T mod p");
      VeryLongModular::set_default_modulus(p);
      if (verbose())
      {
          std::cout << "p = " << p << std::endl;
      }
      std::vector<VeryLongModular> cc;
      cc.resize(d + 1);
      for (int jj = 0; jj < d + 1; jj++)
      {
         VeryLong v = T.coefficient(jj) % p;
         cc[jj] = VeryLongModular(v);
      }
      Polynomial<VeryLongModular> T_(cc);

      std::vector<Polynomial<VeryLongModular > > ifactors;
      factor_over_F_p<VeryLong, VeryLong, VeryLongModular>(T, p, ifactors);
//timing_stop();
//timing_start("Point 3");
      std::vector<std::pair<Polynomial<VeryLongModular>, int> > dfactors;

      Polynomial<VeryLongModular> curr;
      Polynomial<VeryLongModular> prev;
      int multiplicity = 0;
      for (size_t i = 0; i < ifactors.size(); i++)
      {
         curr = ifactors[i];
         if (curr != prev && i != 0)
         {
            dfactors.push_back(std::pair<Polynomial<VeryLongModular>, int>(prev, multiplicity));
            multiplicity = 0;
         }
         prev = curr;
         multiplicity++;
      }
      dfactors.push_back(std::pair<Polynomial<VeryLongModular>, int>(prev, multiplicity));

      Polynomial<VeryLongModular> g_(1L);
      for (auto& df: dfactors)
      {
         g_ *= df.first;
      }

      Polynomial<VeryLongModular> h_ = T_ / g_;

      Polynomial<VeryLong> g = monic_lift<VeryLong, VeryLongModular>(g_);
      Polynomial<VeryLong> h = monic_lift<VeryLong, VeryLongModular>(h_);

      Polynomial<VeryLong> f = (g * h - T) / p;
      Polynomial<VeryLongModular> f_ = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(f, p);

      Polynomial<VeryLongModular> Z_ = gcd(f_, g_);
      if (verbose())
      {
         std::cout << "Z_ = " << Z_ << std::endl;
         std::cout << "h_ = " << h_ << std::endl;
      }
      Z_ = gcd(Z_, h_);
      if (verbose())
      {
         std::cout << "Z_ = " << Z_ << std::endl;
      }

      Polynomial<VeryLongModular> U_ = T_ / Z_;
      Polynomial<VeryLong> U = monic_lift<VeryLong, VeryLongModular>(U_);
      Polynomial<VeryLong> Z = monic_lift<VeryLong, VeryLongModular>(Z_);
      int m = Z.deg();
//timing_stop();
      // Step 5. [Apply Dedekind]
      if (m == 0)
      {
//timing_start("Point 5");
         while (F % p == 0L) F /= p;
         F_factors.erase(F_factors.begin());
//timing_stop();
      }
      else
      {
//timing_start("Point 6");
         Matrix<VeryLong> nu(d, d + m);
         Matrix<Quotient<VeryLong > > qnu(d, d + m);
         AlgebraicNumber U_at_theta(U, theta);

         VeryLong p_m_plus_1 = p;
         VeryLong lcm(1L);
         AlgebraicNumber tmp;
         Quotient<VeryLong> x;
         VeryLong y;
         for (int i = 0; i < m; i++)
         {
            tmp = U_at_theta;
            tmp *= omega[i];
            // Put coefficients of tmp (w.r.t. to theta) into column i of nu.
            // The way we have done the computations,
            // the coefficients of tmp are in
            // terms of alpha, but since omega[1] = theta = c_d * alpha
            // it's easy enough to convert.
            cf = qone;
            auto tmpIter = tmp.coefficients().begin();
            for (int r = 0; r < d; r++)
            {
               x = *tmpIter;
               x /= cf;
               qnu(r,i) = x;
               y = ::gcd(lcm, x.denominator());
               lcm *= x.denominator();
               lcm /= y;
               cf *= c_d;
               ++tmpIter;
            }
            p_m_plus_1 *= p;
         }
//timing_stop();
//timing_start("Point 7");
         for (int i = m; i < d + m; i++)
         {
            tmp = omega[i - m] * p;
            cf = qone;
            auto tmpIter = tmp.coefficients().begin();
            Quotient<VeryLong> x;
            for (int r = 0; r < d; r++)
            {
               x = *tmpIter;
               x /= cf;
               qnu(r,i) = x;
               y = ::gcd(lcm, x.denominator());
               lcm *= x.denominator();
               lcm /= y;
               cf *= c_d;
               ++tmpIter;
            }
         }
         for (int i = 0; i < d; i++)
         {
            for (size_t j = 0; j < qnu.columns(); j++)
            {
               x = qnu(i, j);
               nu(i, j) = x.numerator();
               nu(i, j) *= lcm;
               nu(i, j) /= x.denominator();
            }
         }
//timing_stop();
//timing_start("Point 8");

         // Apply Hermite reduction to nu
         //Matrix<VeryLong> H = HNF_mod_D(nu, D1*pd);
         Matrix<VeryLong> H = HNF(nu);

         // new omega[i]
         // make sure we keep omega in terms of alpha
         for (int i = 0; i < d; i++)
         {
            cf = qone;
            for (int j = 0; j < d; j++)
            {
               omega[i].set_coefficient(j, Quotient<VeryLong>(H(j,i), p * lcm) * cf);
               cf *= c_d;
            }
         }

//timing_stop();
         // Step 6. [Is new order p-maximal?]
         if (F % p_m_plus_1 != 0L)
         {
//timing_start("Point 8");
            while (F % p == 0L) F /= p;
            F_factors.erase(F_factors.begin());
//timing_stop();
         }
         else
         {
            int done7 = 0;
            while (!done7)
            {
//timing_start("Point 9");
               // Step 7. [Compute radical]
//               cout << "Step 7" << endl;
               VeryLong q = p;
               while (q < VeryLong((long)d)) q *= p;
//	       cout << "Step 7. q = " << q << endl;

               Matrix<Quotient<VeryLong> > W(d, d);
               Quotient<VeryLong> x;
               for (int c = 0; c < d; c++)
               {
                  int r = 0;
                  for (auto& r1: omega[c].coefficients())
                  {
                     x = r1;
                     W(r,c) = x;
                     r++;
                  }
               }

//timing_stop();

               AlgebraicNumber_in_O_pO::set_basis(W, omega, p);
               Matrix<VeryLongModular> Beta = AlgebraicNumber_in_O_pO::make_beta(q);

               int l = Beta.columns();
//timing_stop();
               // Step 8. [Compute new basis mod p]
               // Supplement basis of kernel of A

//timing_start("Point 12");
               Matrix<VeryLongModular> Beta_ = supplement(Beta);

//timing_stop();
//timing_start("Point 13");
               std::vector<AlgebraicNumber> beta;
               beta.resize(d);
               // Beta_ gives coefficients of beta[i] in terms of the omega[i]
               // Calculate beta[i] explicitly
               for (int i = 0; i < d; i++)
               {
                  beta[i] = AlgebraicNumber(VeryLong(0L));
                  for (int j = 0; j < d; j++)
                  {
                     beta[i] += Beta_(j,i).get_very_long() * omega[j];
                  }
               }
//timing_stop();

               // Step 9. [Compute big matrix]
               // Calculate alph[i] from lift of beta[i]
//timing_start("Point 14");
               std::vector<AlgebraicNumber> alph;
               alph.resize(d);
               for (int i = 0; i < d; i++)
               {
                  alph[i] = beta[i];
                  alph[i] += p * omega[i];
                  if (i + 1 > l)
                  {
                     alph[i] = p * alph[i];
                  }
               }
               Matrix<Quotient<VeryLong> > Alph(d, d);
               for (int i = 0; i < d; i++)
               {
                  for (int j = 0; j < d; j++)
                  {
                     Alph(i,j) = alph[j].coefficients()[i];
                  }
               }
               // Alph is coefficients of alph in terms of alpha[i]
//timing_stop();
//timing_start("Point 15");
               Matrix<Quotient<VeryLong > > Alph_inv(d,d);
               invert(Alph, Alph_inv);
//	       cout << "Alph = " << endl << Alph;
//	       cout << "Alph_inv = " << endl << Alph_inv;
               // Alph_inv is theta in terms of alph
//timing_stop();
//timing_start("Point 16.1");
               // Now calculate big matrix
               Matrix<VeryLongModular> C(d*d, d);
               Matrix<Quotient<VeryLong > > Cq(d*d, d);
//timing_stop();
               AlgebraicNumber tmp;
               for (int j = 0; j < d; j++)
               {
//timing_start("Point 16.2");
                  for (int k = 0; k < d; k++)
                  {
                     tmp = omega[k];
                     tmp *= alph[j];
                     // Get coefficients in terms of tmp in terms of
                     // alph[i]
                     Matrix<Quotient<VeryLong> > Cjk(d, 1);
                     Quotient<VeryLong> x;
                     int i = 0;
                     for (auto& i1: tmp.coefficients())
                     {
                        x = i1;
                        Cjk(i,0) = x;
                        i++;
                     }
                     Cjk = Alph_inv * Cjk;
                     // Cjk should now be coefficients of omega[j] * alpha[k] in
                     // terms of alph
                     for (int i = 0; i < d; i++)
                     {
                        Cq(i*d + j,k) = Cjk(i,0);
                        if (Cq(i*d + j,k).denominator() % p != 0L)
                        {
                           C(i*d + j,k) = VeryLongModular(Cq(i*d + j,k).numerator() % p) / VeryLongModular(Cq(i*d + j,k).denominator() % p);
                        }
                     }
                  }
//timing_stop();
               }
//timing_start("Point 17");

               // Step 10. [Compute new order]
               Matrix<VeryLongModular> Gamma = kernel(C);

//               cout << "Step 10. Kernel of C = " << endl << Gamma;

               int mm = Gamma.columns();
               // Gamma is in terms of omega[i]
               // so calculate gamma[i];
               std::vector<AlgebraicNumber> gamma;
               gamma.resize(mm);
               for (int i = 0; i < mm; i++)
               {
                  gamma[i] = AlgebraicNumber(VeryLong(0L));
                  for (int j = 0; j < d; j++)
                  {
                     gamma[i] += Gamma(j,i).get_very_long() * omega[j];
                  }
               }
//timing_stop();
//timing_start("Point 18");

               Matrix<VeryLong> nu1(d, d + mm);
               Matrix<Quotient<VeryLong > > qnu1(d, d + mm);
               VeryLong lcm(1L);
               VeryLong y;
               for (int i = 1; i <= mm; i++)
               {
                  cf = qone;
                  int r = 0;
                  for (auto& r1: gamma[i-1].coefficients())
                  {
                     x = r1;
                     x /= cf;
                     qnu1(r,i-1) = x;
                     y = ::gcd(lcm, x.denominator());
                     lcm *= x.denominator();
                     lcm /= y;
                     cf *= c_d;
                     r++;
                  }
               }
               for (int j = 1; j <= d; j++)
               {
                  cf = qone;
                  int r = 0;
                  for (auto& r1: omega[j-1].coefficients())
                  {
                     x = r1 * p;
                     x /= cf;
                     qnu1(r,mm+j-1) = x;
                     y = ::gcd(lcm, x.denominator());
                     lcm *= x.denominator();
                     lcm /= y;
                     cf *= c_d;
                     r++;
                  }
               }

               for (int i = 0; i < d; i++)
               {
                  for (int j = 0; j < d + mm; j++)
                  {
                     nu1(i,j) = qnu1(i,j).numerator() * (lcm / qnu1(i,j).denominator());
                  }
               }

               //Matrix<VeryLong> H = HNF_mod_D(nu1, D1*pd);
               Matrix<VeryLong> H = HNF(nu1);
//timing_stop();
//timing_start("Point 19");
               std::vector<AlgebraicNumber> omega1;
               omega1.resize(d);

               for (int i = 0; i < d; i++)
               {
                  omega1[i] = omega[i];
                  cf = qone;
                  for (int j = 0; j < d; j++)
                  {
                     omega1[i].set_coefficient(j, Quotient<VeryLong>(H(j,i), p * lcm) * cf);
                     cf *= c_d;
                  }
               }

               // Step 11. [Finished with p?]
               int ii = 0;
               while (ii < d && omega[ii] == omega1[ii]) ii++;
               if (ii < d)
               {
                  for (int i = 0; i < d; i++)
                  {
                     omega[i] = omega1[i];
                  }

               }
               else
               {
                  while (F % p == 0L) F /= p;
                  F_factors.erase(F_factors.begin());
                  done7 = 1;
               }
//timing_stop();
            }
         }
      }
   } // Step 3 loop
}

VeryLong NumberField::idealBound() const
{
   // this function calculates the constant C, which depends only
   // on the number field, and which satisfies
   // N(<delta>) <= C * N(I)
   // where delta is an algebraic number contained in the integral ideal
   // I, and <delta> is the principal ideal generated by delta.
   // For more details see:
   // Phong Nguyen, "A Montgomery-like Square Root for the Number Field Sieve",
   // Proceedings of ANTS-III

   VeryLong C;
   int d = degree();
   const std::vector<AlgebraicNumber>& omega = AlgebraicNumber::integralBasis();
   //long double best_sum = 0.0;
   double best_sum = 0.0;
   for (int j = 0; j < d; j++)
   {
      //long double sum = 0.0;
      double sum = 0.0;
      for (auto& om: omega)
      {
         sum += om.mod_sigma_2(j);
      }
      if (sum > best_sum) best_sum = sum;
   }

   //std::cout << "best_sum = " << best_sum << std::endl;
   //VeryLong C1((long double)sqrt(best_sum));
   VeryLong C1((double)sqrt(best_sum));
   //std::cout << "C1 = " << C1 << std::endl;
   VeryLong C3(1L);
   C3 += C1 * C1 * (long int)(d);
   C3 = C3.nth_root(2);

   VeryLong dd((long int)d);
   const VeryLong two(2L);
   VeryLong C2 = pow<VeryLong, int>(two, d + 1);
   //std::cout << "C2 = " << C2 << std::endl;
   C2 *= pow<VeryLong, int>(dd, d);
   //std::cout << "C2 = " << C2 << std::endl;
   VeryLong tmp = pow<VeryLong, int>(two, d * (d - 1) / 2);
   //std::cout << "tmp = " << tmp << std::endl;
   C2 *= tmp.nth_root(2);
   //std::cout << "C2 = " << C2 << std::endl;

   C = tmp * C2 * pow<VeryLong>(C3, dd);
   if (C1 > 1L)
   {
      C *= pow<VeryLong>(C1, dd);
   }
   //std::cout << "NumberField::idealBound() : C = " << C << std::endl;
   return C;
}

void NumberField::factorise_monic_min_poly_over_p(const VeryLong& p, std::vector<std::pair<Polynomial<VeryLongModular>, int> >& dfactors) const
{
   VeryLongModular::set_default_modulus(p);

   std::vector<Polynomial<VeryLongModular > > ifactors;
   factor_over_F_p<VeryLong, VeryLong, VeryLongModular>(monic_min_poly_, p, ifactors);
   Polynomial<VeryLongModular> curr;
   Polynomial<VeryLongModular> prev;
   int multiplicity = 0;
   for (size_t i = 0; i < ifactors.size(); i++)
   {
      curr = ifactors[i];
      if (curr != prev && i != 0)
      {
         dfactors.push_back(std::pair<Polynomial<VeryLongModular>, int>(prev, multiplicity));
         multiplicity = 0;
      }
      prev = curr;
      multiplicity++;
   }
   dfactors.push_back(std::pair<Polynomial<VeryLongModular>, int>(prev, multiplicity));
}
