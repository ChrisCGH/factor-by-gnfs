#include "NumberField.h"
#include "AlgebraicNumber.h"
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
#else
#ifdef linux
#include <math.h>
#else
#include <ieeefp.h>
#ifdef sun
int isinf(double x) { return !finite(x) && x==x; }
#endif
#endif
#endif

const NumberField* AlgebraicNumber::numberField_ = 0;

AlgebraicNumber::AlgebraicNumber() : ibc_defined_(false), M_(AlgebraicNumber::degree(),AlgebraicNumber::degree()), matrix_defined_(false)
{
   c_.resize(AlgebraicNumber::degree());
}

AlgebraicNumber::AlgebraicNumber(std::vector<Quotient<VeryLong > >& coeff)
      : c_(coeff), ibc_defined_(false), M_(AlgebraicNumber::degree(), AlgebraicNumber::degree()), matrix_defined_(false)
{
   // what if coeff.size() > AlgebraicNumber::degree() ?
}

AlgebraicNumber::AlgebraicNumber(const VeryLong& v)
      : ibc_defined_(false), M_(AlgebraicNumber::degree(), AlgebraicNumber::degree()), matrix_defined_(false)
{
   c_.resize(AlgebraicNumber::degree());
   c_[0] = Quotient<VeryLong>(v);
   for (int i = 1; i < AlgebraicNumber::degree(); i++)
   {
      c_[i] = Quotient<VeryLong>(0L);
   }
}

AlgebraicNumber::AlgebraicNumber(const VeryLong& a, const VeryLong& b)
      : ibc_defined_(false), M_(AlgebraicNumber::degree(), AlgebraicNumber::degree()), matrix_defined_(false)
{
   c_.resize(AlgebraicNumber::degree());
   c_[0] = Quotient<VeryLong>(a);
   c_[1] = -Quotient<VeryLong>(b);
   for (int i = 2; i < AlgebraicNumber::degree(); i++)
   {
      c_[i] = Quotient<VeryLong>(0L);
   }
}

AlgebraicNumber::AlgebraicNumber(const Matrix<VeryLong>& M, const VeryLong& denominator, int column)
      : ibc_defined_(false), M_(AlgebraicNumber::degree(), AlgebraicNumber::degree()), matrix_defined_(false)
{
   c_.resize(M.rows());
   for (size_t i = 0; i < M.rows(); i++)
   {
      c_[i] = Quotient<VeryLong>(M(i, column), denominator);
   }
}

AlgebraicNumber::AlgebraicNumber(const Matrix<Quotient<VeryLong> >& M, int column)
      : ibc_defined_(false), M_(AlgebraicNumber::degree(), AlgebraicNumber::degree()), matrix_defined_(false)
{
   c_.resize(M.rows());
   for (size_t i = 0; i < M.rows(); i++)
   {
      c_[i] = M(i, column);
   }
}

AlgebraicNumber::AlgebraicNumber(const AlgebraicNumber& a)
      : c_(a.c_), ibc_(a.ibc_), ibc_defined_(a.ibc_defined_), M_(a.M_), matrix_defined_(a.matrix_defined_)
{}

AlgebraicNumber::AlgebraicNumber(const Polynomial<VeryLong>& f, const AlgebraicNumber& a)
      : ibc_defined_(false), M_(1,1), matrix_defined_(false)
{
   // Algebraic number got by evaluating f at a

   int d = f.deg();
   const AlgebraicNumber zero(VeryLong(0L));
   *this = zero;
   for (int i = d; i >= 0; --i)
   {
      *this *= a;
      *this += AlgebraicNumber(f.coefficient(i));
   }
}

AlgebraicNumber::AlgebraicNumber(const Polynomial<Quotient<VeryLong> >& f)
: ibc_defined_(false), M_(AlgebraicNumber::degree(),AlgebraicNumber::degree()), matrix_defined_(false)
{
    int d = f.deg(); 
    const AlgebraicNumber zero(VeryLong(0L));
    *this = zero;
    AlgebraicNumber a = alpha();
    for (int i = d; i >= 0; --i)
    {
       *this *= a;
       *this += AlgebraicNumber(f.coefficient(i));
    }
}

AlgebraicNumber::AlgebraicNumber(const Quotient<VeryLong>& v)
: ibc_defined_(false), M_(AlgebraicNumber::degree(),AlgebraicNumber::degree()), matrix_defined_(false)
{
   c_.resize(AlgebraicNumber::degree());
   c_[0] = v;
   for (int i = 1; i < AlgebraicNumber::degree(); i++)
   {
      c_[i] = Quotient<VeryLong>(0L);
   }
}

AlgebraicNumber::~AlgebraicNumber()
{}

AlgebraicNumber& AlgebraicNumber::operator=(const AlgebraicNumber& a)
{
   if (this != &a)
   {
      c_ = a.c_;
      ibc_defined_ = a.ibc_defined_;
      ibc_ = a.ibc_;
      matrix_defined_ = a.matrix_defined_;
      M_ = a.M_;
   }
   return *this;
}

AlgebraicNumber operator*(const AlgebraicNumber& a1,
                          const AlgebraicNumber& a2)
{
   // First compute c_k = Sum(i+j=k) a_i b_j
   std::vector<Quotient<VeryLong > > c;
   int d = AlgebraicNumber::degree();
   size_t csize = 2 * d - 1;
   if (c.size() != csize)
   {
      c.resize(csize);
   }
   const Quotient<VeryLong> zero(0L);
   Quotient<VeryLong> x;
   for (auto& i: c)
   {
      i = zero;
   }

   int i = 0;
   for (auto& i1: a1.c_)
   {
      int j = 0;
      for (auto& j1: a2.c_)
      {
         x = i1;
         x *= j1;
         c[i + j] += x;
         j++;
      }
      i++;
   }

   // now calculate coefficients for a1 * a2
   // as z_k = c_k + sum(i=0..n-2) r_k_i c_n+i

   AlgebraicNumber z;
   if ((int)z.c_.size() != d)
   {
      z.c_.resize(d);
      z.M_.set_size(d, d);
   }

   int k = 0;
   std::vector<Quotient<VeryLong> >::const_iterator c_iter = c.begin();
   for (auto& k1: z.c_)
   {
      k1 = *c_iter;
      std::vector<Quotient<VeryLong> >::const_iterator cc_iter = c.begin() + d;
      for (int j = 0; j < d - 1; j++)
      {
         x = (*cc_iter);
         x *= AlgebraicNumber::numberField_->structureMatrix()(j + d, k);
         k1 += x;
         ++cc_iter;
      }
      k++;
      ++c_iter;
   }

   return z;
}

AlgebraicNumber& AlgebraicNumber::operator*=(const VeryLong& v)
{
   for (auto& cc: c_)
   {
      cc *= v;
   }
   ibc_defined_ = false;
   return *this;
}

AlgebraicNumber& AlgebraicNumber::operator*=(long int p)
{
   for (auto& cc: c_)
   {
      cc *= p;
   }
   ibc_defined_ = false;
   return *this;
}

AlgebraicNumber& AlgebraicNumber::multiply_by_alpha_minus_r(long int r)
{
   // multiply by c_d (alpha - r)
   size_t d = AlgebraicNumber::degree();
   size_t csize = d + 1;
   std::vector<Quotient<VeryLong > > c(csize);

   int i = 0;
   for (auto& i1: c_)
   {
      Quotient<VeryLong> x(i1);
      x *= AlgebraicNumber::c_d();
      c[i+1] += x;
      x *= r;
      c[i] -= x;
      i++;
   }

   if (c_.size() != d)
   {
      c_.resize(d);
   }
   if (M_.rows() != d || M_.columns() != d)
   {
      M_.set_size(d, d);
   }
   for (size_t k = 0; k < c_.size(); ++k)
   {
      c_[k] = c[k] + c[d] * AlgebraicNumber::numberField_->structureMatrix()(d, k);
   }

   ibc_defined_ = false;
   return *this;
}

AlgebraicNumber& AlgebraicNumber::multiply(const VeryLong& a, const VeryLong& b)
{
   // multiply by a - b alpha
   size_t d = AlgebraicNumber::degree();
   size_t csize = d + 1;
   std::vector<Quotient<VeryLong > > c(csize);

   int i = 0;
   for (auto& i1: c_)
   {
      Quotient<VeryLong> x(i1);
      x *= a;
      c[i] += x;
      x = i1;
      x *= b;
      c[i+1] -= x;
      i++;
   }

   if (c_.size() != d)
   {
      c_.resize(d);
   }
   if (M_.rows() != d || M_.columns() != d)
   {
      M_.set_size(d, d);
   }
   for (size_t k = 0; k < c_.size(); ++k)
   {
      c_[k] = c[k] + c[d] * AlgebraicNumber::numberField_->structureMatrix()(d, k);
   }

   ibc_defined_ = false;
   return *this;
}

AlgebraicNumber& AlgebraicNumber::operator*=(const AlgebraicNumber& a)
{
   std::vector<Quotient<VeryLong > > c;
   size_t d = AlgebraicNumber::degree();
   size_t csize = 2 * d - 1;
   if (c.size() != csize)
   {
      c.resize(csize);
   }
   const Quotient<VeryLong> zero(0L);
   Quotient<VeryLong> x;
   for (auto& i1: c)
   {
      i1 = zero;
   }

   int i = 0;
   for (auto& i1: c_)
   {
      int j = 0;
      for (auto& j1: a.c_)
      {
         x = i1;
         x *= j1;
         c[i + j] += x;
         j++;
      }
      i++;
   }

   if (c_.size() != d)
   {
      c_.resize(d);
   }
   if (M_.rows() != d || M_.columns() != d)
   {
      M_.set_size(d, d);
   }

   int k = 0;
   std::vector<Quotient<VeryLong> >::const_iterator c_iter = c.begin();
   for (auto& k1: c_)
   {
      k1 = *c_iter;
      std::vector<Quotient<VeryLong> >::const_iterator cc_iter = c.begin() + d;
      for (size_t j = 0; j < d - 1; j++)
      {
         x = (*cc_iter);
         x *= AlgebraicNumber::numberField_->structureMatrix()(j + d, k);
         k1 += x;
         ++cc_iter;
      }
      k++;
      ++c_iter;
   }

   ibc_defined_ = false;
   return *this;
}

AlgebraicNumber operator*(const Quotient<VeryLong>& x,
                          const AlgebraicNumber& a)
{
   return a * x;
}

AlgebraicNumber operator*(const AlgebraicNumber& a,
                          const Quotient<VeryLong>& x)
{
   AlgebraicNumber c;
   int d = AlgebraicNumber::degree();
   if ((int)c.c_.size() != d) c.c_.resize(d);
   std::vector<Quotient<VeryLong> >::iterator iter = c.c_.begin();
   for (auto& i1: a.c_)
   {
      *iter = i1 * x;
      ++iter;
   }
   return c;
}

AlgebraicNumber operator+(const AlgebraicNumber& a1,
                          const AlgebraicNumber& a2)
{
   int n = AlgebraicNumber::degree();
   AlgebraicNumber c;
   if ((int)c.c_.size() != n) c.c_.resize(n);

   std::vector<Quotient<VeryLong> >::const_iterator a1_iter = a1.c_.begin();
   std::vector<Quotient<VeryLong> >::const_iterator a2_iter = a2.c_.begin();
   size_t i = 0;
   for (auto& i1: c.c_)
   {
      i1 = Quotient<VeryLong>(0L);
      if (i < a1.c_.size())
      {
         i1 += *a1_iter;
         ++a1_iter;
      }
      if (i < a2.c_.size())
      {
         i1 += *a2_iter;
         ++a2_iter;
      }
      i++;
   }

   return c;
}

AlgebraicNumber& AlgebraicNumber::operator+=(const AlgebraicNumber& a)
{
   std::vector<Quotient<VeryLong> >::const_iterator a_iter = a.c_.begin();
   int s = a.c_.size();
   int i = 0;
   for (auto& i1: c_)
   {
      if (i < s)
      {
         i1 += *a_iter;
         ++a_iter;
      }
      i++;
   }

   ibc_defined_ = false;
   return *this;
}

AlgebraicNumber operator/(const AlgebraicNumber& a1,
                          const AlgebraicNumber& a2)
{
   if (!a1.matrix_defined_) a1.defineMatrix();
   if (!a2.matrix_defined_) a2.defineMatrix();
   Matrix<Quotient<VeryLong > > m(a1.M_);
   Matrix<Quotient<VeryLong > > m2(a1.M_.rows(), a1.M_.columns());
   invert(a2.M_, m2);
   m = a1.M_ * m2;

   //cout << "a2.M_ = " << endl << a2.M_;
   //cout << "a2.M_^(-1) = " << endl << m2;
   //cout << endl << a2.M_ * m2;

   std::vector<Quotient<VeryLong > > c;
   c.resize(a1.c_.size());

   for (size_t i = 0; i < c.size(); i++)
   {
      c[i] = m(0, i);
   }

   return AlgebraicNumber(c);
}

AlgebraicNumber operator-(const AlgebraicNumber& a1,
                          const AlgebraicNumber& a2)
{
   return a1 + Quotient<VeryLong>(-1L) * a2;
}

Quotient<VeryLong > AlgebraicNumber::norm() const
{
   defineMatrix();
   return determinant(M_);
}

Quotient<VeryLong > AlgebraicNumber::trace() const
{
   defineMatrix();
   return M_.trace();
}

void AlgebraicNumber::defineMatrix() const
{
   if (c_.size() == 0) return;
   size_t d = AlgebraicNumber::degree();
   if (M_.rows() != d || M_.columns() != d) M_.set_size(d, d);
   AlgebraicNumber alp = alpha();
   AlgebraicNumber tmp = *this;
   for (size_t i = 0; i < d; i++)
   {
      for (size_t j = 0; j < d; j++)
      {
         M_(i, j) = tmp.c_[j];
      }

      if (i < d - 1) tmp = tmp * alp;
      //if (i < d - 1) tmp *= alp;
   }

//   cout << "M_ =" << endl << M_;
   matrix_defined_ = true;
}

AlgebraicNumber& AlgebraicNumber::alpha()
{
   static int first_time = 1;
   static AlgebraicNumber* alpha_ = 0;
   if (first_time)
   {
      first_time = 0;
      std::vector<Quotient<VeryLong > > c;
      c.resize(AlgebraicNumber::degree());
      c[0] = Quotient<VeryLong>(0L);
      c[1] = Quotient<VeryLong>(1L);
      for (int i = 2; i < AlgebraicNumber::degree(); i++)
      {
         c[i] = Quotient<VeryLong>(0L);
      }
      alpha_ = new AlgebraicNumber(c);
   }
   return *alpha_;
}

const std::vector<AlgebraicNumber>& AlgebraicNumber::integralBasis()
{
   static int first_time = 1;
   static std::vector<AlgebraicNumber> ib;
   if (first_time)
   {
      first_time = 0;
      for (int i = 0; i < AlgebraicNumber::degree(); i++)
      {
         ib.push_back(AlgebraicNumber(numberField_->w(), i));
      }
   }

   return ib;
}

int AlgebraicNumber::operator==(const AlgebraicNumber& a) const
{
   if (c_.size() != a.c_.size()) return 0;
   for (size_t i = 0; i < c_.size(); i++)
   {
      if (c_[i] != a.c_[i]) return 0;
   }
   return 1;
}

int AlgebraicNumber::operator!=(const AlgebraicNumber& a) const
{
   return !(*this == a);
}

void AlgebraicNumber::ln_sigma(int j, long double& ln_re, long int& re_sign,
                               long double& ln_im, long int& im_sign)
{
   if (!numberField_)
   {
      throw std::string("AlgebraicNumber::ln_sigma() : numberField_ not set");
   }
   const NumberField& nf = *AlgebraicNumber::numberField_;
   complex<long double> alpha_j = nf.conjugate(j);
   complex<long double> sigma = (long double)0.0;
   complex<long double> alpha_power = (long double)1.0;
   VeryLong coeff;
   VeryLong max_coeff(0L);
   try
   {
      for (auto& cc: c_)
      {
         VeryLong n = cc.numerator();
         VeryLong d = cc.denominator();
         long double sign = 1.0;
         if (n < 0L)
         {
            sign = -1.0;
            n = -n;
         }

         coeff = n / d;
         if (coeff > max_coeff) max_coeff = coeff;
         long double coeff_ld = sign * coeff.get_long_double();
         sigma += coeff_ld * alpha_power;
         if (isnan(sigma.real()) || isnan(sigma.imag())) throw std::overflow_error("NaN in complex multiply");
         if (isinf(sigma.real()) || isinf(sigma.imag())) throw std::overflow_error("Inf in complex multiply");
         alpha_power *= alpha_j;
      }
      re_sign = 1L;
      if (sigma.real() < 0) re_sign = -1L;
      ln_re = log(fabs(sigma.real()));
      im_sign = 1L;
      if (sigma.real() < 0) im_sign = -1L;
      ln_im = log(fabs(sigma.imag()));
      if (isnan(ln_im)) throw std::overflow_error("NaN in ln_im");

      return;
   }

   catch (std::overflow_error& oe)
   {
      std::cout << "Caught overflow exception: " << oe.what() << std::endl;
      double lg_coeff = log10(max_coeff);
      int power = (int)lg_coeff - 50;
      const VeryLong ten(10L);
      VeryLong scale = pow<VeryLong, int>(ten, power);
      long double log_scale_ld = ln(scale);
      complex<long double> sigma = (long double)0.0;
      complex<long double> alpha_power = (long double)1.0;
      for (auto& cc: c_)
      {
         VeryLong n = cc.numerator();
         VeryLong d = cc.denominator();
         long double sign = 1.0;
         if (n < 0L)
         {
            sign = -1;
            n = -n;
         }

         coeff = n / d;
         coeff /= scale;
         long double coeff_ld = sign * coeff.get_long_double();
         sigma += coeff_ld * alpha_power;
         alpha_power *= alpha_j;
      }
      re_sign = 1L;
      if (sigma.real() < 0) re_sign = -1L;
      ln_re = log(fabs(sigma.real())) + log_scale_ld;
      im_sign = 1L;
      if (sigma.real() < 0) im_sign = -1L;
      ln_im = log(fabs(sigma.imag())) + log_scale_ld;
      return;
   }
}

long double AlgebraicNumber::ln_sigma(int j)
{
   if (!numberField_)
   {
      throw std::string("AlgebraicNumber::ln_sigma() : numberField_ not set");
   }
   const NumberField& nf = *AlgebraicNumber::numberField_;
   complex<long double> alpha_j = nf.conjugate(j);
   complex<long double> sigma = (long double)0.0;
   complex<long double> alpha_power = (long double)1.0;
   VeryLong coeff;
   VeryLong max_coeff(0L);
   try
   {
      for (auto& cc: c_)
      {
         VeryLong n = cc.numerator();
         VeryLong d = cc.denominator();
         long double sign = 1.0;
         if (n < 0L)
         {
            sign = -1.0;
            n = -n;
         }

         coeff = n / d;
         if (coeff > max_coeff) max_coeff = coeff;
         long double coeff_ld = sign * coeff.get_long_double();
         sigma += coeff_ld * alpha_power;
         if (isnan(sigma.real()) || isnan(sigma.imag())) throw std::overflow_error("NaN in complex multiply");
         if (isinf(sigma.real()) || isinf(sigma.imag())) throw std::overflow_error("Inf in complex multiply");
         alpha_power *= alpha_j;
      }
      long double mod2 = std::norm(sigma);
      if (isinf(mod2)) throw std::overflow_error("Inf in modulus_squared");
      if (isnan(mod2)) throw std::overflow_error("NaN in modulus_squared");
      return (long double)log((double)mod2) / 2.0;
   }

   catch (std::overflow_error& oe)
   {
      std::cout << "Caught overflow exception: " << oe.what() << std::endl;
      double lg_coeff = log10(max_coeff);
      int power = (int)lg_coeff - 50;
      const VeryLong ten(10L);
      VeryLong scale = pow<VeryLong, int>(ten, power);
      long double log_scale_ld = ln(scale);
      complex<long double> sigma = (long double)0.0;
      complex<long double> alpha_power = (long double)1.0;
      for (auto& cc: c_)
      {
         VeryLong n = cc.numerator();
         VeryLong d = cc.denominator();
         long double sign = 1.0;
         if (n < 0L)
         {
            sign = -1;
            n = -n;
         }

         coeff = n / d;
         coeff /= scale;
         long double coeff_ld = sign * coeff.get_long_double();
         sigma += coeff_ld * alpha_power;
         alpha_power *= alpha_j;
      }
      long double mod2 = std::norm(sigma);
      return (long double)log((double)mod2) / 2.0 + log_scale_ld;
   }
}

long double AlgebraicNumber::mod_sigma_2(int j) const
{
   if (!numberField_)
   {
      throw std::string("AlgebraicNumber::mod_sigma_2() : numberField_ not set");
   }
   const NumberField& nf = *AlgebraicNumber::numberField_;
   complex<long double> alpha_j = nf.conjugate(j);
   complex<long double> sigma = (long double)0.0;
   complex<long double> alpha_power = (long double)1.0;
   for (auto& cc: c_)
   {
      long double coeff = cc.numerator().get_long_double() / cc.denominator().get_long_double();
      sigma += coeff * alpha_power;
      alpha_power *= alpha_j;
   }
   long double mod2 = std::norm(sigma);
   return mod2;
}

Polynomial<VeryLong> AlgebraicNumber::minimalPolynomial() const
{
   //cout << "Calculating minimal polynomial of " << *this << endl;
   int d = degree();
   // Now compute successive powers to find minimal polynomial
   int power = 2;
   Matrix<Quotient<VeryLong> > AA(d, power);
   const VeryLong one(1L);
   const VeryLong zero(0L);
   AA(0,0) = one;
   for (int row = 1; row < d; row++) AA(row,0) = zero;
   for (int row = 0; row < d; row++)
   {
      AA(row,1) = c_[row];
   }
   AlgebraicNumber x_pow = *this;
   Matrix<Quotient<VeryLong> > kerAA = kernel(AA);
   while (kerAA.columns() == 0)
   {
      x_pow = x_pow * (*this);
      power++;
      AA.add_column();
      for (int row = 0; row < d; row++)
      {
         AA(row,power-1) = x_pow.c_[row];
      }
      kerAA = kernel(AA);
   }

   //cout << "AA = " << endl << AA;
   //cout << "kerAA = " << endl << kerAA;
   // columns of kerAA now give coefficients of relations between powers of alph
   // find the minimal polynomial by looking for the column with least highest
   // non-zero row
   size_t least_highest_row = kerAA.rows();
   size_t min_col = 0;
   for (size_t col = 0; col < kerAA.columns(); col++)
   {
      size_t row = kerAA.rows() - 1;
      while (kerAA(row, col).is_zero()) --row;
      if (row < least_highest_row)
      {
         least_highest_row = row;
         min_col = col;
      }
   }

   // Construct minimal polynomial
   std::vector<VeryLong> coeff;
   coeff.resize(least_highest_row+1);
   VeryLong max_denominator(0L);
   for (size_t row = 0; row < least_highest_row + 1; row++)
   {
      if (kerAA(row, min_col).denominator() > max_denominator)
         max_denominator = kerAA(row, min_col).denominator();
   }
   for (size_t row = 0; row < least_highest_row + 1; row++)
   {
      Quotient<VeryLong> tmp = kerAA(row, min_col) * max_denominator;
      coeff[row] = tmp.numerator();
   }
   Polynomial<VeryLong> m(coeff);
   m.make_primitive();
   //cout << "Minimal polynomial is " << m << endl;
   AlgebraicNumber check(m, *this);
   //cout << "check = " << check << endl;
   return m;
}

AlgebraicNumber AlgebraicNumber::sqrt() const
{
   Polynomial<VeryLong> A = minimalPolynomial();
   std::vector<VeryLong> c;
   c.resize(3);
   c[0] = 0L;
   c[1] = 0L;
   c[2] = 1L;
   Polynomial<VeryLong> X2(c);
   Polynomial<VeryLong> AX2 = A.evaluate(X2);

   //cout << "A(X^2) = " << AX2 << endl;

   std::vector<Polynomial<VeryLong> > factors;
   VeryLong cont;
   Polynomial<VeryLong>::factor(AX2, factors, cont);

   //cout << "Factors of A(X^2) are : " << endl;
   //for (size_t i = 0; i < factors.size(); i++)
   //{
   //cout << factors[i] << endl;
   //}

   if (factors.size() == 2)
   {
      std::vector<AlgebraicNumber> s;
      s.resize(factors[0].deg() + 1);
      for (int i = 0; i <= factors[0].deg(); i++)
      {
         s[i] = AlgebraicNumber(factors[0].coefficient(i));
      }
      Polynomial<AlgebraicNumber> S(s);
      //cout << "S = " << S << endl;
      s.resize(3);
      s[0] = - *this;
      s[1] = AlgebraicNumber(0L);
      s[2] = AlgebraicNumber(1L);
      Polynomial<AlgebraicNumber> X2_minus_x(s);
      //cout << "X^2 - x = " << X2_minus_x << endl;
      Polynomial<AlgebraicNumber> Q;
      Polynomial<AlgebraicNumber> R;
      euclidean_division(S, X2_minus_x, Q, R);
      //cout << "result of euclidean division :" << endl;
      //cout << "Q = " << Q << endl;
      //cout << "R = " << R << endl;
      //cout << "Q * (X^2 - x) + R = " << Q * X2_minus_x + R << endl;
      if (R.deg() != 1)
      {
         std::cout << "Problem : R is not degree 1 : " << R << std::endl;
      }

      AlgebraicNumber a = R.coefficient(1);
      AlgebraicNumber b = R.coefficient(0);
      return -b / a;
   }

   return AlgebraicNumber(VeryLong(0L));
}

//
// if beta = Sum c_i alpha ^ (i - 1) with the c_i rational
// then Phi is the homomorphism to Z / NZ which maps alpha to m mod N
//
VeryLongModular AlgebraicNumber::Phi(const VeryLong& N, const VeryLong& m) const
{
   VeryLongModular::set_default_modulus(N);
   int d = degree();
   VeryLongModular result(0L);
   VeryLongModular zero(0L);
   VeryLongModular m_(m);
   for (int i = d - 1; i >= 0; --i)
   {
      Quotient<VeryLong> q = c_[i];
      VeryLongModular numer(q.numerator());
      VeryLongModular denom(q.denominator());
      //VeryLongModular vlm = VeryLongModular(q.numerator()) / VeryLongModular(q.denominator());
      if (denom == zero)
      {
         std::cout << "Problem in Phi(" << *this << ")" << std::endl;
         std::cout << "c_[" << i << "] = " << c_[i] << std::endl;
      }
      VeryLongModular vlm = numer / denom;
      result *= m_;
      result += vlm;
   }
   return result;
}

void AlgebraicNumber::createSpecialBasis(const VeryLong& p, std::vector<AlgebraicNumber>& j)
{
   if (!numberField_)
   {
      throw std::string("AlgebraicNumber::createSpecialBasis() : numberField_ not set");
   }
   j.push_back(AlgebraicNumber(p));
   AlgebraicNumber tmp(AlgebraicNumber::c_d());
   for (int i = AlgebraicNumber::degree() - 1; i >= 1; --i)
   {
      tmp = tmp * AlgebraicNumber::alpha() + AlgebraicNumber(AlgebraicNumber::numberField_->min_poly().coefficient(i));
      j.push_back(tmp);
   }
}
