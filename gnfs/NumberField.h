#ifndef _NUMBERFIELD_H
#define _NUMBERFIELD_H
#include "VeryLong.h"
#include "VeryLongModular.h"
#include "gcd.h"
#include "Quotient.h"
#include "Polynomial.h"
#include "Matrix.h"
#include <complex>
using std::complex;
#include "FactorBase.h"
#include <vector>

class NumberField
{
   public:
      NumberField();
      NumberField(const Polynomial<VeryLong>& poly, const char* fbFile = 0);
#if 0
private:
      NumberField(const NumberField& nf);
      NumberField& operator=(const NumberField& nf);
public:
#endif
      ~NumberField();

      int conjugates() const;
      complex<long double > conjugate(int r) const;
      long double ln_sigma(int j, const VeryLong& a, const VeryLong& b) const;

      int degree() const
      {
         return min_poly_.deg();
      }
      VeryLong c_d() const
      {
         return min_poly_.coefficient(min_poly_.deg());
      }

      const Matrix<Quotient<VeryLong > >& structureMatrix() const
      {
         return structureMatrix_;
      }

      FactorBase& factorBase() const
      {
         return *factorBase_;
      }

      const VeryLong index() const
      {
         return index_;
      }

      const Polynomial<VeryLong>& min_poly() const
      {
         return min_poly_;
      }

      const Polynomial<VeryLong>& monic_min_poly() const
      {
         return monic_min_poly_;
      }

      VeryLong discriminant() const
      {
         return discriminant_;
      }
      VeryLong fieldDiscriminant() const
      {
         return fieldDiscriminant_;
      }

      const Matrix<Quotient<VeryLong> >& w() const
      {
         return integralBasisAlpha_;
      }
      const Matrix<Quotient<VeryLong> >& winv() const
      {
         return integralBasisAlphaInv_;
      }
      const Matrix<Quotient<VeryLong> >& W() const
      {
         return integralBasisTheta_;
      }
      const Matrix<Quotient<VeryLong> >& Winv() const
      {
         return integralBasisThetaInv_;
      }

      VeryLong idealBound() const;

      void factorise_monic_min_poly_over_p(const VeryLong& p, std::vector<std::pair<Polynomial<VeryLongModular>, int> >& dfactors) const;
      enum { MAX_DEGREE = 10 };
   private:
      Matrix<Quotient<VeryLong > > structureMatrix_;
      Polynomial<VeryLong> min_poly_;
      Polynomial<VeryLong> monic_min_poly_;
      std::vector<complex<long double > > roots_;
      VeryLong discriminant_;
      VeryLong fieldDiscriminant_;
      VeryLong index_;
      FactorBase* factorBase_;
      Matrix<Quotient<VeryLong > > integralBasisAlpha_;
      Matrix<Quotient<VeryLong > > integralBasisAlphaInv_;
      Matrix<Quotient<VeryLong > > integralBasisTheta_;
      Matrix<Quotient<VeryLong > > integralBasisThetaInv_;
      void Round2();
};

#endif
