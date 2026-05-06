#ifndef _NUMBERFIELD_H
#define _NUMBERFIELD_H
#include "VeryLong.h"
#include "VeryLongModular.h"
#include "gcd.h"
#include "Quotient.h"
#include "Polynomial.h"
#include "Matrix.h"
#include <complex>
#include <memory>
using std::complex;
#include "FactorBase.h"
#include <vector>

class NumberField
{
public:
    NumberField();
    NumberField(const Polynomial<VeryLong>& poly, const char* fbFile = 0);
    NumberField(const NumberField&) = delete;
    NumberField& operator=(const NumberField&) = delete;
    ~NumberField();

    int conjugates() const;
    complex<long double > conjugate(int r) const;
    long double ln_sigma(int j, const VeryLong& a, const VeryLong& b) const;
    void ln_sigma_all(const VeryLong& a, const VeryLong& b, std::vector<long double>& sigma_out) const;

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

    VeryLong index() const
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
    std::unique_ptr<FactorBase> factorBase_;
    Matrix<Quotient<VeryLong > > integralBasisAlpha_;
    Matrix<Quotient<VeryLong > > integralBasisAlphaInv_;
    Matrix<Quotient<VeryLong > > integralBasisTheta_;
    Matrix<Quotient<VeryLong > > integralBasisThetaInv_;
    void Round2();
};

#endif
