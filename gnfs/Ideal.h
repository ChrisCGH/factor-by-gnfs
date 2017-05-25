#ifndef __IDEAL_H
#define __IDEAL_H

#include "NumberField.h"
#include "AlgebraicNumber.h"
#include "AlgebraicNumber_in_O_pO.h"
#include "FactorBase.h"
#include "timings.h"

class Ideal
{
public:
    Ideal();
    Ideal(const VeryLong& p, const VeryLong& q);
    Ideal(const VeryLong& p, const AlgebraicNumber& beta);
    Ideal(const AlgebraicNumber& a);
    Ideal(const Matrix<VeryLong>& basis);
    Ideal(const Ideal& id);
    Ideal(const std::vector<AlgebraicNumber>& generator);
    Ideal(const Matrix<VeryLong>& hnf, const VeryLong& denom);
    ~Ideal();
    Ideal& operator=(const Ideal& id);
    int operator==(const Ideal& id) const;
    int operator!=(const Ideal& id) const
    {
        return !(*this == id);
    }

    friend Ideal operator*(const Ideal& I, const Ideal& J);
    friend Ideal operator/(const Ideal& I, const Ideal& J);
    Ideal invert() const;
    friend Ideal intersection(const Ideal& I, const Ideal& J);
    friend ostream& operator<<(ostream& os, const Ideal& I);

    Ideal& operator*=(const VeryLong& v);
    Ideal& operator/=(const VeryLong& v);

    Matrix<VeryLong> reducedBasisOmega() const;
    Quotient<VeryLong> norm() const;
    bool isPrincipal() const
    {
        return isPrincipal_;
    }

    const VeryLong& denominator() const
    {
        return denominator_;
    }
    const Matrix<VeryLong>& hnf_basis() const
    {
        return hnf_basis_;
    }

    static void integralPart(const Ideal& I, VeryLong& lcm, Matrix<VeryLong>& AA);
    bool isIntegral() const;

protected:
    Matrix<VeryLong> hnf_basis_;
    VeryLong denominator_;
    bool isPrincipal_;
#ifdef XYZ
    // a_ and b_ are set for principal ideals
    VeryLong a_;
    VeryLong b_;
#endif

private:
    void create(const std::vector<AlgebraicNumber>& generator);
    static Timing* timing_;
};

// class to represent ideal of the form I/pO where
// p is a prime and O is the ring of integers (given by
// an integral basis omega which has already been
// calculated by the Round2 algorithm.
class Ideal_mod_pO
{
public:
    Ideal_mod_pO(const Matrix<VeryLongModular>& F_p_basis);
    Ideal_mod_pO(const Ideal_mod_pO& IpO);
    ~Ideal_mod_pO();

    Ideal_mod_pO& operator=(const Ideal_mod_pO& IpO);

    bool operator==(const Ideal_mod_pO& rhs);

    friend Ideal_mod_pO operator*(const Ideal_mod_pO& IpO1,
                                  const Ideal_mod_pO& IpO2);
    friend Ideal_mod_pO operator/(const Ideal_mod_pO& IpO1,
                                  const Ideal_mod_pO& IpO2);
    Ideal makeIdeal();
    static void set_basis(const VeryLong& p);
    Matrix<VeryLongModular> basis() const
    {
        return F_p_basis_;
    }
    int rank()
    {
        return F_p_basis_.columns();
    }

private:
    static VeryLong p_;
    Matrix<VeryLongModular> F_p_basis_;

};

class SeparableAlgebraElement
{
public:
    SeparableAlgebraElement(const std::vector<VeryLongModular>& F_p_basis) : F_p_basis_(F_p_basis)
    {}
    SeparableAlgebraElement(const SeparableAlgebraElement& sa) : F_p_basis_(sa.F_p_basis_)
    {}
    ~SeparableAlgebraElement()
    {}

    SeparableAlgebraElement& operator=(const SeparableAlgebraElement& sa);

    friend SeparableAlgebraElement operator*(const SeparableAlgebraElement& sa1,
            const SeparableAlgebraElement& sa2);
    SeparableAlgebraElement& operator*=(const SeparableAlgebraElement& sa);
    friend SeparableAlgebraElement operator-(const SeparableAlgebraElement& sa1,
            const SeparableAlgebraElement& sa2);
    friend SeparableAlgebraElement operator+(const SeparableAlgebraElement& sa1,
            const SeparableAlgebraElement& sa2);
    friend SeparableAlgebraElement operator+(const SeparableAlgebraElement& sa1,
            const VeryLongModular& vlm);

    static void set_multiplication_table(const Matrix<VeryLongModular>& A)
    {
        A_ = A;
    }

    friend std::ostream& operator<< (std::ostream& os, const SeparableAlgebraElement& sae)
    {
        for (size_t row = 0; row < sae.F_p_basis_.size(); ++row)
        {
            os << sae.F_p_basis_[row].get_very_long() << " ";
        }
        return os;
    }
    const std::vector<VeryLongModular>& basis() const
    {
        return F_p_basis_;
    }

private:
    static Matrix<VeryLongModular> A_;
    std::vector<VeryLongModular> F_p_basis_;

};

class PrimeIdealRep;

template <class Z>
class PrimeIdealT : public Ideal
{
    // A prime ideal is an Ideal, but can also be represented as
    // a quintuplet (p_, alpha_, e_, f_, beta_)
    // where the PrimeIdeal is above the prime p_,
    // (p_, alpha_) is the 2-element representation,
    // e_ is the ramification index and f_ the residual index
public:
    PrimeIdealT() : beta_computed_(0) {}
    PrimeIdealT(const Z& p, const AlgebraicNumber& a)
        : Ideal(p, a), e_(0), beta_computed_(0), p_(p)
    {
        alpha_ = a;
        Z n = norm().numerator();
        f_ = 0;
        while (n % p == 0L)
        {
            n /= p;
            f_++;
        }
    }

    PrimeIdealT(const Z& p, const Z& q)
        : Ideal(p, q), e_(0), beta_computed_(0), p_(p)
    {
        VeryLong c_d = AlgebraicNumber::c_d();
        alpha_ = (AlgebraicNumber::alpha() - AlgebraicNumber(q)) * c_d;
        Z n = norm().numerator();
        f_ = 0;
        while (n % p == 0L)
        {
            n /= p;
            f_++;
        }
    }

    PrimeIdealT(const Z& p, const std::vector<AlgebraicNumber>& gamma)
        : Ideal(gamma), beta_computed_(0), p_(p)
    {
        std::vector<AlgebraicNumber> gamma1;
        gamma1.resize(hnf_basis_.columns());
        for (size_t i = 0; i < hnf_basis_.columns(); i++)
        {
            gamma1[i] = AlgebraicNumber(hnf_basis_, denominator_, i);
        }
        computeTwoElementRep(gamma1);
    }

    PrimeIdealT(const PrimeIdealT& pi)
        : Ideal(pi)
    {
        p_ = pi.p_;
        alpha_ = pi.alpha_;
        e_ = pi.e_;
        f_ = pi.f_;
        beta_ = pi.beta_;
        beta_omega_ = pi.beta_omega_;
        beta_computed_ = pi.beta_computed_;
    }

    ~PrimeIdealT() {}

    // Based on Theorem 4.8.13, decompose pZ into product of prime ideals above p
    static void primeDecompositionEasy(const Z& p,
                                       std::vector<std::pair<PrimeIdealT, int> >& primeIdeal)
    {
        static VeryLongModular zero(0L);
        static VeryLongModular one(1L);
        // For the moment, only the easy case
        if (AlgebraicNumber::index() % p == 0L)
        {
            std::cout << "p = " << p << " divides the index = " << AlgebraicNumber::index() << std::endl;
            return;
        }

        AlgebraicNumber theta = AlgebraicNumber::c_d() * AlgebraicNumber::alpha();

        // Factor T mod p into distinct irreducible polynomials
        std::vector<std::pair<Polynomial<VeryLongModular>, int> > dfactors;
        const NumberField& nf = AlgebraicNumber::nf();
        nf.factorise_monic_min_poly_over_p(p, dfactors);

        for (size_t i = 0; i < dfactors.size(); i++)
        {
            //cout << dfactors[i].first << " : " << dfactors[i].second << endl;
            Polynomial<VeryLong> Ti = monic_lift<VeryLong, VeryLongModular>(dfactors[i].first);
            //cout << "Ti = " << Ti << endl;
            AlgebraicNumber Ti_theta(Ti, theta);
            //cout << "Ti(theta) = " << endl << Ti_theta << endl;
            Ideal I(p, Ti_theta);
            //cout << "I(" << p << ", Ti_theta) = " << endl << I << endl;
            //cout << "N(I) = " << I.norm() << endl;

            PrimeIdealT PI(p, Ti_theta);
            PI.e_ = dfactors[i].second;
            primeIdeal.push_back(std::pair<PrimeIdealT, int>(PI, dfactors[i].second));
        }
        return;
    }

    // Algorithm 6.2.9 (Prime Decomposition)
    static void primeDecomposition(const Z& p,
                                   std::vector<std::pair<PrimeIdealT, int> >& primeIdeal)
    {
        static VeryLongModular zero(0L);
        static VeryLongModular one(1L);
        // Step 1. [Check if easy]
        if (AlgebraicNumber::index() % p != 0L)
        {
            PrimeIdealT::primeDecompositionEasy(p, primeIdeal);
            return;
        }

        bool debug = false;
        VeryLongModular::set_default_modulus(p);
        //std::cerr << "primeDecomposition(" << p << "), p | index = " << AlgebraicNumber::index() << std::endl;
        int degree = AlgebraicNumber::degree();
        const std::vector<AlgebraicNumber>& omega = AlgebraicNumber::integralBasis();

        // Step 2. [Compute radical]
        // note this code is the same as step 7 of NumberField::Round2() -
        // maybe we can factor it out?
        Z q = p;
        while (q < Z((long)degree)) q *= p;

        AlgebraicNumber_in_O_pO::set_basis(p);
        Matrix<VeryLongModular> Beta = AlgebraicNumber_in_O_pO::make_beta(q);
        //std::cerr << "Beta = " << Beta << std::endl;

        // Step 3. [Compute K_i]
        if (debug) std::cout << "Step 3 - Compute K_i" << std::endl;
        Ideal_mod_pO::set_basis(p);
        std::vector<Ideal_mod_pO> K;
        Ideal_mod_pO K_i(Beta);
        K.push_back(K_i);
        int i = 1;
        while (K_i.rank() != 0)
        {
            i++;
            K_i = K[0] * K[i-2];
            //std::cerr << "K_i = " << K_i.basis() << std::endl;
            K.push_back(K_i);
        }

        if (debug) std::cout << "K has " << K.size() << " entries" << std::endl;
        if (debug)
        {
            for (size_t ii = 0; ii < K.size(); ii++)
            {
                std::cout << "K[" << ii << "] = " << std::endl << K[ii].basis();
                K[ii].makeIdeal();
            }
        }

        // Step 4. [Compute J_j]
        std::vector<Ideal_mod_pO> J;
        J.push_back(K[0]);
        for (int j = 2; j <= i; j++)
        {
            Ideal_mod_pO J_j = K[j-1] / K[j-2];
            J.push_back(J_j);
        }

        if (debug) std::cout << "J has " << J.size() << " entries" << std::endl;
        if (debug)
        {
            for (size_t ii = 0; ii < J.size(); ii++)
            {
                std::cout << "J[" << ii << "] = " << std::endl << J[ii].basis();
                J[ii].makeIdeal();
            }
        }

        // Step 5. [Compute H_j]
        std::vector<Ideal_mod_pO> H;
        for (int j = 1; j <= i - 1; j++)
        {
            Ideal_mod_pO H_j = J[j-1] / J[j];
            H.push_back(H_j);
        }
        H.push_back(J[i-1]);

        if (debug) std::cout << "H has " << H.size() << " entries" << std::endl;
        if (debug)
        {
            for (size_t ii = 0; ii < H.size(); ii++)
            {
                std::cout << "H[" << ii << "] = " << std::endl << H[ii].basis();
                H[ii].makeIdeal();
            }
        }

        // Step 6. [Initialize loop]
        int j = 0;
        int c = 0;
        std::vector<Ideal_mod_pO> L;

        while (1)
        {
            // Step 7. [Finished?]
            //cout << "Step 7" << endl;
            int done7 = 0;
            while (!done7)
            {
                if (debug) std::cout << "c = " << c << ", i = " << i << ", j = " << j << std::endl;
                if (c == 0)
                {
                    if (j == i) return;
                    j++;
                    if (H[j-1].rank() < degree)
                    {
                        L.clear();
                        L.push_back(H[j-1]);
                        c = 1;
                        done7 = 1;
                    }
                }
                else done7 = 1;
            }

            int done8 = 0;
            while (!done8)
            {
                // Step 8. [Compute separable algebra A]
                if (debug) std::cout << "Step 8 : L has " << L.size() << " elements" << std::endl;
                Ideal_mod_pO H_ = L[0];
                Matrix<VeryLongModular> H_basis = H_.basis();
                if (debug) std::cout << "H_basis = " << std::endl << H_basis;
                int r = H_basis.columns();
                Matrix<VeryLongModular> B(degree, r + 1);
                for (int row = 0; row < degree; row++)
                {
                    for (int col = 0; col < r; col++)
                    {
                        B(row,col) = H_basis(row,col);
                    }
                    B(row, r) = zero;
                }
                B(0,r) = one;

                Matrix<VeryLongModular> B1 = supplement(B);
                int f = degree - r;
                if (debug) std::cout << "B1 = " << std::endl << B1;
                if (debug) std::cout << "f = " << f << ", r = " << r << std::endl;

                Matrix<VeryLongModular> Gamma(degree, f);
                for (int row = 0; row < degree; row++)
                {
                    for (int col = 0; col < f; col++)
                    {
                        Gamma(row,col) = B1(row, r + col);
                    }
                }
                if (debug) std::cout << "Gamma = " << std::endl << Gamma;
                // Gamma is matrix whose columns define F_p basis for
                // separable algebra A = O/H_ = (O/pO)/(H_/pO)

                // Step 9. [Compute multiplication table]
                // Construct a degree x f^2 matrix V which gives the
                // multiplication table of the columns of Gamma (basis for the
                // separable algebra A) in terms of the integral basis omega.
                // We use the already computed multiplication table of the omega,
                // W_mult_, to do this.
                Matrix<VeryLongModular> V(degree, f*f);
                for (int m = 0; m < degree; m++)
                {
                    int ik = 0;
                    for (int ii = 0; ii < f; ii++)
                    {
                        for (int k = 0; k < f; k++)
                        {
                            Quotient<VeryLong> V_mik(0L);
                            int jl = 0;
                            for (int jj = 0; jj < degree; jj++)
                            {
                                for (int l = 0; l < degree; l++)
                                {
                                    V_mik += Gamma(jj, ii).get_very_long() * Gamma(l,k).get_very_long() * AlgebraicNumber_in_O_pO::W_mult()(m,jl);
                                    jl++;
                                }
                            }
                            V(m,ik) = VeryLongModular(V_mik.numerator()) / VeryLongModular(V_mik.denominator());
                            ik++;
                        }
                    }
                }

                if (debug) std::cout << "V = " << std::endl << V;
                // To get the multiplication table in terms of the columns of Gamma
                // and the result of the columns of B1, calculate the inverse
                // image matrix of V, i.e. Y such that V = B1 * Y
                // (remember that B1 contains Gamma as a sub-matrix, and spans the
                // whole of O/pO)

                Matrix<VeryLongModular> Y = inverse_image_matrix(B1, V);

                // check result of inverse_image_matrix
                Matrix<VeryLongModular> checkV = B1 * Y;
                if (checkV != V)
                {
                    std::cout << "Problem: inverse_image_matrix incorrect" << std::endl;
                }
                if (debug) std::cout << "Y = " << std::endl << Y;

                // we have V = B1 * Y so B1 is degree x degree, V is degree x f^2
                // and so Y is degree x f^2
                // Pull out the multiplication table of the columns of Gamma,
                // ignoring the components along the basis of H_

                Matrix<VeryLongModular> A(f, f*f);
                for (int jj = 0; jj < f; jj++)
                {
                    int ik = 0;
                    for (int ii = 0; ii < f; ii++)
                    {
                        for (int k = 0; k < f; k++)
                        {
                            A(jj, ik) = Y(r+jj, ik);
                            ik++;
                        }
                    }
                }
                if (debug) std::cout << "A = " << std::endl << A;
                SeparableAlgebraElement::set_multiplication_table(A);

                // A is multiplication table for the columns of Gamma:
                // gamma[i] * gamma[k] = Sum(j = 1 ... f) A(j,ik) * gamma[j] in O/H = (O/pO) / (H/pO)

                // Step 10. [Compute V = ker(phi)]
                // First compute the matrix for the map phi: x -> x^p - x
                // where x is in the separable algebra A
                // by computing gamma[i]^p - gamma[i]
                Matrix<VeryLongModular> M(f, f);

                for (int ii = 0; ii < f; ii++)
                {
                    std::vector<VeryLongModular> v;
                    v.resize(f);
                    for (int jj = 0; jj < f; jj++) v[jj] = zero;
                    v[ii] = one;
                    SeparableAlgebraElement gam(v);
                    SeparableAlgebraElement phigam = pospow(gam, p) - gam;
                    if (debug) std::cout << "gamma[" << ii << "] = " << gam << std::endl;
                    if (debug) std::cout << "p = " << p << ", gamma[" << ii << "]^" << p << ") = " << pospow(gam, p) << std::endl;
                    if (debug) std::cout << "p = " << p << ", phi(gamma[" << ii << "]) = " << phigam << std::endl;
                    for (int jj = 0; jj < f; jj++)
                    {
                        M(jj, ii) = phigam.basis()[jj];
                    }
                }
                if (debug) std::cout << "M = " << std::endl << M;

                Matrix<VeryLongModular> M1 = kernel(M);
                if (debug) std::cout << "M1 = " << std::endl << M1;

                // Step 11. [Do we have a field]
                if (debug) std::cout << "Step 11." << std::endl;
                if (M1.columns() <= 1)
                {
                    if (debug) std::cout << "We have a field" << std::endl;
                    std::vector<AlgebraicNumber> gamma;
                    gamma.push_back(AlgebraicNumber(p));
                    if (debug) std::cout << "H_basis = " << std::endl << H_basis;
                    for (size_t col = 0; col < H_basis.columns(); col++)
                    {
                        AlgebraicNumber gam = AlgebraicNumber(0L);
                        for (int row = 0; row < degree; row++)
                        {
                            gam = gam + H_basis(row,col).get_very_long() * omega[row];
                        }
                        gamma.push_back(gam);
                    }
                    for (int ii = 1; ii < degree; ii++)
                    {
                        gamma.push_back(p*omega[ii]);
                    }
                    if (debug)
                    {
                        for (size_t ii = 0; ii < gamma.size(); ii++)
                        {
                            std::cout << "gamma[" << ii << "] = " << gamma[ii] << std::endl;
                            std::cout << "N(gamma[" << ii << "]) = " << gamma[ii].norm() << std::endl;
                        }
                    }
                    PrimeIdealT pi(p, gamma);
                    pi.e_ = j;
                    pi.f_ = f;
                    if (debug)
                    {
                        std::cout << "New prime ideal pi = " << std::endl << pi << std::endl;
                    }
                    primeIdeal.push_back(std::pair<PrimeIdealT, int>(pi, j));
                    L.erase(L.begin());
                    c--;
                    // .. go to step 7
                    done8 = 1;
                }
                else
                {
                    // Step 12. [Find m(X)]
                    // M1 has more than one column
                    // Find first column of M1 not proportional to transpose(1,0,0,...0)
                    size_t col = 0;
                    bool done = false;
                    while (col < M1.columns() && !done)
                    {
                        size_t row = 1;
                        while (row < M1.rows() && M1(row,col).is_zero())
                        {
                            row++;
                        }
                        if (row < M1.rows()) done = true;
                        else col++;
                    }
                    if (col >= M1.columns())
                    {
                        std::cerr << "Problem: can't find column in M1: " << std::endl << M1;
                        return;
                    }
                    std::vector<VeryLongModular> v;
                    v.resize(M1.rows());
                    for (size_t row = 0; row < M1.rows(); row++) v[row] = M1(row,col);
                    SeparableAlgebraElement alph(v);
                    if (debug) std::cout << "alph.F_p_basis_ = " << std::endl;
                    if (debug)
                    {
                        for (size_t rr = 0; rr < v.size(); rr++)
                        {
                            std::cout << alph.basis()[rr] << std::endl;
                        }
                    }

                    // Now compute successive powers of alph to find its minimal polynomial
                    int power = 2;
                    Matrix<VeryLongModular> AA_(f, power);
                    AA_(0,0) = VeryLongModular(1L);
                    for (int row = 1; row < f; row++)
                    {
                        AA_(row, 0) = VeryLongModular(0L);
                    }
                    for (int row = 0; row < f; row++)
                    {
                        AA_(row, 1) = alph.basis()[row];
                    }
                    SeparableAlgebraElement alph_pow = alph;
                    if (debug) std::cout << "alpha^" << power - 1 << " = " << alph_pow << std::endl;
                    Matrix<VeryLongModular> kerAA_ = kernel(AA_);

                    while (kerAA_.columns() == 0)
                    {
                        alph_pow = alph_pow * alph;
                        power++;
                        if (debug) std::cout << "alpha^" << power - 1 << " = " << alph_pow << std::endl;
                        AA_.add_column();
                        for (int row = 0; row < f; row++)
                        {
                            AA_(row, power-1) = alph_pow.basis()[row];
                        }
                        kerAA_ = kernel(AA_);
                    }

                    if (debug) std::cout << "AA_ = " << std::endl << AA_;
                    if (debug) std::cout << "kerAA_ = " << std::endl << kerAA_;
                    // columns of kerAA now give coefficients of relations between powers of alph
                    // find the minimal polynomial by looking for the column with least highest
                    // non-zero row
                    size_t least_highest_row = kerAA_.rows();
                    size_t min_col = 0;
                    for (size_t col = 0; col < kerAA_.columns(); col++)
                    {
                        size_t row = kerAA_.rows() - 1;
                        while (kerAA_(row,col).is_zero()) --row;
                        if (row < least_highest_row)
                        {
                            least_highest_row = row;
                            min_col = col;
                        }
                    }

                    // Construct minimal polynomial
                    std::vector<VeryLong> coeff;
                    coeff.resize(least_highest_row+1);
                    for (size_t row = 0; row < least_highest_row + 1; row++) coeff[row] = kerAA_(row, min_col).get_very_long();
                    Polynomial<VeryLong> m(coeff);
                    if (debug) std::cout << "Minimal polynomial is " << m << std::endl;

                    // Step 13. [Factor m(X)]
                    std::vector<Polynomial<VeryLongModular > > factors;
                    factor_over_F_p<VeryLong, VeryLong, VeryLongModular>(m, p, factors);

                    // Step 14. [Split H]
                    int k = factors.size();

                    int r = H_basis.columns();
                    if (debug) std::cout << "Step 14: k = " << k << ", r = " << r << std::endl;

                    for (int s = 1; s <= k; s++)
                    {
                        // assume each polynomial is linear
                        if (factors[s-1].deg() != 1)
                        {
                            std::cerr << "Problem with final split: factors[" << s << "] = " << factors[s-1] << std::endl;
                        }
                        factors[s-1].make_monic();
                        if (debug) std::cout << "factor " << s << " = " << factors[s-1] << std::endl;
                        SeparableAlgebraElement beta = alph + factors[s-1].coefficient(0);
                        Matrix<VeryLongModular> M(degree, r + degree);
                        for (int col = 0; col < r; col++)
                        {
                            for (int row = 0; row < degree; row++)
                            {
                                M(row,col) = H_basis(row,col);
                            }
                        }

                        for (int ii = 0; ii < degree; ii++)
                        {
                            for (int m = 0; m < degree; m++)
                            {
                                Quotient<VeryLong> tmp(0L);
                                for (int jj = 0; jj < degree; jj++)
                                {
                                    for (int l = 0; l < f; l++)
                                    {
                                        tmp += beta.basis()[l].get_very_long() *
                                               Gamma(jj, l).get_very_long() *
                                               AlgebraicNumber_in_O_pO::W_mult()(m, ii*degree + jj);
                                    }
                                }
                                M(m,ii + r) = VeryLongModular(tmp.numerator()) / VeryLongModular(tmp.denominator());
                            }
                        }

                        if (debug) std::cout << "M = " << std::endl << M;
                        Matrix<VeryLongModular> H_s = image(M);
                        if (debug) std::cout << "H_s = " << std::endl << H_s;
                        if (H_s.columns() < static_cast<size_t>(degree))
                        {
                            L.push_back(Ideal_mod_pO(H_s));
                            c++;
                        }
                    }

                    // Step 15. [Update list]
                    // Remove H from list
                    L.erase(L.begin());
                    c--;
                    // .. go to step 8
                }
            } // done8
        } // main loop

    }

    // Algorithm 4.8.17 (Valuation at a Prime Ideal)
    // p is assumed to be a prime ideal
    static int padicValuation(const Z& p, const PrimeIdealT& pi, const Ideal& I)
    {
        const NumberField& nf = AlgebraicNumber::nf();
        int degree = AlgebraicNumber::degree();
        // (N, d) is HNF of the maximal order Z_k
        VeryLong d1 = I.denominator();
        //   cout << "PrimeIdeal::padicValuation : d1 = " << d1 << endl;
        Matrix<Quotient<VeryLong> > Mq(degree, degree);

        // degree of the number field

        Matrix<VeryLong> hnf_basis = I.hnf_basis();
        for (int i = 0; i < degree; i++)
        {
            for (int j = 0; j < degree; j++)
            {
                Mq(i,j) = Quotient<VeryLong>(hnf_basis(i,j), d1);
            }
        }

        // Also need the integral basis as a list of AlgebraicNumber objects
        const std::vector<AlgebraicNumber>& omega = AlgebraicNumber::integralBasis();
        VeryLongModular::set_default_modulus(p);

        const AlgebraicNumber& beta = pi.beta_;

        // Step 3. [Compute N(I)]
        Matrix<Quotient<VeryLong> > Aq = nf.winv() * Mq;
        //   cout << "Step 3: Aq = " << endl << Aq;
        // check if we have an integral ideal
        VeryLong lcm(1L);
        for (int i = 0; i < degree; i++)
        {
            for (int j = 0; j < degree; j++)
            {
                lcm = lcm * Aq(i,j).denominator() / gcd(lcm, Aq(i,j).denominator());
            }
        }
        const VeryLong one(1L);
        const VeryLong zero(0L);
        if (lcm != one)
        {
            //      cout << "not an integral ideal\n" << endl;
            Ideal J(I);
            J *= lcm;
            int val = PrimeIdealT::padicValuation(p, pi, J);
            int e = pi.ramificationIndex();
            int v = 0;
            VeryLong x = lcm;
            while (x % p == zero)
            {
                x /= p;
                v++;
            }
            val -= e * v;
            return val;
        }

        Matrix<VeryLong> AA(degree, degree);
        VeryLong P(1L);
        for (int i = 0; i < degree; i++)
        {
            for (int j = 0; j < degree; j++)
            {
                AA(i,j) = Aq(i,j).numerator();
            }
            P *= AA(i,i);
        }
        // Check
        if (Quotient<VeryLong>(P) != I.norm())
        {
            std::cout << "Problem, P = " << P << ", N(I) = " << I.norm() << std::endl;
        }

        if (P % p != 0L) return 0;

        int nu = 0;

        int done4 = 0;
        while (!done4)
        {
            // Step 4. [Multiply]
            for (int c = 0; c < degree; c++)
            {
                AlgebraicNumber a(VeryLong(0L));
                for (int r = 0; r < degree; r++)
                {
                    a += AA(r,c) * omega[r];
                }
                a *= beta;
                // Now find coefficients of a in terms of omega
                Matrix<Quotient<VeryLong> > at(degree, 1);
                for (int r = 0; r < degree; r++)
                {
                    at(r,0) = a.coefficients()[r];
                }
                at = nf.winv() * at;

                for (int r = 0; r < degree; r++)
                {
                    if (at(r,0).denominator() != VeryLong(1L))
                    {
                        std::cout << "Problem: not an integer : at[" << r << "] = " << at(r,0) << std::endl;
                    }
                    AA(r,c) = at(r,0).numerator();
                }
            }

            // Step 5. [Simple test]
            if (AA(degree-1, degree-1) % p != 0L)
            {
                return nu;
            }
            if (AlgebraicNumber::index() % p != 0L)
            {
                nu++;
                AA = AA / p;
            }
            else
            {
                for (int i = 0; i < degree; i++)
                {
                    for (int j = 0; j < degree; j++)
                    {
                        if (AA(i,j) % p != 0L) return nu;
                    }
                }
                AA = AA / p;
                nu++;
            }
        }

        return 0;
    }

    // Algorithm 4.8.17 (Valuation at a Prime Ideal)
    // p is assumed to be a prime ideal
    // matrix A representing an integral ideal wrt integral basis
    static int padicValuation(const Z& p, const PrimeIdealT& pi, const Matrix<VeryLong>& A)
    {
        const VeryLong zero(0L);
        const VeryLong one(1L);
        const NumberField& nf = AlgebraicNumber::nf();
        int degree = AlgebraicNumber::degree();

        int divisibleByp = 0;
        for (int i = 0; !divisibleByp && i < degree; i++)
        {
            if (A(i,i) % p == zero) divisibleByp = 1;
        }

        if (!divisibleByp) return 0;

        pi.beta_omega();

        VeryLongModular::set_default_modulus(p);

        static Matrix<VeryLong> AA(degree, degree);
        AA = A;

        int nu = 0;

        int done4 = 0;
        static Matrix<Quotient<VeryLong> > at(degree, degree);
        static Quotient<VeryLong> ac;
        static Quotient<VeryLong> q(zero);
        while (!done4)
        {
            // Step 4. [Multiply]
            for (int c = 0; c < degree; c++)
            {
                for (int r = 0; r < degree; r++)
                {
                    ac = zero;
                    for (int rr = 0; rr < degree; rr++)
                    {
                        q = pi.beta_omega_[rr].coefficients()[r];
                        q *= AA(rr,c);
                        ac += q;
                    }
                    at(r,c) = ac;
                }
            }

            // Now find coefficients in terms of omega
            at = nf.winv() * at;

            for (int c = 0; c < degree; c++)
            {
                for (int r = 0; r < degree; r++)
                {
                    AA(r,c) = at(r,c).numerator();
                }
            }

            // Step 5. [Simple test]
            //AA = HNF(AA); // correction in Errata
            if (AA(degree-1, degree-1) % p != 0L)
            {
                done4 = 1;
            }
            else if (AlgebraicNumber::index() % p != 0L)
            {
                nu++;
                //AA = AA / p;
                AA /= p;
            }
            else
            {
                for (int i = 0; i < degree; i++)
                {
                    for (int j = 0; j < degree; j++)
                    {
                        if (AA(i,j) % p != 0L) return nu;
                    }
                }
                //         AA = AA / p;
                AA /= p;
                nu++;
            }
        }
        return nu;
    }

    // Algorithm 4.8.17 (Valuation at a Prime Ideal)
    // p is assumed to be a prime ideal
    // matrix A representing an integral ideal wrt integral basis
    // lcm gives corresponding denominator to give a possibly fractional ideal
    static int padicValuation(const VeryLong& p, const PrimeIdealT& pi,
                              const VeryLong& lcm, const Matrix<VeryLong>& A)
    {
        int nu = padicValuation(p, pi, A);
        const VeryLong one(1L);
        if (lcm != one)
        {
            int e = pi.ramificationIndex();
            int v = 0;
            VeryLong x = lcm;
            const VeryLong zero(0L);
            while (x % p == zero)
            {
                x /= p;
                v++;
            }
            nu -= e * v;
        }

        return nu;
    }

    friend ostream& operator<<(ostream& os, const PrimeIdealT& PI)
    {
        os << "Prime Ideal : denominator = " << PI.denominator_ << ", HNF basis = " << std::endl << PI.hnf_basis_;
        os << "p_ = " << PI.p_ << ", alpha_ = " << PI.alpha_ << std::endl;
        os << ", e_ = " << PI.e_ << ", f_ = " << PI.f_ << ", beta_ = " << PI.beta_ << std::endl;
        os << ", norm = " << PI.norm() << std::endl;

        return os;
    }


    PrimeIdealT& operator=(const PrimeIdealT& pi)
    {
        if (&pi != this)
        {
            hnf_basis_ = pi.hnf_basis_;
            denominator_ = pi.denominator_;
            isPrincipal_ = pi.isPrincipal_;
            p_ = pi.p_;
            alpha_ = pi.alpha_;
            e_ = pi.e_;
            f_ = pi.f_;
            beta_ = pi.beta_;
            beta_omega_ = pi.beta_omega_;
            beta_computed_ = pi.beta_computed_;
        }
        return *this;
    }

    friend Ideal& operator*=(Ideal& I, const PrimeIdealT& pi)
    {
        // More efficient multiplication of an ideal by a prime ideal, using
        // the 2-element representation of the prime ideal
        int d = pi.hnf_basis_.rows();
        Matrix<VeryLong> basis(d, 2*d);
        Matrix<Quotient<VeryLong> > qbasis(d, 2*d);

        int k = 0;
        AlgebraicNumber p_an(pi.p_);
        VeryLong lcm(1L);
        for (int i = 0; i < d; i++)
        {
            AlgebraicNumber u(I.hnf_basis(), I.denominator(), i);
            AlgebraicNumber v = u * p_an;
            for (int l = 0; l < d; l++)
            {
                qbasis(l,k) = v.coefficient(l);
                lcm = lcm * v.coefficient(l).denominator() / gcd(lcm, v.coefficient(l).denominator());
            }
            k++;
            v = u * pi.alpha_;
            for (int l = 0; l < d; l++)
            {
                qbasis(l,k) = v.coefficient(l);
                lcm = lcm * v.coefficient(l).denominator() / gcd(lcm, v.coefficient(l).denominator());
            }
            k++;
        }

        for (int i = 0; i < d; i++)
        {
            for (int j = 0; j < 2 * d; j++)
            {
                basis(i,j) = qbasis(i,j).numerator() * (lcm / qbasis(i,j).denominator());
            }
        }
        I = Ideal(HNF1(basis), lcm);

        return I;
    }

    friend Ideal& operator*=(Ideal& I, const PrimeIdealRep& pi);

    VeryLong p() const
    {
        return p_;
    }
    int ramificationIndex() const
    {
        return e_;
    }
    int residualIndex() const
    {
        return f_;
    }
    const std::vector<AlgebraicNumber>& beta_omega() const
    {
        if (!beta_computed_) computeBeta();
        return beta_omega_;
    }
    AlgebraicNumber alpha() const
    {
        return alpha_;
    }

private:
    void computeBeta() const
    {
        const NumberField& nf = AlgebraicNumber::nf();
        // degree of the number field
        int degree = AlgebraicNumber::degree();
        // (N, d) is HNF of the maximal order Z_k

        // Also need the integral basis as a list of AlgebraicNumber objects
        const std::vector<AlgebraicNumber>& omega = AlgebraicNumber::integralBasis();
        VeryLongModular::set_default_modulus(p_);

        std::vector<AlgebraicNumber> gamma;
        gamma.resize(degree);
        for (int i = 0; i < degree; i++)
        {
            gamma[i] = AlgebraicNumber(hnf_basis_, denominator_, i);
        }

        // Step 1. [Compute structure constants]
        Matrix<VeryLongModular> A(degree*degree, degree);
        const Matrix<Quotient<VeryLong> >& Nq = nf.w();
        for (int i = 0; i < degree; i++)
        {
            for (int j = 0; j < degree; j++)
            {
                AlgebraicNumber tmpij = omega[i] * gamma[j];
                //         cout << "omega[" << i << "] * gamma[" << j << "] = " << tmpij << ", N(tmpij) = " << tmpij.norm() << endl;
                // Bij is coefficients of tmpij in terms of alpha
                std::vector<Quotient<VeryLong> > Aij;
                std::vector<Quotient<VeryLong> > Bij;
                Aij.resize(degree);
                Bij.resize(degree);
                for (int r = 0; r < degree; r++)
                {
                    Bij[r] = tmpij.coefficients()[r];
                }
                for (int k = degree; k > 0; --k)
                {
                    Quotient<VeryLong> S(0L);
                    for (int l = k + 1; l <= degree; l++)
                    {
                        Quotient<VeryLong> tmp = Aij[l-1] * Nq(k-1,l-1);
                        S += tmp;
                    }
                    Aij[k-1] = (Bij[k-1] - S) / Nq(k-1,k-1);
                    A(j + (k-1)*degree, i) = Aij[k-1].numerator() % p_;
                    if (Aij[k-1].denominator() != VeryLong(1L))
                    {
                        std::cout << "Problem: Aijk not integer = " << Aij[k-1] << std::endl;
                    }
                }
                // check
                AlgebraicNumber tmpij1(0L);
                for (int k = 0; k < degree; k++)
                {
                    tmpij1 = tmpij1 + Aij[k] * omega[k];
                }
                if (tmpij1 != tmpij)
                {
                    std::cout << "Problem: tmpij1 = " << tmpij1 << std::endl;
                    std::cout << "         tmpij  = " << tmpij << std::endl;
                }
            }
        }

        // Step 2. [Compute beta]
        //cout << "A = " << endl << A;
        Matrix<VeryLongModular> Beta = kernel(A);
        if (Beta.columns() < 1)
        {
            std::cout << "Problem: kernel of A = " << std::endl << A;
            std::cout << " is empty" << std::endl;
            return;
        }
        //cout << "ker(A) = " << endl << Beta;

        // Explicitly calculate beta in the kernel of A
        AlgebraicNumber beta(VeryLong(0L));
        for (int j = 0; j < degree; j++)
        {
            beta = beta + Beta(j,0).get_very_long() * omega[j];
        }
        beta_ = beta;
        beta_omega_.resize(degree);
        for (int j = 0; j < degree; j++)
        {
            beta_omega_[j] = beta_ * omega[j];
        }
        beta_computed_ = 1;
    }

    // Algorithm 4.7.10 (Two-Element Representation of a Prime Ideal)
    void computeTwoElementRep(const std::vector<AlgebraicNumber>& gamma)
    {
        const VeryLong zero(0L);
        //   cout << "PrimeIdeal::computeTwoElementRep" << endl;
        VeryLong p_f = norm().numerator();
        int k = gamma.size();
        // Step 1. [Initialize]
        VeryLong R(1L);

        // First of all check the gamma, to see if one of them will do
        for (int i = 1; i < k; i++)
        {
            alpha_ = gamma[i];
            //      cout << "Trying alpha_ = " << alpha_ << ", N(alpha_) = " << alpha_.norm() << endl;
            Quotient<VeryLong> n = alpha_.norm() / Quotient<VeryLong>(p_f);
            if (n.numerator() % p_ != zero) return;
            AlgebraicNumber a = alpha_ + AlgebraicNumber(p_);
            //      cout << "a.norm() = " << a.norm() << endl;
            n = a.norm() / Quotient<VeryLong>(p_f);
            //      cout << "n = " << n << endl;
            if (n.numerator() % p_ != zero) return;
        }

        while (1)
        {
            // Step 2. [Set coefficients]
            std::vector<VeryLong> lambda;
            lambda.resize(k);
            for (int i = 2; i <= k; i++) lambda[i-1] = R;

            int done3 = 0;
            while (!done3)
            {
                // Step 3. [Compute alpha_ and check]
                alpha_ = AlgebraicNumber(zero);
                for (int i = 2; i <= k; i++)
                {
                    alpha_ = alpha_ + lambda[i-1] * gamma[i-1];
                }
                //         cout << "Trying alpha_ = " << alpha_ << ", N(alpha_) = " << alpha_.norm() << endl;

                Quotient<VeryLong> n = alpha_.norm() / Quotient<VeryLong>(p_f);
                //         cout << "n = " << n << endl;
                if (n.numerator() % p_ != zero) return;

                AlgebraicNumber a = alpha_ + AlgebraicNumber(p_);
                //         cout << "a.norm() = " << a.norm() << endl;
                n = a.norm() / Quotient<VeryLong>(p_f);
                //         cout << "n = " << n << endl;
                if (n.numerator() % p_ != zero) return;

                // Step 4. [Decrease coefficients]
                int j = k;
                while (j >= 1 && lambda[j-1] == -R) j--;

                lambda[j-1] -= 1L;

                for (int i = j+1; i <= k; i++) lambda[i-1] = R;

                // Step 5. [Search for first non-zero]
                j = 1;
                while (j <= k && lambda[j-1].is_zero()) j++;
                if (j > k)
                {
                    R += 1L;
                    done3 = 1;
                }
            }
        }
    }

    int e_;
    int f_;
    mutable int beta_computed_;
    VeryLong p_;
    AlgebraicNumber alpha_;
    mutable AlgebraicNumber beta_;
    mutable std::vector<AlgebraicNumber> beta_omega_;
};

typedef PrimeIdealT<VeryLong> PrimeIdeal;

Ideal& operator*=(Ideal& I, const PrimeIdealRep& pi);

Ideal& operator*=(Ideal& I, const std::vector<PrimeIdealRep*>& pirs);

class PrimeIdealRep
{
public:
    PrimeIdealRep(long int p, long int r)
        : p_(p), r_(r), pi_(0), alpha_(0), special_(0), norm_(0)
    {}
    PrimeIdealRep(PrimeIdeal* pi)
        : p_(0L), r_(0L), pi_(pi), alpha_(0), special_(1), norm_(0)
    {}
    ~PrimeIdealRep()
    {
        delete norm_;
        delete pi_;
    }
    PrimeIdeal* getPrimeIdeal() const
    {
        if (!pi_)
        {
            pi_ = new PrimeIdeal(VeryLong(p_), VeryLong(r_));
        }
        return pi_;
    }
    AlgebraicNumber* getAlpha() const
    {
        if (!alpha_)
        {
            alpha_ = new AlgebraicNumber(-VeryLong(r_)*AlgebraicNumber::c_d(), -AlgebraicNumber::c_d());
        }

        return alpha_;
    }

    friend ostream& operator<<(ostream& os, const PrimeIdealRep& PIR);
    VeryLong norm() const
    {
        if (norm_) return *norm_;

        if (pi_)
        {
            norm_ = new VeryLong(pi_->norm().numerator());
            return *norm_;
        }
        else return p_;
    }
    void clearPrimeIdeal()
    {
        if (special_) return;
        delete pi_;
        delete alpha_;
        pi_ = 0;
        alpha_ = 0;
    }

    bool isSpecial() const
    {
        return (special_ == 1);
    }

    long int p() const
    {
        return p_;
    }

    long int r() const
    {
        return r_;
    }
private:
    long int p_;
    long int r_;
    mutable PrimeIdeal* pi_;
    mutable AlgebraicNumber* alpha_;
    char special_;
    mutable VeryLong* norm_;
};

#endif
