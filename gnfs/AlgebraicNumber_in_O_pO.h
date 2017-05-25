#ifndef ALGEBRAICNUMBER_IN_O_PO_H
#define ALGEBRAICNUMBER_IN_O_PO_H
#include <string>
#include <sstream>
template <class INTEGER, class INTEGER2, class MODULAR_INTEGER> class AlgebraicNumber_in_O_pO_
{
public:
    AlgebraicNumber_in_O_pO_()
    {
        int d = AlgebraicNumber::degree();
        if (d >= MAX_DEGREE)
        {
            std::ostringstream oss;
            oss << "degree of number field (" << d << ") too big, must be < " << MAX_DEGREE;
            throw std::string(oss.str());
        }
        const MODULAR_INTEGER zero(0L);
        for (int i = 0; i < d; i++) Fp_basis_[i] = zero;
    }
    AlgebraicNumber_in_O_pO_(long int li)
    {
        int d = AlgebraicNumber::degree();
        if (d >= MAX_DEGREE)
        {
            std::ostringstream oss;
            oss << "degree of number field (" << d << ") too big, must be < " << MAX_DEGREE;
            throw std::string(oss.str());
        }
        const MODULAR_INTEGER zero(0L);
        for (int i = 0; i < d; i++) Fp_basis_[i] = zero;
        Fp_basis_[0] = MODULAR_INTEGER(li);
    }
    AlgebraicNumber_in_O_pO_(int dummy, int basisElement)
    {
        int d = AlgebraicNumber::degree();
        if (d >= MAX_DEGREE)
        {
            std::ostringstream oss;
            oss << "degree of number field (" << d << ") too big, must be < " << MAX_DEGREE;
            throw std::string(oss.str());
        }
        if (basisElement < 0 || basisElement >= d)
        {
            std::ostringstream oss;
            oss << "basisElement (" << basisElement << ") out of range, must be in [0," << d - 1 << "]";
            throw std::string(oss.str());
        }
        const MODULAR_INTEGER zero(0L);
        const MODULAR_INTEGER one(1L);
        for (int i = 0; i < d; i++) Fp_basis_[i] = zero;
        Fp_basis_[basisElement] = one;
    }

    AlgebraicNumber_in_O_pO_(const std::vector<MODULAR_INTEGER>& c)
    {
        int d = AlgebraicNumber::degree();
        if (d >= MAX_DEGREE)
        {
            std::ostringstream oss;
            oss << "degree of number field (" << d << ") too big, must be < " << MAX_DEGREE;
            throw std::string(oss.str());
        }
        if (static_cast<int>(c.size()) > d)
        {
            std::ostringstream oss;
            oss << "too many elements supplied, must be <= " << d;
            throw std::string(oss.str());
        }
        for (size_t i = 0; i < c.size(); ++i)
        {
            Fp_basis_[i] = c[i];
        }
    }

    AlgebraicNumber_in_O_pO_(const AlgebraicNumber& a)
    {
        const std::vector<Quotient<INTEGER2> >& c1 = a.ib_coefficients();
        for (size_t i = 0; i < c1.size(); i++)
        {
            // c1 gives coefficients of a in terms of omega
            Fp_basis_[i] = MODULAR_INTEGER(c1[i].numerator() % p_) / MODULAR_INTEGER(c1[i].denominator() % p_);
        }
    }

    AlgebraicNumber_in_O_pO_(long long int a, long int b)
    {
        if (!optimisation_ok_)
        {
            throw std::string("Problem: can't use optimisation to create AlgebraicNumber_in_O_pO_");
        }
        // corresponding to a - b alpha
        Fp_basis_[0] = MODULAR_INTEGER(a) - MODULAR_INTEGER(b) * w01_;
        Fp_basis_[1] = MODULAR_INTEGER(-b) * w11_;

        for (size_t i = 2; i < static_cast<size_t>(AlgebraicNumber::degree()); ++i)
        {
            Fp_basis_[i] = MODULAR_INTEGER(0L);
        }
    }

    AlgebraicNumber_in_O_pO_& operator=(const AlgebraicNumber& a)
    {
        *this = AlgebraicNumber_in_O_pO_(a);
        return *this;
    }

    AlgebraicNumber_in_O_pO_(const AlgebraicNumber_in_O_pO_& a)
    {
        int d = AlgebraicNumber::degree();
        for (int i = 0; i < d; i++) Fp_basis_[i] = a.Fp_basis_[i];
    }
    AlgebraicNumber_in_O_pO_& operator=(const AlgebraicNumber_in_O_pO_& a)
    {
        if (this != &a)
        {
            for (int i = 0; i < AlgebraicNumber::degree(); i++) Fp_basis_[i] = a.Fp_basis_[i];
        }
        return *this;
    }

    ~AlgebraicNumber_in_O_pO_() {}

    AlgebraicNumber_in_O_pO_ operator*(const AlgebraicNumber_in_O_pO_& b)
    {
        AlgebraicNumber_in_O_pO_ x;
        int degree = AlgebraicNumber::degree();
        const MODULAR_INTEGER zero(0L);
        for (int k = 0; k < degree; k++)
        {
            //typename std::vector<MODULAR_INTEGER>::const_iterator ij_iter = M_[k].begin();
            int ij = 0;
            x.Fp_basis_[k] = zero;
            for (int i = 0; i < degree; i++)
            {
                for (int j = 0; j < degree; j++)
                {
                    x.Fp_basis_[k] += Fp_basis_[i] * b.Fp_basis_[j] * M_(k, ij);
                    ++ij;
                }
            }
        }
        return x;
    }
    AlgebraicNumber_in_O_pO_& operator*=(const AlgebraicNumber_in_O_pO_& b)
    {
        int degree = AlgebraicNumber::degree();
        MODULAR_INTEGER zero(0L);
        static MODULAR_INTEGER* tmp = 0;
        static int tmp_size = 0;
        if (!tmp || tmp_size != degree)
        {
            tmp = new MODULAR_INTEGER[degree];
            tmp_size = degree;
        }

        static MODULAR_INTEGER x;
        static MODULAR_INTEGER y;
        for (int k = 0; k < degree; k++)
        {
            int ij = 0;
            tmp[k] = zero;
#if 1
            for (int i = 0; i < degree; ++i)
            {
                ij = i * (degree + 1);
                tmp[k].add_product(Fp_basis_[i], b.Fp_basis_[i], M_(k, ij));
            }
            for (int i = 0; i < degree; ++i)
            {
                for (int j = i + 1; j < degree; ++j)
                {
                    ij = i * degree + j;
                    x = Fp_basis_[i];
                    x *= b.Fp_basis_[j];
                    y = Fp_basis_[j];
                    y *= b.Fp_basis_[i];
                    x += y;
                    x *= M_(k, ij);
                    tmp[k] += x;
                }
            }
#else
            for (int i = 0; i < degree; i++)
            {
                for (int j = 0; j < degree; j++)
                {
#if 0
                    x = Fp_basis_[i];
                    x *= b.Fp_basis_[j];
                    x *= M_(k, ij);
                    tmp[k] += x;
#else
                    tmp[k].add_product(Fp_basis_[i], b.Fp_basis_[j], M_(k, ij));
#endif
                    ++ij; // ij = i * degree + j
                }
            }
#endif
        }

        for (int k = 0; k < degree; k++)
        {
            Fp_basis_[k] = tmp[k];
        }
        return *this;
    }
    AlgebraicNumber_in_O_pO_& operator*=(const AlgebraicNumber& b)
    {
        int degree = AlgebraicNumber::degree();
        MODULAR_INTEGER zero(0L);
        static MODULAR_INTEGER* tmp = 0;
        static int tmp_size = 0;
        if (!tmp || tmp_size != degree)
        {
            tmp = new MODULAR_INTEGER[degree];
            tmp_size = degree;
        }

        MODULAR_INTEGER x;
        const std::vector<Quotient<INTEGER2> >& c = b.ib_coefficients();
        for (int k = 0; k < degree; k++)
        {
            int ij = 0;
            tmp[k] = zero;
            for (int i = 0; i < degree; i++)
            {
                for (int j = 0; j < degree; j++)
                {
                    x = Fp_basis_[i];
                    x *= c[j].numerator() % p_;
                    x /= c[j].denominator() % p_;
                    x *= M_(k, ij);
                    tmp[k] += x;
                    ++ij;
                }
            }
        }

        for (int k = 0; k < degree; k++)
        {
            Fp_basis_[k] = tmp[k];
        }
        return *this;
    }
    AlgebraicNumber_in_O_pO_ operator/(const AlgebraicNumber_in_O_pO_& b)
    {
        int d = AlgebraicNumber::degree();
        AlgebraicNumber_in_O_pO_ x;
        Matrix<MODULAR_INTEGER> V(d, d);
        const MODULAR_INTEGER zero(0L);

        for (int k = 0; k < d; k++)
        {
            int ij = 0;
            for (int j = 0; j < d; j++)
            {
                V(k,j) = zero;
                for (int i = 0; i < d; i++)
                {
                    V(k,j) += b.Fp_basis_[i] * M_(k,ij);
                    ij++; // ij = j*d + i
                }
            }
        }
        Matrix<MODULAR_INTEGER> Vinv(d, d);
        invert(V, Vinv);
        for (int j = 0; j < d; j++)
        {
            x.Fp_basis_[j] = zero;
            for (int k = 0; k < d; k++)
            {
                x.Fp_basis_[j] += Vinv(j,k) * Fp_basis_[k];
            }
        }

        return x;
    }
    int operator==(const AlgebraicNumber_in_O_pO_& a)
    {
        if (this == &a) return 1;
        for (int i = 0; i < AlgebraicNumber::degree(); i++)
        {
            if (Fp_basis_[i] != a.Fp_basis_[i]) return 0;
        }
        return 1;
    }
    int operator!=(const AlgebraicNumber_in_O_pO_& a)
    {
        return !(*this == a);
    }

    static void set_basis(const Matrix<Quotient<INTEGER2> >& W,
                          const std::vector<AlgebraicNumber>& omega)
    {
        size_t d = AlgebraicNumber::degree();
        if (M_.rows() != d || M_.columns() != d*d) M_.set_size(d, d * d);
        // Calculate multiplication table for W
        // W gives coefficients for integral basis omega, in
        // terms of alpha (root of non-monic polynomial)
        static Matrix<Quotient<VeryLong> > V(1, 1);
        if (V.rows() != d || V.columns() != d*d) V.set_size(d, d * d);
        int ij = 0;
        AlgebraicNumber x;
        for (auto& i: omega)
        {
            for (auto& j: omega)
            {
                x = i * j;
                // get coefficients in terms of theta
                int m = 0;
                for (auto& x_co: x.coefficients())
                {
                    V(m,ij) = x_co;
                    m++;
                }
                ij++;
            }
        }
        // W_mult_ is inverse image matrix
        W_mult_ = inverse_image_matrix(W, V);
    }

    static void set_basis(const INTEGER& p)
    {
        p_ = p;
        MODULAR_INTEGER::set_default_modulus(p_);
        long int d = AlgebraicNumber::degree();
        for (int i = 0; i < d; i++)
        {
            for (int j = 0; j < d*d; j++)
            {
                M_(i, j) = W_mult()(i, j).numerator() % p_;
            }
        }
        MODULAR_INTEGER w01n = AlgebraicNumber::nf().winv()(0,1).numerator() % p_;
        MODULAR_INTEGER w01d = AlgebraicNumber::nf().winv()(0,1).denominator() % p_;
        MODULAR_INTEGER w11n = AlgebraicNumber::nf().winv()(1,1).numerator() % p_;
        MODULAR_INTEGER w11d = AlgebraicNumber::nf().winv()(1,1).denominator() % p_;
        // If p_ doesn't divide the denominator if winv()(0,1) or winv()(1,1) then
        // we can pre-calculate some numbers (mod p_) which are used when constructing
        // a - b alpha (mod p_) in O/pO
        // In practice, this is only needed when p_ is an inert prime, which will
        // always satisfy this condition
        optimisation_ok_ = false;
        if (w01d != MODULAR_INTEGER(0L) && w11d != MODULAR_INTEGER(0L))
        {
            w01_ = w01n / w01d;
            w11_ = w11n / w11d;
            optimisation_ok_ = true;
        }
    }

    static void set_basis(const Matrix<Quotient<INTEGER2> >& W,
                          const std::vector<AlgebraicNumber>& omega,
                          const INTEGER& p)
    {
        set_basis(W, omega);
        set_basis(p);
    }

    static Matrix<MODULAR_INTEGER> make_beta(const INTEGER2& q)
    {
        //std::cerr << "make_beta : q = <" << q << ">, modulus = <" << MODULAR_INTEGER::get_default_modulus() << ">" << std::endl;
        int d = AlgebraicNumber::degree();
        Matrix<MODULAR_INTEGER> Ap(d, d);
        for (int i = 0; i < d; i++)
        {
            AlgebraicNumber_in_O_pO_ tmp2(0, i);
            tmp2 = pow<AlgebraicNumber_in_O_pO_>(tmp2, q);
            for (int j = 0; j < d; j++)
            {
                // We want to find the kernel of the transpose,
                // so transpose now
                Ap(j, i) = tmp2.basis()[j];
            }
        }
        // We have Ap, so find kernel of A over F_p
        return kernel(Ap);
    }

    friend ostream& operator<<(ostream& os, const AlgebraicNumber_in_O_pO_ & a)
    {
        int d = AlgebraicNumber::degree();
        const MODULAR_INTEGER zero(0L);
        int first = 1;
        for (int i = 0; i < d; i++)
        {
            MODULAR_INTEGER value = a.Fp_basis_[i];
            if (value != zero)
            {
                if (first) first = 0;
                else os << " + ";
                os << value;
                os << " omega_" << i + 1;
            }
        }

        return os;
    }
    MODULAR_INTEGER coefficient(int i) const
    {
        int d = AlgebraicNumber::degree();
        if (i < 0 || i >= d)
        {
            std::ostringstream oss;
            oss << "index (" << i << ") out of range, must be in [0," << d - 1 << "]";
            throw std::string(oss.str());
        }
        return Fp_basis_[i];
    }

    const MODULAR_INTEGER* basis() const
    {
        return Fp_basis_;
    }
    static const Matrix<Quotient<INTEGER2> >& W_mult()
    {
        return W_mult_;
    }

private:
    // W_mult_ is a dxd^2 matrix giving the coefficients of omega[i] x omega[j] in terms of omega.
    // ??? I'm not sure why this isn't in AlgebraicNumber ???
    // Note that this should really be a Matrix<INTEGER2>
    static Matrix<Quotient<INTEGER2> > W_mult_;
    // M_ is W_mult_ reduced modulo p_
    static Matrix<MODULAR_INTEGER> M_;
    static INTEGER p_;
    static MODULAR_INTEGER w01_;
    static MODULAR_INTEGER w11_;
    static bool optimisation_ok_;

    enum { MAX_DEGREE = NumberField::MAX_DEGREE };
    MODULAR_INTEGER Fp_basis_[MAX_DEGREE];
};
typedef AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular> AlgebraicNumber_in_O_pO;
typedef AlgebraicNumber_in_O_pO_<long int, VeryLong, LongModular> AlgebraicNumber_in_O_pO_1;

#endif
