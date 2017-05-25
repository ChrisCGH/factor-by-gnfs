#ifndef ALGEBRAICNUMBER_H
#define ALGEBRAICNUMBER_H

#include "Polynomial.inl"
class AlgebraicNumber
{
public:
    friend class AlgebraicNumberTest;
    AlgebraicNumber();
    explicit AlgebraicNumber(const VeryLong& v);
    AlgebraicNumber(const VeryLong& a, const VeryLong& b);
    AlgebraicNumber(std::vector<Quotient<VeryLong > >& coeff);
    AlgebraicNumber(const AlgebraicNumber& a);
    AlgebraicNumber(const Matrix<Quotient<VeryLong> >& M, int column);
    AlgebraicNumber(const Matrix<VeryLong>& M, const VeryLong& denominator, int column);
    AlgebraicNumber(const Polynomial<VeryLong>& f, const AlgebraicNumber& a);
    AlgebraicNumber(const Polynomial<Quotient<VeryLong> >& f);
    AlgebraicNumber(const Quotient<VeryLong>& f);
    ~AlgebraicNumber();
    AlgebraicNumber& operator=(const AlgebraicNumber& a);
    int operator==(const AlgebraicNumber& a) const;
    int operator!=(const AlgebraicNumber& a) const;
    int operator<(const AlgebraicNumber& a) const
    {
        return 0;
    }
    AlgebraicNumber operator-() const
    {
        AlgebraicNumber out(VeryLong(0L));
        return out - *this;
    }

    friend AlgebraicNumber operator*(const AlgebraicNumber& a1,
                                     const AlgebraicNumber& a2);
    AlgebraicNumber& operator*=(long int p);
    AlgebraicNumber& operator*=(const VeryLong& v);
    AlgebraicNumber& operator*=(const AlgebraicNumber& a);
    AlgebraicNumber& multiply(const VeryLong& a, const VeryLong& b);
    AlgebraicNumber& multiply_by_alpha_minus_r(long int r);

    friend AlgebraicNumber operator*(const AlgebraicNumber& a,
                                     const Quotient<VeryLong>& x);

    friend AlgebraicNumber operator*(const Quotient<VeryLong>& x,
                                     const AlgebraicNumber& a);

    friend AlgebraicNumber operator+(const AlgebraicNumber& a1,
                                     const AlgebraicNumber& a2);
    AlgebraicNumber& operator+=(const AlgebraicNumber& a);

    friend AlgebraicNumber operator/(const AlgebraicNumber& a1,
                                     const AlgebraicNumber& a2);

    friend AlgebraicNumber operator-(const AlgebraicNumber& a1,
                                     const AlgebraicNumber& a2);
    void display()
    {
        const std::vector<Quotient<VeryLong> >& c = ib_coefficients();
        for (size_t i = 0; i < c.size(); ++i)
        {
            std::cout << c[i] << " ";
        }
        std::cout << std::endl;
    }
    static void setNumberField(const NumberField& nf)
    {
        numberField_ = &nf;
    }
    static void clearNumberField()
    {
        numberField_ = 0;
    }
    static int degree()
    {
        if (!numberField_)
        {
            throw std::string("AlgebraicNumber::degree() : numberField_ not set");
        }
        return numberField_->degree();
    }
    static const VeryLong index()
    {
        if (!numberField_)
        {
            throw std::string("AlgebraicNumber::index() : numberField_ not set");
        }
        return numberField_->index();
    }
    static const VeryLong c_d()
    {
        if (!numberField_)
        {
            throw std::string("AlgebraicNumber::c_d() : numberField_ not set");
        }
        return numberField_->c_d();
    }

    static void createSpecialBasis(const VeryLong& p, std::vector<AlgebraicNumber>& j);

    Quotient<VeryLong> norm() const;

    Quotient<VeryLong> trace() const;

    static AlgebraicNumber read_algebraic_number(const char* an_str)
    {
        return Polynomial<Quotient<VeryLong> >::read_polynomial(an_str);
    }

    friend void read(const std::string& str, AlgebraicNumber& a)
    {
        a = Polynomial<Quotient<VeryLong> >::read_polynomial(str.c_str());
    }

    friend std::istream& operator>>(std::istream& is, AlgebraicNumber& a)
    {
        std::string str;
        std::getline(is, str);
        read(str, a);
        return is;
    }

    friend ostream& operator<<(ostream& os, const AlgebraicNumber& a)
    {
        int i = 0;
        const long int zero = 0L;
        for (auto& v: a.c_)
        {
            Quotient<VeryLong> value = v;
            int sign = 1;
            if (value < Quotient<VeryLong>(zero))
            {
                sign = -1;
                value = Quotient<VeryLong>(-1L) * value;
            }
            if (value != Quotient<VeryLong>(zero))
            {
                if (i > 0 && sign > 0) os << " + ";
                if (i > 0 && sign < 0) os << " - ";
                if (i == 0 && sign < 0) os << "-";
                os << value;
                if (i == 1) os << " alpha";
                else if (i > 1) os << " alpha^" << i;
            }
            i++;
        }

        return os;
    }

    static AlgebraicNumber& alpha();
    static const std::vector<AlgebraicNumber>& integralBasis();
    static const NumberField& nf()
    {
        if (!numberField_)
        {
            throw std::string("AlgebraicNumber::nf() : numberField_ not set");
        }
        return *numberField_;
    }

    const Quotient<VeryLong>& coefficient(int i) const
    {
        if (i < 0 || i >= static_cast<int>(c_.size()))
        {
            throw std::string("AlgebraicNumber::coefficient(): index out of range");
        }
        return c_[i];
    }
    const std::vector<Quotient<VeryLong> >& coefficients() const
    {
        return c_;
    }
    const std::vector<Quotient<VeryLong> >& ib_coefficients() const
    {
        make_ibc();
        return ibc_;
    }
    void set_coefficient(int j, const Quotient<VeryLong>& c)
    {
        if (j < 0 || j >= static_cast<int>(c_.size()))
        {
            throw std::string("AlgebraicNumber::set_coefficient(): index out of range");
        }
        c_[j] = c;
    }

    void ln_sigma(int j, long double& ln_re, long int& re_sign,
                  long double& ln_im, long int& im_sign);
    long double ln_sigma(int j);
    long double mod_sigma_2(int j) const;

    Polynomial<VeryLong> minimalPolynomial() const;

    AlgebraicNumber sqrt() const;

    VeryLongModular Phi(const VeryLong& N, const VeryLong& m) const;

private:
    void defineMatrix() const;
    void make_ibc() const
    {
        if (ibc_defined_) return;
        if (!numberField_)
        {
            throw std::string("AlgebraicNumber::make_ibc() : numberField_ not set");
        }
        ibc_ = numberField_->winv() * c_;
        ibc_defined_ = true;
    }
    static const NumberField* numberField_;
    // Coefficients of number in terms of alpha
    std::vector<Quotient<VeryLong > > c_;
    // Coefficients of number in terms of omega, the integral basis
    mutable std::vector<Quotient<VeryLong > > ibc_;
    mutable bool ibc_defined_;
    mutable Matrix<Quotient<VeryLong > > M_;
    mutable bool matrix_defined_;
};
#endif
