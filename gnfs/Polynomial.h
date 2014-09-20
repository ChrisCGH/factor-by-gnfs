#ifndef _POLYNOMIAL_H
#define _POLYNOMIAL_H
#include <vector>
#include <math.h>
#include <complex>
#include <algorithm>
#include <cstring>
using std::complex;
using std::norm;
#include "pow.h"
using std::vector;
#include "Quotient.h"
#include "crt.h"
#include "gcd.h"
#include "mt19937int.h"

template <class F> class Polynomial
{
   public:
      Polynomial() : degree(-1)
      {
         _coefficients.reserve(100L);
      }
      /*
      explicit Polynomial(long int li)
      {
         _coefficients.reserve(100L);
         _coefficients.resize(1);
         _coefficients[0] = F(li);
         check_degree();
      }
      */
      Polynomial(const F& f)
      {
         _coefficients.reserve(100L);
         _coefficients.resize(1);
         _coefficients[0] = f;
         check_degree();
      }
      Polynomial(int count, F* coefficient)
      {
         _coefficients.reserve(100L);
         _coefficients.resize(count);
         int i = 0;
         for (i = 0; i < count; ++i)
         {
            _coefficients[i] = coefficient[i];
         }
         check_degree();
      }
      Polynomial(const std::vector<F>& coefficient)
      {
         _coefficients.reserve(100L);
         _coefficients.resize(coefficient.size());
         int i = 0;
         for (typename std::vector<F>::const_iterator iter = coefficient.begin();
               iter != coefficient.end();
               ++iter)
         {
            _coefficients[i] = *iter;
            i++;
         }
         check_degree();
      }
   private:
      static int get_max_power(const std::vector<std::pair<int, std::string> >& terms)
      {
         int max_power = 0;
         for (size_t i = 0; i < terms.size(); ++i)
         {
            if (terms[i].first > max_power)
            {
               max_power = terms[i].first;
            }
         } 
         return max_power;
      }
   public:
      Polynomial(const std::vector<std::pair<int, std::string> >& terms)
      {
         _coefficients.resize(get_max_power(terms) + 1); 
         for (size_t i = 0; i < terms.size(); ++i)
         {
             _coefficients[terms[i].first] += F(terms[i].second);
         }
         check_degree();
      }
      ~Polynomial()
      {}
      Polynomial(const Polynomial& p)
      {
         _coefficients.reserve(100L);
         _coefficients = p._coefficients;
         check_degree();
      }

      Polynomial& operator= (const Polynomial& p)
      {
         if (this == &p) return *this;
         _coefficients.clear();
         _coefficients = p._coefficients;
         check_degree();
         return *this;
      }

      Polynomial operator- () const
      {
         Polynomial tmp(*this);
         for (int i = 0; i < deg() + 1; i++)
         {
            tmp._coefficients[i] = - _coefficients[i];
         }
         return tmp;
      }

      int operator< (const F& scalar) const
      {
         return 0;
      }
      int operator< (const Polynomial& p) const
      {
         if (deg() == p.deg())
         {
            for (int i = deg(); i >= 0; i--)
            {
               if (_coefficients[i] < p._coefficients[i]) return 1;
               if (_coefficients[i] > p._coefficients[i]) return 0;
            }
            return 0;
         }
         return (deg() < p.deg());
      }

/*
      int operator> (const Polynomial& p) const
      {
         return (deg() > p.deg());
      }
*/

      int operator== (const Polynomial& p) const
      {
         if (is_zero() && p.is_zero()) return 1;
         if (deg() != p.deg()) return 0;
         if (deg() == -1 && p.deg() == -1) return 1;
         typename std::vector<F>::const_iterator iter_p = p._coefficients.begin();
         for (typename std::vector<F>::const_iterator iter = _coefficients.begin();
               iter != _coefficients.end();
               ++iter)
         {
            if (*iter != *iter_p) return 0;
            ++iter_p;
         }
         return 1;
      }
      int operator== (const F& scalar) const
      {
         if (deg() == 0 && scalar == _coefficients[0]) return 1;
         if (deg() == -1 && scalar == F(0L)) return 1;
         return 0;
      }
      int is_zero() const
      {
         if (deg() == 0 && F(0L) == _coefficients[0]) return 1;
         if (deg() == -1) return 1;
         return 0;
      }
      int operator!= (const Polynomial& p)
      {
         return !(*this == p);
      }

      Polynomial& operator += (const Polynomial& p)
      {
         int at_end = 0;
         typename std::vector<F>::iterator iter = _coefficients.begin();
         for (typename std::vector<F>::const_iterator iter_p = p._coefficients.begin();
               iter_p != p._coefficients.end();
               ++iter_p)
         {
            if (at_end || (iter == _coefficients.end()))
            {
               at_end = 1;
               _coefficients.push_back(*iter_p);
            }
            else
            {
               (*iter) += (*iter_p);
               ++iter;
            }
         }
         check_degree();
         return *this;
      }

      Polynomial& operator -= (const Polynomial& p)
      {
         int at_end = 0;
         typename std::vector<F>::iterator iter = _coefficients.begin();
         for (typename std::vector<F>::const_iterator iter_p = p._coefficients.begin();
               iter_p != p._coefficients.end();
               ++iter_p)
         {
            if (at_end || (iter == _coefficients.end()))
            {
               at_end = 1;
               F value = 0L - (*iter_p);
               _coefficients.push_back(value);
            }
            else
            {
               (*iter) -= (*iter_p);
               ++iter;
            }
         }
         check_degree();
         return *this;
      }

      Polynomial& operator*= (const Polynomial& p)
      {
         Polynomial pp;
         if (deg() < 0 || p.deg() < 0)
         {
            *this = pp;
            return *this;
         }
         int degree = deg() + p.deg();
         pp._coefficients.resize(degree + 1, F(0L));
         int i = 0;
         int j = 0;

         for (typename std::vector<F>::const_iterator iter_p1 = _coefficients.begin();
               iter_p1 != _coefficients.end();
               ++iter_p1)
         {
            j = 0;
            for (typename std::vector<F>::const_iterator iter_p2 = p._coefficients.begin();
                  iter_p2 != p._coefficients.end();
                  ++iter_p2)
            {
               F res = *iter_p1 * *iter_p2;
               pp._coefficients[i + j] = pp._coefficients[i + j] + res;
               j++;
            }
            i++;
         }
         *this = pp;
#if 0
         Polynomial p1 = p * (*this);
         *this = p1;
#endif
         check_degree();
         return *this;
      }

      friend Polynomial operator+(const Polynomial& p1, const Polynomial& p2)
      {
         Polynomial p = p1;
         int at_end = 0;

         typename std::vector<F>::iterator iter = p._coefficients.begin();
         for (typename std::vector<F>::const_iterator iter_p2 = p2._coefficients.begin();
               iter_p2 != p2._coefficients.end();
               ++iter_p2)
         {
            if (at_end || (iter == p._coefficients.end()))
            {
               at_end = 1;
               p._coefficients.push_back(*iter_p2);
            }
            else
            {
               (*iter) = (*iter) + (*iter_p2);
               ++iter;
            }

         }
         p.check_degree();
         return p;
      }

      friend Polynomial operator-(const Polynomial& p1, const Polynomial& p2)
      {
         Polynomial p = p1;
         int at_end = 0;

         typename std::vector<F>::iterator iter = p._coefficients.begin();
         for (typename std::vector<F>::const_iterator iter_p2 = p2._coefficients.begin();
               iter_p2 != p2._coefficients.end();
               ++iter_p2)
         {
            if (at_end || (iter == p._coefficients.end()))
            {
               at_end = 1;
               p._coefficients.push_back(F(0L)-(*iter_p2));
            }
            else
            {
               *iter = *iter - *iter_p2;
               ++iter;
            }

         }
         p.check_degree();
         return p;
      }

      Polynomial& operator-= (const F& value)
      {
         if (deg() == -1)
         {
            _coefficients.push_back(-value);
            degree = 0;
         }
         else
         {
            _coefficients[0] -= value;
            if (deg() == 0 && _coefficients[0] == F(0L)) degree = -1;
         }
         return *this;
      }
      Polynomial& operator+= (const F& value)
      {
         if (deg() == -1)
         {
            _coefficients.push_back(value);
            degree = 0;
         }
         else
         {
            _coefficients[0] += value;
            if (deg() == 0 && _coefficients[0] == F(0L)) degree = -1;
         }
         return *this;
      }

      friend Polynomial operator*(const Polynomial& p1, const F& scalar)
      {
         return p1 * Polynomial<F> (scalar);
      }

      friend Polynomial operator*(const F& scalar, const Polynomial& p1)
      {
         return p1 * Polynomial<F> (scalar);
      }

      friend Polynomial<F> operator*(const Polynomial<F>& p1, const Polynomial<F>& p2)
      {
         Polynomial p;

         if (p1.deg() < 0 || p2.deg() < 0) return p;
         int degree = p1.deg() + p2.deg();
         p._coefficients.resize(degree + 1, F(0L));
         int i = 0;
         int j = 0;

         for (typename std::vector<F>::const_iterator iter_p1 = p1._coefficients.begin();
               iter_p1 != p1._coefficients.end();
               ++iter_p1)
         {
            j = 0;
            for (typename std::vector<F>::const_iterator iter_p2 = p2._coefficients.begin();
                  iter_p2 != p2._coefficients.end();
                  ++iter_p2)
            {
               F res = *iter_p1 * *iter_p2;
               p._coefficients[i + j] = p._coefficients[i + j] + res;
               j++;
            }
            i++;
         }
         p.check_degree();
         return p;
      }

      friend Polynomial operator/(const Polynomial& p1, const F& value)
      {
         Polynomial p = p1;
         if (p.deg() < 0) return p;

         for (typename std::vector<F>::iterator iter = p._coefficients.begin();
               iter != p._coefficients.end();
               ++iter)
         {
            *iter /= value;
         }
         p.check_degree();
         return p;
      }

      Polynomial& operator/=(const Polynomial& p2)
      {
         *this = *this / p2;
         return *this;
      }

      friend Polynomial operator/(const Polynomial& p1, const Polynomial& p2)
      {
         int d1 = p1.deg();
         int d2 = p2.deg();
         if (d1 < d2) return Polynomial();
         F b = p2._coefficients[d2];
         F c;
         Polynomial<F> result;
         int index = d1 - d2;
         result._coefficients.resize(d1 - d2 + 1, F(0L));
         result.degree = d1 - d2;
         Polynomial tmp1 = p1;

         while (index >= 0)
         {
            d1 = tmp1.deg();
            c = tmp1._coefficients[d1] / b;
            result._coefficients[index] = c;
            typename std::vector<F>::const_reverse_iterator iter2 = p2._coefficients.rbegin();
            for (typename std::vector<F>::reverse_iterator iter = tmp1._coefficients.rbegin();
                  iter != tmp1._coefficients.rend() &&
                  iter2 != p2._coefficients.rend();
                  ++iter,++iter2)
            {
               *iter = *iter - ((*iter2) * c);
            }
            tmp1._coefficients[d1] = F(0L);
            tmp1.check_degree();
            index = tmp1.deg() - d2;
         }

         return result;
      }

      // Algorithm 3.1.2 (Pseudo-Division)
      friend void pseudo_divide(const Polynomial& A,
                                const Polynomial& B,
                                Polynomial& Q,
                                Polynomial& R)
      {
         int m = A.deg();
         int n = B.deg();
         F d = B._coefficients[n];
         // Assume m >= n
         // Step 1. [Initialize]
         F zero(0L);
         R = A;
         Q = Polynomial<F>(zero);
         int e = m - n + 1;

         while (R.deg() >= B.deg())
         {
            // Step 3. [Find coefficient]
            int c = R.deg() - B.deg();
            std::vector<F> s;
            s.resize(c + 1);
            s[c] = R._coefficients[R.deg()];
            Polynomial<F> S(s);
            Q = d * Q + S;
            R = d * R - S * B;
            e--;
         }

         // Step 2. [Finished]
         F q = pow<F, int>(d, e);
         Q = q * Q;
         R = q * R;
      }

/*
      friend Polynomial pseudo_divide(const Polynomial& p1, const Polynomial& p2)
      {
         int d1 = p1.deg();
         int d2 = p2.deg();
         if (d1 < d2) return Polynomial();
         F b = p2._coefficients[d2];
         F c;
         Polynomial<F> result;
         int index = d1 - d2;
         result._coefficients.resize(d1 - d2 + 1, 0L);
         result.degree = d1 - d2;
         Polynomial tmp1 = p1;

         while (index >= 0)
         {
            d1 = tmp1.deg();
            F q = gcd(tmp1._coefficients[d1], b);
            if (q != b)
            {
               tmp1 = tmp1 * b;
               for (int i = 0; i < index; i++) result._coefficients[i] *= b;
            }
            c = tmp1._coefficients[d1] / b;
            result._coefficients[index] = c;
            typename std::vector<F>::const_reverse_iterator iter2 = p2._coefficients.rbegin();
            for (typename std::vector<F>::reverse_iterator iter = tmp1._coefficients.rbegin();
                  iter != tmp1._coefficients.rend() &&
                  iter2 != p2._coefficients.rend();
                  ++iter,++iter2)
            {
               *iter = *iter - ((*iter2) * c);
            }
            tmp1.check_degree();
            index = tmp1.deg() - d2;
         }

         return result;
      }
*/

      Polynomial& divide_by_X()
      {
         if (deg() < 0) 
         {
            return *this;
         }
         if (deg() == 0) 
         {
            --degree;
            _coefficients.clear();
            return *this; 
         }

         int i = 0;

         typename std::vector<F>::const_iterator iter = _coefficients.begin();
         ++iter; // skip lowest coefficient
         for (;
               iter != _coefficients.end();
               ++iter)
         {
            _coefficients[i] = *iter;
            i++;
         }
         degree--;
         _coefficients.resize(degree+1);
         return *this;
      }

      friend Polynomial operator%(const Polynomial& p1, const Polynomial& p2)
      {
         int d1 = p1.deg();
         int d2 = p2.deg();
         if (d1 < d2) return p1;
         F b = p2._coefficients[d2];
         F c;
         int index = d1 - d2;
         Polynomial tmp1 = p1;

         while (index >= 0)
         {
            d1 = tmp1.deg();
            c = tmp1._coefficients[d1] / b;
            if (c * b != tmp1._coefficients[d1])
            {
               throw std::string("Polynomial operator%(p1, p2) : coefficient of highest power in p2 must divide coefficient of highest power in p1");
            }
            typename std::vector<F>::const_reverse_iterator iter2 = p2._coefficients.rbegin();
            for (typename std::vector<F>::reverse_iterator iter = tmp1._coefficients.rbegin();
                  iter != tmp1._coefficients.rend() &&
                  iter2 != p2._coefficients.rend();
                  ++iter,++iter2)
            {
               if (*iter2 != F(0L))
               {
                  *iter -= ((*iter2) * c);
               }
            }
            tmp1.check_degree();
            index = tmp1.deg() - d2;
         }

         return tmp1;
      }

      friend Polynomial remainder(const Polynomial& p1, const Polynomial& p2)
      {
         long int d1 = p1.deg();
         long int d2 = p2.deg();
         if (d1 < d2) return p1;
         F multiplier = pow<F, long int>(p2.coefficient(d2), (long int)(d1 - d2 + 1));
         return (p1 * multiplier) % p2;
      }

      Polynomial& make_primitive()
      {
         F cont = content();
         if (cont != F(1L))
         {
            *this = *this / cont;
         }
         return *this;
      }

      Polynomial& make_monic()
      {
         if (deg() < 0) return *this;
         F c = _coefficients[deg()];
         if (c != F(1L))
         {
            *this = *this / c;
         }
         return *this;
      }

      Polynomial derivative() const
      {
         if (deg() == 0) return Polynomial();

         Polynomial result;
         result._coefficients.resize(deg());
         typename std::vector<F>::const_iterator iter = _coefficients.begin();
         ++iter;
         long int i = 1;
         F tmp;
         for (;
               iter != _coefficients.end();
               ++iter)
         {
            tmp = *iter * F(i);
            result._coefficients[i - 1] = tmp;
            i++;
         }
         result.check_degree();
         return result;
      }

      F discriminant() const
      {
         // ??? NOT YET IMPLEMENTED ???
         return F(0L);
      }

      int deg() const
      {
         return degree;
      }

      F content() const
      {
         // content is the GCD of the coefficients
         typename std::vector<F>::const_iterator iter = _coefficients.begin();
         if (iter == _coefficients.end())
         {
             return F(0L);   
         }
         F c(*iter);
         for (;
               iter != _coefficients.end();
               ++iter)
         {
            if (*iter != F(0L)) c = gcd<F>(c, *iter);
         }
         return c;
      }

      F coefficient(int i) const
      {
         if (i < 0) 
         {
            throw std::string("Polynomial::coefficient() : index out of range");
         }
         if (i > deg())
         {
             return F(0L);
	 }
         return _coefficients[i];
      }

      friend std::ostream& operator<< (std::ostream& os, const Polynomial<F>& p)
      {
         int i = 0;
         const long zero = 0;
         for (typename std::vector<F>::const_iterator iter = p._coefficients.begin();
               iter != p._coefficients.end();
               ++iter)
         {
            F value = *iter;
            int sign = 1;
            if (value < F((const long)zero))
            {
               sign = -1;
               value = -value;
            }
            if (value != F((const long)zero))
            {
               if (i > 0 && sign > 0) os << " + ";
               if (i > 0 && sign < 0) os << " - ";
               if (i == 0 && sign < 0) os << "-";
               os << value;
               if (i == 1) os << " X";
               else if (i > 1) os << " X^" << i;
            }
            i++;
         }

         return os;
      }

      F evaluate(const F& value) const
      {
         // evaluate the polynomial
         F result(0L);

         for (typename std::vector<F>::const_reverse_iterator iter = _coefficients.rbegin();
               iter != _coefficients.rend();
               ++iter)
         {
            result *= value;
            result += *iter;
         }
         return result;
      }

      F evaluate_homogeneous(const F& a, const F& b) const
      {
         // evaluate the homogeneous polynomial F(X,Y)
         // corresponding to f(X) at (a,b)
         F result = 0L;
         F temp = 1L;
         for (typename std::vector<F>::const_iterator iter = _coefficients.begin();
               iter != _coefficients.end();
               ++iter)
         {
            result *= b;
            result += temp * (*iter);
            temp *= a;
         }
         return result;
      }

      F evaluate_homogeneous_1(const F& a, const F& b) const
      {
         F result = 0L;
         int i = 0;
         F y = 1L;
         for (i = 0; i < deg(); i++) y *= b;
         F x = 1L;
         result = coefficient(0) * y * x;
         for (i = 1; i <= deg(); i++)
         {
            y /= b;
            x *= a;
            result += coefficient(i) * y * x;
         }
         return result;
      }

      Polynomial evaluate(const Polynomial& value) const
      {
         // evaluate the polynomial
         Polynomial result = F(0L);

         for (typename std::vector<F>::const_reverse_iterator iter = _coefficients.rbegin();
               iter != _coefficients.rend();
               ++iter)
         {
            result *= value;
            result += *iter;
         }
         return result;
      }

      void set_coefficient(long int i, const F& value)
      {
         if (i < 0 || i > deg()) return;
         _coefficients[i] = value;
      }

      F max_coefficient() const
      {
         F max_coeff = _coefficients[0];
         if (max_coeff < F(0L)) max_coeff = -max_coeff;
         for (int i = 1; i <= deg(); i++)
         {
            F coeff = _coefficients[i];
            if (coeff < F(0L)) coeff = -coeff;
            if (coeff > max_coeff) max_coeff = coeff;
         }
         return max_coeff;
      }
      static Polynomial<F> read_polynomial(const char* poly_str)
      {
         // assumes format
         // 7111977108918472837611244937303928 - 2523764323493368834122550403406 X - 100831972740733548080518208
         // X^2 + 7467483273731340429431 X^3 + 150036498813232980 X^4 + 2859632424000 X^5
         // BNF definition:
         // polynomial := term |
         //               term + whitespace + polynomial
         // term       := [ sign ] + [ whitespace ] + [ number + [ whitespace ] + [ "*" ] ] + whitespace + letter + [ "^" + number ]
         // letter     := a-zA-Z
         // number     := [ { digit } ] + [ "." ] + { digit } + [ "e" + { digit } ]
         // digit      := 0-9
         // whitespace := { space | tab }
         //
      
         try
         {
             std::string s(poly_str);
             std::vector<std::pair<int, std::string> > terms;
          
             std::string new_s;
             std::pair<int, std::string> term;
             bool first_term = true;
             char variable = '\0';
             char term_variable = '\0';
             while (get_next_term(s, new_s, term, first_term, term_variable))
             {
                first_term = false;
                if (variable != term_variable)
                {
                   if (variable == '\0')
                   {
                      variable = term_variable;
                   }
                   else
                   {
                      throw std::string("Only mono-variate polynomials supported");
                   }
                }
                
                s = new_s;
                terms.push_back(term);
             }
            
             return Polynomial<F>(terms); 
         }
         catch (const std::string& err)
         {
            std::ostringstream oss;
            oss << "Polynomial::read_polynomial() : Failed to parse input : " << err;
            throw oss.str();
         }
         
         return Polynomial<F>();
      }
      static Polynomial<F> read_polynomial(std::istream& is);
      template <typename DOUBLE>
      static Polynomial<DOUBLE> convert_to_double(const Polynomial<F>& poly)
      {
        std::vector<DOUBLE > dcoefficients;
        dcoefficients.resize(poly.deg() + 1);

        int k = 0;
        for (k = 0; k < poly.deg() + 1; k++)
        {
            F c = poly.coefficient(k);
            dcoefficients[k] = c;
        }
        Polynomial<DOUBLE> dpoly(dcoefficients);
        return dpoly;
      }
      static void factor(const Polynomial<F>& AA, std::vector<Polynomial<F> >& factors, F& cont);
   private:
      static bool parse_number(const char*& c)
      {
          // e.g. 1234, 1234.0, 1234.0e10, 1234.0e+10, 1234.0e-10, .023
          // 
          char* tail;
          const char* d = c;
          int slash_count = 0;
          const char* slash = 0;
          while (c && *c != '\0' && (isdigit(*c) || *c == '.' || *c == 'e' || *c == 'E' || *c == '+' || *c == '-' || *c == '/'))
          {
              if (*c == '/')
              {
                 ++slash_count;
                 slash = c;
              }
              ++c;
          }

          if (slash_count == 0)
          {
              // what if d is something like " Y" ???
              if (!isdigit(*d))
              {
                 return false;
              }
              ::strtod(d, &tail);
              if (tail != c)
              {
                 if (*c == '^' && tail == c - 1)
                 {
                    c -= 1;
                    return true;
                 }
                 return false;
              }
          }
          else
          {
              if (slash_count > 1)
              {
                 return false;
              }

              // slash_count is 1, so we may have a quotient of 2 numbers
              std::string num1(d, slash - d);
              const char* c1 = num1.c_str();
              std::string num2(slash + 1);
              const char* c2 = num2.c_str();
              if (parse_number(c1) && parse_number(c2))
              {
                 //std::cout << "num1 = [" << num1 << "], num2 = [" << num2 << "]" << std::endl;
                 return true;
              }
              return false;
          }
          return true;
      }
      static bool get_next_term(const std::string& s, std::string& new_s, std::pair<int, std::string>& term, bool first_term, char& variable)
      {
          //std::cout << "get_next_term, s = [" << s << "]" << std::endl;
          if (s.empty())
          {
              return false;
          }
      
          const char* c = s.c_str();
          // skip any leading whitespace
          while (c && *c != '\0' && (*c == ' ' || *c == '\t'))
          {
              ++c;
          }
      
          if (!c || *c == '\0')
          {
             return false; 
          }
      
          // check for sign
          if (!first_term)
          {
              if (*c != '-' && *c != '+')
              {
                 std::ostringstream err;
                 err << "Bad polynomial term : [" << s << "]";
                 throw err.str();
              }
          }
      
          int sign = 1;
          if (*c == '-')
          {
              sign = -1;
          }
         
          // skip past sign and any more white space 
          while (c && *c != '\0' && (*c == '+' || *c == '-' || *c == ' ' || *c == '\t'))
          {
              ++c;
          }
          
          if (!c || *c == '\0' || (!isdigit(*c) && !isalpha(*c) && *c != '.'))
          {
             std::ostringstream err;
             err << "Bad polynomial term : [" << s << "]";
             throw err.str();
          }
      
          term.first = 0;
          term.second = "0";
          new_s = "";
      
          if (isdigit(*c) || *c == '.')
          {
              const char* d = c;
              if (!parse_number(c))
              {
                  std::ostringstream err;
                  err << "Bad polynomial term : [" << s << "]";
                  throw err.str();
              }

              term.second = std::string(d, c - d); 
              if (sign < 0)
              {
                  term.second = "-" + term.second;
              }
      
              // skip any more whitespace
              while (c && *c != '\0' && (*c == ' ' || *c == '\t'))
              {
                 ++c;
              }
      
              if (!c || *c == '\0')
              {
                 //std::cout << "1. Term : " << term.first << " : " << term.second << ", new_s = [" << new_s << "]" << std::endl;
                 return true; 
              }
      
              if (!isalpha(*c))
              {
                  new_s = c;
                  //std::cout << "2. Term : " << term.first << " : " << term.second << ", new_s = [" << new_s << "]" << std::endl;
                  return true;
              }
          }
      
          variable = *c;
          if (term.second == "0")
          {
              // we hit 'X' without first hitting a number, set to +/- 1 depending on sign
              term.second = "1";
              if (sign < 0)
              {
                  term.second = "-1";
              }
          }

          if (std::strlen(c) >= 5 && std::string(c, 5) == std::string("alpha"))
          {
              c += 5;
          }  
          else
          {
              ++c;
          }
          
          if (!c || *c == '\0')
          {
              term.first = 1;
              //std::cout << "3. Term : " << term.first << " : " << term.second << ", new_s = [" << new_s << "]" << std::endl;
              return true;
          }
      
          if (*c != '^')
          {
              term.first = 1;
              new_s = c;
              //std::cout << "4. Term : " << term.first << " : " << term.second << ", new_s = [" << new_s << "]" << std::endl;
              return true;
          }
      
          ++c;
      
          if (!c || *c == '\0' || !isdigit(*c))
          {
             std::ostringstream err;
             err << "Bad polynomial term : [" << s << "]";
             throw err.str();
          }
      
         
          const char* d = c;
          while (c && isdigit(*c))
          {
              ++c;
          } 
      
          term.first = ::atoi(std::string(d, c - d).c_str());
          new_s = c;
          //std::cout << "5. Term : " << term.first << " : " << term.second << ", new_s = [" << new_s << "]" << std::endl;
          return true;
      }
      vector<F> _coefficients;
      int degree;
      void check_degree()
      {
         // check that degree is correct and adjust accordingly
         degree = static_cast<int>(_coefficients.size() - 1);
         int d = degree;
         while (degree >= 0 && _coefficients[degree] == F(0L))
         {
            degree--;
         }
         if (d > degree) _coefficients.resize(degree + 1);
      }
};

// Algorithm 3.6.6 (Complex Roots)
template <class T > int find_roots_over_C(const Polynomial<complex<T > >& P, std::vector<complex<T > >& roots)
{
   // Quote from H. Cohen "A Course in Computational Algebraic Number Theory", p. 146 :
   // "Given a squarefree polynomial P, this algorithm outputs its complex roots (in a random order).
   // In quite rare cases the algorithm may fail. On the other hand it is absolutely necessary that
   // the polynomial be squarefree (this can be achieved by replacing P by P/(P,P'))."

   // Step 1. [Initializations]
   const T ACCURACY = T(1L) / (T(100000000L) * T(100000000L) * T(100000000L));
   Polynomial<complex<T > > Q = P;
   Polynomial<complex<T > > P1 = P.derivative();
   Polynomial<complex<T > > Q1 = P1;
   //std::cout << "find_roots_over_C of poly " << P << std::endl;
   //std::cout << "P'(X) = " << P1 << std::endl;
   int n = P.deg();
   // f = 1 if P has real coefficients, f = 0 otherwise
   int f = 1;
   for (int i = 0; i < n; i++)
   {
      complex<T > c = P.coefficient(i);
      if (!c.real()) f = 0;
   }
   //if (f) std::cout << "P has real coefficients" << std::endl;
   //else std::cout << "P has complex coefficients" << std::endl;

   while (n > 1)
   {
      // Step 2. [Initialize root finding]
      complex<T> x = complex<T>(T(13L) / T(10L), T(314159L)/T(1000000L));
      complex<T> v = Q.evaluate(x);
      T m = norm(v);
      //std::cout << "Step 2., x = " << x << ", v = Q(x) = " << v << ", m = |v|^2 = " << v << std::endl;

      int done = 0;
      while (!done)
      {
         // Step 3. [Initialize recursion]
         int c = 0;
         complex<T> dx = v / Q1.evaluate(x);
         //std::cout << "Step 3., dx = v / Q'(x) = " << dx << std::endl;
         if (norm(dx) > ACCURACY)
         {
            int done1 = 0;
            while (!done1)
            {
               // Step 4. [Try a lambda]
               complex<T> y = x - dx;
               complex<T> v1 = Q.evaluate(y);
               T m1 = norm(v1);
               //std::cout << "Step 4., y = " << y << ", v1 = " << v1 << ", m1 = " << m1 << ", m = " << m << std::endl;
               if (m1 <= m)
               {
                  x = y;
                  v = v1;
                  m = m1;
                  done1 = 1;
               }
               else
               {
                  c++;
                  dx /= T(4L);
                  if (c >= 20)
                  {
                     // failed
                     std::cerr << "failed to find roots" << std::endl;
                     return 0;
                  }
               }
            }
         }
         else done = 1;
      }

      // Step 5. [Polish root]
      //std::cout << "Step 5., x = " << x << std::endl;
      x = x - P.evaluate(x) / P1.evaluate(x);
      x = x - P.evaluate(x) / P1.evaluate(x);
      //std::cout << "after Step 5., x = " << x << std::endl;

      // Step 6. [Divide]
      //std::cout << "Step 6." << std::endl;
      if (f == 0 || (f == 1 &&
                     ((x.imag() > T(0L) && x.imag() < ACCURACY) ||
                      (x.imag() < T(0L) && T(-1L)*x.imag() < ACCURACY))))
      {
         //std::cout << "one real root, x = " << x << std::endl;
         x = complex<T>(x.real(), T(0L));
         roots.push_back(x);
         std::vector<complex<T > > co;
         co.resize(2);
         co[1] = complex<T>(1L);
         co[0] = complex<T>(0L) - x;
         Polynomial<complex<T > > X_minus_x(co);
         Q = Q / X_minus_x;
         n--;
      }
      else
      {
         //std::cout << "two conjugate roots, x = " << x << std::endl;
         roots.push_back(x);
         roots.push_back(conj(x));
         std::vector<complex<T > > co;
         co.resize(3);
         co[2] = complex<T>(T(1L));
         co[1] = complex<T>(T(1L)) - T(2L) * x.real();
         co[0] = complex<T>(norm(x));
         Polynomial<complex<T > > poly(co);
         //std::cout << "poly = " << poly << std::endl;
         Q = Q / poly;
         //std::cout << "new Q = " << Q << std::endl;
         //std::cout << "check: Q * poly = " << Q * poly << std::endl;
         n -= 2;
      }
   }

   if (n == 1)
   {
      // just a linear equation to solve
      // Q(X) = z X + w = 0 => X = -w / z
      complex<T> x = complex<T>(0L) - Q.coefficient(0) / Q.coefficient(1);
      roots.push_back(x);
   }


   return P.deg();
}

template <class FLOAT>
FLOAT minimize_size_over_Re(const Polynomial<complex<FLOAT > >& P, complex<FLOAT>& x)
{
   FLOAT s = x.real();
   FLOAT last_s = 0.0;
   const long int SAMPLE_SIZE = 10;
   const FLOAT MAX_SAMPLE_RANGE = FLOAT(1e15);
   FLOAT s_min = -MAX_SAMPLE_RANGE;
   FLOAT s_max = MAX_SAMPLE_RANGE;
   FLOAT s_sample_size = SAMPLE_SIZE;
   FLOAT s_delta = (s_max - s_min) / s_sample_size;
   FLOAT s_try = s_min;
   FLOAT s_diff = 0.0;
   FLOAT value;
   FLOAT min_value = norm(P.evaluate(complex<FLOAT>(s, x.imag())));
   // return s;
   int done = 0;
   int iterations = 0;

   while (!done)
   {
      s_try = s_min;
      for (int l = 0; l < s_sample_size; l++)
      {
         value = norm(P.evaluate(complex<FLOAT>(s_try, x.imag())));
         if (value < min_value)
         {
            min_value = value;
            last_s = s;
            s = s_try;
         }
         s_try += s_delta;
      }
      s_min = s - s_sample_size*s_delta / FLOAT(10.0);
      s_max = s + s_sample_size*s_delta / FLOAT(10.0);
      s_delta = (s_max - s_min) / s_sample_size;
      s_diff = fabs(last_s - s);

      last_s = s;
      const FLOAT epsilon(1e-60);
      iterations++;
      if (min_value < FLOAT(1e-10) || (s_diff < epsilon && s_delta < epsilon) || iterations > 1000) done = 1;
   }
   value = norm(P.evaluate(complex<FLOAT>(s, x.imag())));
   x = complex<FLOAT>(s, x.imag());
   //cout << "minimize_size_over_Re: x = " << x << ", value = " << value << endl;
   return value;
}

template <class FLOAT>
FLOAT minimize_size_over_Im(const Polynomial<complex<FLOAT > >& P, complex<FLOAT>& x)
{
   FLOAT s = x.imag();
   FLOAT last_s = 0.0;
   const long int SAMPLE_SIZE = 10;
   const FLOAT MAX_SAMPLE_RANGE = FLOAT(1e15);
   FLOAT s_max = MAX_SAMPLE_RANGE;
   FLOAT s_min = -MAX_SAMPLE_RANGE;
   FLOAT s_sample_size = SAMPLE_SIZE;
   FLOAT s_delta = (s_max - s_min) / s_sample_size;
   FLOAT s_try = s_min;
   FLOAT s_diff = 0.0;
   FLOAT value;
   FLOAT min_value = norm(P.evaluate(complex<FLOAT>(x.real(),s)));
   // return s;
   int done = 0;
   int iterations = 0;

   while (!done)
   {
      s_try = s_min;
      for (int l = 0; l < s_sample_size; l++)
      {
         value = norm(P.evaluate(complex<FLOAT>(x.real(),s_try)));
         if (value < min_value)
         {
            min_value = value;
            last_s = s;
            s = s_try;
         }
         s_try += s_delta;
      }
      s_min = s - s_sample_size*s_delta / FLOAT(10.0);
      s_max = s + s_sample_size*s_delta / FLOAT(10.0);
      s_delta = (s_max - s_min) / s_sample_size;
      s_diff = fabs(last_s - s);

      last_s = s;
      if (iterations > 1000)
      {
         //cout << s_diff << endl;
      }
      const FLOAT epsilon(1e-60);
      iterations++;
      if (min_value < FLOAT(1e-10) || (s_diff < epsilon && s_delta < epsilon) || iterations > 1000) done = 1;
   }
   value = norm(P.evaluate(complex<FLOAT>(x.real(),s)));
   x = complex<FLOAT>(x.real(), s);
   //cout << "minimize_size_over_Im: x = " << x << ", value = " << value << endl;
   return value;
}

template <class FLOAT>
FLOAT minimize_size(const Polynomial<complex<FLOAT > >& P, complex<FLOAT>& x)
{
   const long int SAMPLE_SIZE = 10;
   const long int sample_size = SAMPLE_SIZE;
   const FLOAT r(10.0); // length of vector;
   static FLOAT s_vec[sample_size];
   static FLOAT t_vec[sample_size];
   static int first_time = 1;

   if (first_time)
   {
      first_time = 0;
      FLOAT theta_inc = FLOAT(2.0) * FLOAT(M_PI) / FLOAT(sample_size);
      FLOAT theta = 0.0;
      for (int i = 0; i < sample_size; i++)
      {
         s_vec[i] = (FLOAT)cos (theta);
         t_vec[i] = (FLOAT)sin (theta);
         theta += theta_inc;
      }
   }
   FLOAT s = x.real();
   FLOAT t = x.imag();

   FLOAT min_value = norm(P.evaluate(x));
   FLOAT initial_min_value = min_value;
   FLOAT min_s = x.real();
   FLOAT min_t = x.imag();
   int min_i = -1;
   FLOAT s_try;
   FLOAT t_try;
   for (int i = 0; i < sample_size; i++)
   {
      s_try = s + r * s_vec[i];
      t_try = t + r * t_vec[i];
      FLOAT value = norm(P.evaluate(complex<FLOAT>(s_try, t_try)));
      if (value < min_value)
      {
         min_value = value;
         min_i = i;
      }
   }
   if (min_value < initial_min_value)
   {
      // use sampling to get minimum point in this direction
      int done = 0;
      FLOAT r_min(0.0);
      FLOAT r_max(1e15);//MAX_SAMPLE_RANGE;
      FLOAT r_sample_size(SAMPLE_SIZE);
      FLOAT r_delta = (r_max - r_min) / r_sample_size;
      FLOAT r_try = r_min;
      FLOAT r_diff(0.0);
      FLOAT rr(0.0);
      FLOAT value;
      FLOAT last_r;
      while (!done)
      {
         r_try = r_min;
         t_try = t + r_try * t_vec[min_i];
         s_try = s + r_try * s_vec[min_i];
         for (int l = 0; l < r_sample_size; l++)
         {
            value = norm(P.evaluate(complex<FLOAT>(s_try, t_try)));
            if (value < min_value)
            {
               min_value = value;
               last_r = rr;
               rr = r_try;
               min_s = s_try;
               min_t = t_try;
               //cout << "(s,t) = (" << s_try << ", " << t_try << "), value = " << value << endl;
            }
            r_try += r_delta;
            t_try += r_delta * t_vec[min_i];
            s_try += r_delta * s_vec[min_i];
         }
         r_min = rr - r_sample_size*r_delta / FLOAT(10.0);
         if (r_min < FLOAT(0.0)) r_min = FLOAT(0.0);
         r_max = rr + r_sample_size*r_delta / FLOAT(10.0);
         r_delta = (r_max - r_min) / r_sample_size;
         r_diff = last_r - rr;
         if (r_diff < FLOAT(0.0)) r_diff = -r_diff;
         last_r = rr;

         const FLOAT epsilon(1e-60);
         if (r_diff < epsilon && r_delta < epsilon) done = 1;
      }
      // min point found
      s = min_s;
      t = min_t;
   }

   x = complex<FLOAT>(s, t);
   //cout << "minimize_size: x = " << x << ", value = " << min_value << endl;
   return min_value;
}

template <class FLOAT>
int find_roots_over_C_q1(const Polynomial<complex<FLOAT > >& P, std::vector<complex<FLOAT > >& roots)
{
   Polynomial<complex<FLOAT > > Q = P;
   Polynomial<complex<FLOAT > > Q1 = P.derivative();
   const FLOAT ACCURACY(1e-10);
   const FLOAT ACCURACY1(1e-20);
   int roots_found = 0;
   while (roots_found < P.deg())
   {
      complex<FLOAT> x((FLOAT)1.3, (FLOAT)0.314159);
      FLOAT s(1.0);

      int i = 0;
      while (s > ACCURACY && i < 100)
      {
         s = minimize_size_over_Re(Q, x);
         s = minimize_size_over_Im(Q, x);
         s = minimize_size(Q, x);
         i++;
      }
      if (i < 100)
      {
         //cout << "Found root: x = " << x << ", value = " << s << endl;
         if (fabs(x.imag()) < ACCURACY)
         {
            //cout << "one real root, x = " << x << endl;
            x = complex<FLOAT>(x.real(), FLOAT(0.0));
            //Polish root
            int j = 0;
            while (s > ACCURACY1 && j < 100)
            {
               x = x - Q.evaluate(x) / Q1.evaluate(x);
               s = norm(Q.evaluate(x));
               j++;
            }
            roots.push_back(x);
            std::vector<complex<FLOAT> > a;
            int d = Q.deg();
            a.resize(Q.deg());
            a[d-1] = Q.coefficient(d);
            for (int i = d-1; i >= 1; --i)
            {
               a[i-1] = Q.coefficient(i) + x * a[i];
            }
            Polynomial<complex<FLOAT> > QQ(a);
            //cout << "QQ = " << QQ;
            //cout << endl;
            Q = QQ;
            Q1 = Q.derivative();
            roots_found++;
         }
         else
         {
            //cout << "two conjugate roots, x = " << x << endl;
            //Polish root
            int j = 0;
            while (s > ACCURACY1 && j < 100)
            {
               x = x - Q.evaluate(x) / Q1.evaluate(x);
               s = norm(Q.evaluate(x));
               j++;
            }
            roots.push_back(x);
            roots.push_back(conj(x));
            std::vector<complex<FLOAT > > co;
            co.resize(3);
            co[2] = complex<FLOAT>((FLOAT)1.0);
            co[1] = complex<FLOAT>((FLOAT)0.0) - (FLOAT)2.0 * x.real();
            co[0] = complex<FLOAT>(norm(x));
            Polynomial<complex<FLOAT > > poly(co);
            //cout << "poly = " << poly << endl;
            Q = Q / poly;
            Q1 = Q.derivative();
            //cout << "new Q = " << Q << endl;
            //cout << "check: Q * poly = " << Q * poly << endl;
            roots_found += 2;
         }
      }
      else
      {
         //std::cout << "Failed to find root, best x = " << x << ", value = " << s << std::endl;
         return roots_found;
      }
   }
   return P.deg();
}

template <class FLOAT>
int find_roots_over_C_q(const Polynomial<complex<FLOAT > >& P, std::vector<complex<FLOAT > >& roots)
{
   // Quote from H. Cohen "A Course in Computational Algebraic Number Theory", p. 146 :
   // "Given a squarefree polynomial P, this algorithm outputs its complex roots (in a random order).
   // In quite rare cases the algorithm may fail. On the other hand it is absolutely necessary that
   // the polynomial be squarefree (this can be achieved by replacing P by P/(P,P'))."

   const FLOAT ACCURACY = (FLOAT)1e-40;
   // Step 1. [Initializations]
   Polynomial<complex<FLOAT > > Q = P;
   Polynomial<complex<FLOAT > > P1 = P.derivative();
   Polynomial<complex<FLOAT > > Q1 = P1;
   Polynomial<complex<FLOAT > > Q2 = Q1.derivative();
   //cout << "find_roots_over_C of poly " << P << endl;
   //cout << "P'(X) = " << P1 << endl;
   int n = P.deg();
   // f = 1 if P has real coefficients, f = 0 otherwise
   int f = 1;
   for (int i = 0; i < n; i++)
   {
      complex<FLOAT > c = P.coefficient(i);
      if (!c.real()) f = 0;
   }

//   if (f) cout << "P has real coefficients" << endl;
//   else cout << "P has complex coefficients" << endl;

   while (n > 1)
   {
      // Step 2. [Initialize root finding]
      complex<FLOAT> x = complex<FLOAT>((FLOAT)1.3, (FLOAT)0.314159);
      complex<FLOAT> v = Q.evaluate(x);
      FLOAT m = norm(v);
      //cout << "Step 2., x = " << x << ", v = Q(x) = " << v << ", m = |v|^2 = " << v << endl;

      int done = 0;
      complex<FLOAT> lambda_c = Q1.evaluate(x);
      FLOAT lambda = FLOAT(2.0) * norm(lambda_c);
      lambda_c = Q.evaluate(x) * Q2.evaluate(x);
      lambda = lambda / sqrt(norm(lambda_c));
      //cout << "lambda = " << lambda << endl;
      int tries = 1000000;
      while (!done)
      {
         // Step 3. [Initialize recursion]
         int c = 0;
         complex<FLOAT> dx = v / Q1.evaluate(x);

         if (lambda > FLOAT(1.0)) lambda = FLOAT(1.0);
         dx = lambda * dx;

         //if (tries % 10000 == 0) cout << "Step 3., dx = v / Q'(x) = " << dx << endl;

         if (norm(dx) > ACCURACY)
         {
            int done1 = 0;
            while (!done1)
            {
               // Step 4. [Try a lambda]
               complex<FLOAT> y = x - dx;
               complex<FLOAT> v1 = Q.evaluate(y);
               FLOAT m1 = norm(v1);
               //if (tries % 10000 == 0) cout << "Step 4., y = " << y << ", v = " << v << ", v1 = " << v1 << ", m1 = " << m1 << ", m = " << m << ", m1 - m = " << m1 - m << endl;
               if (m1 < m && ((m - m1)/m > FLOAT(1e-10)))
               {
                  x = y;
                  v = v1;
                  m = m1;
                  done1 = 1;
               }
               else
               {
                  c++;
                  dx /= (FLOAT)4.0;
                  if (c >= 20)
                  {
                     // failed
                     FLOAT s(1.0);

                     int i = 0;
                     while (s > FLOAT(1e-10) && i < 100)
                     {
                        s = minimize_size_over_Re(Q, x);
                        s = minimize_size_over_Im(Q, x);
                        s = minimize_size(Q, x);
                        i++;
                     }
                     if (i >= 100)
                     {
                        std::cerr << "failed to find roots" << std::endl;
                        return 0;
                     }
                     done1 = 1;
                     done = 1;
                  }
               }
            }
         }
         else
         {
            done = 1;
         }
         tries--;
      }

      // Step 5. [Polish root]
      //cout << "Step 5., x = " << x << endl;
      x = x - P.evaluate(x) / P1.evaluate(x);
      x = x - P.evaluate(x) / P1.evaluate(x);
      //cout << "after Step 5., x = " << x << endl;

      // Step 6. [Divide]
      //cout << "Step 6." << endl;
      if (f == 0 || (f == 1 && fabs(x.imag()) < ACCURACY))
      {
         //cout << "one real root, x = " << x << endl;
         x = complex<FLOAT>(x.real(), FLOAT(0.0));
         roots.push_back(x);
         // Alternate method from Cohen p. 147, Remark (4)
         std::vector<complex<FLOAT> > a;
         int d = Q.deg();
         a.resize(Q.deg());
         a[d-1] = Q.coefficient(d);
         for (int i = d-1; i >= 1; --i)
         {
            a[i-1] = Q.coefficient(i) + x * a[i];
         }
         Polynomial<complex<FLOAT> > QQ(a);
         //cout << "QQ = " << QQ;
         //cout << endl;
         Q = QQ;
         n--;
      }
      else
      {
         //cout << "two conjugate roots, x = " << x << endl;
         roots.push_back(x);
         roots.push_back(conj(x));
         std::vector<complex<FLOAT > > co;
         co.resize(3);
         co[2] = complex<FLOAT>((FLOAT)1.0);
         co[1] = complex<FLOAT>((FLOAT)0.0) - (FLOAT)2.0 * x.real();
         co[0] = complex<FLOAT>(norm(x));
         Polynomial<complex<FLOAT > > poly(co);
         //cout << "poly = " << poly << endl;
         Q = Q / poly;
         //cout << "new Q = " << Q << endl;
         //cout << "check: Q * poly = " << Q * poly << endl;
         n -= 2;
      }
   }

   if (n == 1)
   {
      // just a linear equation to solve
      // Q(X) = z X + w = 0 => X = -w / z
      complex<FLOAT> x = complex<FLOAT>((FLOAT)0.0) - Q.coefficient(0) / Q.coefficient(1);
      roots.push_back(x);
   }


   return P.deg();
}

// Algorithm 3.4.2 (Squarefree Factorization)
// e.g. squarefree_factorisation<VeryLong, VeryLongModular>()
// or   squarefree_factorisation<long int, LongModular>()

template <class INTEGER> long int get_long(const INTEGER& i);
template <class INTEGER, class MODULAR_INTEGER> INTEGER get_integer(const MODULAR_INTEGER& i);
template <class INTEGER, class INTEGER2> INTEGER2 get_integer2(const INTEGER& i);

template <class INTEGER>
Polynomial<Quotient<INTEGER> > convert_to_quotient_field(const Polynomial<INTEGER>& zpoly)
{
   std::vector<Quotient<INTEGER> > qcoefficients;
   qcoefficients.resize(zpoly.deg() + 1);
   for (int i = 0; i <= zpoly.deg(); i++)
   {
      qcoefficients[i] = Quotient<INTEGER>(zpoly.coefficient(i));
   }
   return Polynomial<Quotient<INTEGER> >(qcoefficients);
}

template <class INTEGER, class INTEGER2, class MODULAR_INTEGER>
Polynomial<MODULAR_INTEGER> convert_to_F_p(const Polynomial<INTEGER>& zpoly, const INTEGER2& modulus)
{
   MODULAR_INTEGER::set_default_modulus(modulus);
   std::vector<MODULAR_INTEGER > fplcoefficients;
   fplcoefficients.resize(zpoly.deg() + 1);

   int k = 0;
   for (k = 0; k < zpoly.deg() + 1; k++)
   {
      INTEGER2 vl = get_integer2<INTEGER, INTEGER2>(zpoly.coefficient(k) % modulus);
      fplcoefficients[k] = MODULAR_INTEGER(vl);
   }
   Polynomial<MODULAR_INTEGER> fppoly(fplcoefficients);
   return fppoly;
}

template <class INTEGER, class MODULAR_INTEGER>
Polynomial<MODULAR_INTEGER> powmodf_vl(const INTEGER& p,
                                       const Polynomial<MODULAR_INTEGER>& g,
                                       const Polynomial<MODULAR_INTEGER>& f0)
{
   const MODULAR_INTEGER one(1L);
   Polynomial<MODULAR_INTEGER> f = f0;

   // xi[0] (X) = X
   Polynomial<MODULAR_INTEGER> xi[2];
   xi[0] = g;

   // power[0] (X) = 1
   Polynomial<MODULAR_INTEGER> power[2];
   power[0] = Polynomial<MODULAR_INTEGER>(one);

   int from = 0;
   int to = 1;
   int powfrom = 0;
   int powto = 1;
   INTEGER tmp = p;

   INTEGER i = 1L;

   for (i = 1L; i <= p; i *= 2L)
   {
      //cout << "i = " << i << endl;
      //cout << "tmp = " << tmp << endl;
      if (tmp % 2L == 1L)
      {
         //cout << "i = " << i << endl;
         power[powto] = (power[powfrom] * xi[from]) % f;
         powfrom = 1 - powfrom;
         powto = 1 - powto;
      }

      //cout << "xi[" << from << "] = " << xi[from] << endl;
      //cout << "f = " << f << endl;
      xi[to] = (xi[from] * xi[from]) % f;

      from = 1 - from;
      to = 1 - to;
      tmp /= 2L;
   }
   return power[powfrom];

}

template <class INTEGER, class MODULAR_INTEGER>
Polynomial<MODULAR_INTEGER> powpowmodf(const INTEGER& p, long int d,
                                       const Polynomial<MODULAR_INTEGER>& T,
                                       const Polynomial<MODULAR_INTEGER>& A)
{
   Polynomial<MODULAR_INTEGER> Tpower = T;
   Polynomial<MODULAR_INTEGER> product = T;
   for (long int i = 1; i < d; ++i)
   {
      Tpower = powmodf_vl(p, Tpower, A);
      product = (product * Tpower) % A;
   }
   INTEGER p1 = p;
   p1 -= 1L;
   p1 /= 2L;

   Polynomial<MODULAR_INTEGER> result = powmodf_vl<INTEGER, MODULAR_INTEGER>(p1, product, A);
   return result;
}

template <class INTEGER, class MODULAR_INTEGER>
void squarefree_factorisation(const Polynomial<MODULAR_INTEGER>& A,
                              const INTEGER& p,
                              std::vector<std::pair<int, Polynomial<MODULAR_INTEGER> > >& Ai)
{
   // Step 1. [Initialize]
   int e = 1;
   Polynomial<MODULAR_INTEGER> T0 = A;
   Polynomial<MODULAR_INTEGER> T;
   Polynomial<MODULAR_INTEGER> V;
   long int pl = get_long<INTEGER>(p);

   bool skip_step_2 = false;
   int k = 0;

   while (1)
   {
      int done = 0;
      while (!done)
      {
         if (!skip_step_2)
         {
            // Step 2. [Initialise e-loop]
//            cout << "Step 2 : T0 = " << T0 << endl;
            if (T0.deg() == 0) return;
            T = gcd(T0, T0.derivative());
            V = T0 / T;
//            cout << "T = " << T << endl;
//            cout << "V = " << V << endl;
            k = 0;
         }
         skip_step_2 = false;
         // Step 3. [Finished e-loop?]
//         cout << "Step 3 : V = " << V << endl;
         if (V.deg() == 0)
         {
            std::vector<MODULAR_INTEGER> c;
            int dd = (T.deg() / pl);
            c.resize(dd + 1);
            for (int i = 0; i < dd + 1; i++)
            {
               c[i] = T.coefficient(i * pl);
            }
            T0 = Polynomial<MODULAR_INTEGER>(c);
            e *= pl;
//            cout << "new T0 = " << T0 << endl;
//            cout << "new e = " << e << endl;
         }
         else done = 1;
      }

      // Step 4. [Special case]
      k++;
      if (k % pl == 0)
      {
//        cout << "Step 4 : T = " << T << ", V = " << V << endl;
         T = T/V;
//        cout << "Step 4 after : T = " << T << ", V = " << V << endl;
         k++;
      }

      // Step 5. [Compute A_ek]
      Polynomial<MODULAR_INTEGER> W = gcd(T, V);
//      cout << "Step 5 : T = " << T << ", V = " << V << ", W = " << W << endl;
      Polynomial<MODULAR_INTEGER> Aek = V / W;
      V = W;
      T = T / V;
//      cout << "Step 5 after : T = " << T << ", V = " << V << ", W = " << W << endl;
      if (Aek.deg() != 0)
      {
         Aek.make_monic();
         Ai.push_back(std::pair<int, Polynomial<MODULAR_INTEGER> >(e * k, Aek));
      }
      skip_step_2 = true;
   }
}

// Algorithm 3.4.3 (Distinct Degree Factorization)
template <class INTEGER, class MODULAR_INTEGER>
void distinct_degree_factorisation(const Polynomial<MODULAR_INTEGER>& A,
                                   const INTEGER& p,
                                   std::vector<Polynomial<MODULAR_INTEGER > >& Ad)
{
   // Given a squarefree polynomial A in F_p[X], this algorithm finds for each d the polynomial Ad which is the
   // product of the irreducible factors A of degree d.
   // Step 1. [Initialise]
   Polynomial<MODULAR_INTEGER> V = A;
   const MODULAR_INTEGER zero(0L);
   const MODULAR_INTEGER one(1L);
   std::vector<MODULAR_INTEGER> x;
   x.resize(2);
   x[0] = zero;
   x[1] = one;
   Polynomial<MODULAR_INTEGER> X(x);
   Polynomial<MODULAR_INTEGER> W = X;
   int d = 0;
   int e = V.deg();
   Ad.resize(e);

   while (1)
   {
      // Step 2. [Finished]
      //cout << "Step 2 : d = " << d << ", e = " << e << ", V = " << V << ", W = " << W << endl;
      e = V.deg();
      if (d + 1 > e / 2)
      {
         if (e > 0)
         {
            //cout << "Finished : d = " << d << ", e = " << e << ", V = " << V << ", W = " << W << endl;
            Ad[e - 1] = V;
            Ad[e - 1].make_monic();
            //cout << "Step 2 : Ad[" << e << "] = " << Ad[e - 1] << endl;
            for (int i = d + 1; i < e; i++)
            {
               Ad[i - 1] = Polynomial<MODULAR_INTEGER>(one);
            }
         }
         return;
      }
      else
      {
         d++;
         // W = W^p mod V
         W = powmodf_vl(p, W, V);
         //cout << "W^" << p << " mod (" << V << ") = " << W << endl;
      }

      // Step 3. [Output Ad]
      Ad[d - 1] = gcd(W - X, V);
      Ad[d - 1].make_monic();

      //cout << "Step 3 : d = " << d << ", e = " << e << ", V = " << V << ", W = " << W << endl;
      //cout << "Step 3 : Ad[" << d << "] = " << Ad[d - 1] << endl;
      if (Ad[d - 1].deg() != 0 || Ad[d - 1].coefficient(0) != 1L)
      {
         V = V / Ad[d - 1];
         W = W % V;
         //cout << "new V = " << V << ", new W = " << W << endl;
      }
   }
}

// Algorithm 3.4.8 (Split for p = 2)
template <class MODULAR_INTEGER>
void split2(const Polynomial<MODULAR_INTEGER>& A, int d, std::vector<Polynomial<MODULAR_INTEGER > >& factors)
{
   // Step 1. [Initialize]
   int k = A.deg() / d;
   if (k == 1)
   {
//      cout << "Irreducible factor = " << A << endl;
      factors.push_back(A);
      return;
   }

   std::vector<MODULAR_INTEGER> c;
   const MODULAR_INTEGER zero(0L);
   const MODULAR_INTEGER one(1L);
   c.resize(2);
   c[0] = zero;
   c[1] = one;
   Polynomial<MODULAR_INTEGER> X(c);
   Polynomial<MODULAR_INTEGER> T = X;

   Polynomial<MODULAR_INTEGER> B;
   int done = 0;
   while (!done)
   {
      // Step 2. [Test T]
      Polynomial<MODULAR_INTEGER> C = T;
      for (int i = 0; i < d - 1; i++)
      {
         C = (T + C * C) % A;
      }
      B = gcd(A,C);
      if (B.deg() == 0 || B.deg() == A.deg())
      {
         T = T * X * X;
      }
      else done = 1;
   }

   // Step 3. [Recurse]
   split2(B, d, factors);
   split2(A/B, d, factors);
}

// Algorithm 3.4.6 (Cantor-Zassenhaus Split)
template <class INTEGER, class MODULAR_INTEGER>
void CZ_split(const Polynomial<MODULAR_INTEGER>& A, const INTEGER& p, int d, std::vector<Polynomial<MODULAR_INTEGER > >& factors)
{
   //cout << "CZ_split : A = " << A << ", p = " << p << ", d = " << d << endl;
   // Step 1. [Initialize]
   int k = A.deg() / d;
   if (k == 1)
   {
      //cout << "Irreducible factor = " << A << endl;
      factors.push_back(A);
      return;
   }

   int done = 0;
   Polynomial<MODULAR_INTEGER> B;
   const MODULAR_INTEGER one(1L);
   while (!done)
   {
      // Step 2. [Try a T]
      // Choose T such that T is monic of degree less than or equal to 2d - 1
      int Td = (genrand() % (long)(2 * d - 1)) + 1;
      std::vector<MODULAR_INTEGER> Tc;
      Tc.resize(Td + 1);
      Tc[Td] = one;
      for (int i = 0; i < Td; i++)
      {
         Tc[i] = static_cast<MODULAR_INTEGER>(static_cast<INTEGER>(genrand()) % p);
      }
      Polynomial<MODULAR_INTEGER> T(Tc);
#if 0
      INTEGER power = 1L;
      for (int j = 0; j < d; j++) power *= p;
      power -= 1L;
      power /= 2L;

      //cout << "CZ_split : T = " << T << ", A = " << A << ", power = " << power << endl;
      B = powmodf_vl(power, T, A);
#endif
      B = powpowmodf(p, d, T, A);
      //cout << "B = " << B << endl;
      B.set_coefficient(0, B.coefficient(0) - one);
      B = gcd(A, B);
      //cout << "after gcd, B = " << B << endl;
      if (B.deg() != 0 && B.deg() != A.deg()) done = 1;
   }

   CZ_split(B, p, d, factors);
   CZ_split(A/B, p, d, factors);
}

template <class INTEGER, class MODULAR_INTEGER>
void final_split(const Polynomial<MODULAR_INTEGER>& A,
                 const INTEGER& p,
                 int d,
                 std::vector<Polynomial<MODULAR_INTEGER > >& factors)
{
   if (p > static_cast<INTEGER>(2L))
   {
      CZ_split(A, p, d, factors);
   }
   else
   {
      split2(A, d, factors);
   }
}

// Algorithm 3.4.1 (Factor in F_p[X])
template <class INTEGER, class INTEGER2, class MODULAR_INTEGER>
void factor_over_F_p(const Polynomial<INTEGER>& poly,
                     const INTEGER2& p,
                     std::vector<Polynomial<MODULAR_INTEGER > >& factors)
{
   MODULAR_INTEGER::set_default_modulus(p);
   // remember l(poly);
   INTEGER l_poly = poly.coefficient(poly.deg());
   MODULAR_INTEGER l_poly_p(get_integer2<INTEGER, INTEGER2>(l_poly));
   Polynomial<MODULAR_INTEGER> A = convert_to_F_p<INTEGER, INTEGER2, MODULAR_INTEGER>(poly, p);
   A.make_monic();
   // Step 1. [Squarefree factorization]
   // Find polynomials A_1, A_2, ... , A_k in F_p[X] s.t.
   // (1) A = A_1 * A_2^2 * ... * A_k^k
   // (2) The A_i are squarefree and coprime

   std::vector<std::pair<int, Polynomial<MODULAR_INTEGER> > > Ai;
   //cout << "A = " << A << endl;
   //cout << "p = " << p << endl;
   squarefree_factorisation(A, p, Ai);
   //cout << "After squarefree_factorisation() : " << endl;
   for (size_t i = 0; i < Ai.size(); i++)
   {
      //cout << "(" << Ai[i].first << ", " << Ai[i].second << ")" << endl;
   }

   // Step 2. [Distinct degree factorization]
   const MODULAR_INTEGER one(1L);
   std::vector<std::vector<Polynomial<MODULAR_INTEGER> > > Ad;
   Ad.resize(Ai.size());
   for (size_t i = 0; i < Ai.size(); i++)
   {
      distinct_degree_factorisation(Ai[i].second, p, Ad[i]);
      Polynomial<MODULAR_INTEGER> check(one);
      //cout << "result of distinct_degree_factorisation for i = " << i << endl;
      for (size_t d = 0; d < Ad[i].size(); d++)
      {
         //cout << "Ad[" << i << "][" << d << "] = " << Ad[i][d] << endl;
         if (Ad[i][d].deg() >= 0) check *= Ad[i][d];
      }
      check.make_monic();
      //cout << "check = " << check << endl;
   }

   // Step 3. [Final Splittings]
   for (size_t i = 0; i < Ad.size(); i++)
   {
      std::vector<Polynomial<MODULAR_INTEGER > > irreducible_factors;
      for (size_t d = 0; d < Ad[i].size(); d++)
      {
         if (Ad[i][d].deg() >= static_cast<int>(d + 1))
         {
            final_split(Ad[i][d], p, static_cast<int>(d + 1), irreducible_factors);
         }
      }
      int multiplicity = Ai[i].first;
      for (int j = 0; j < multiplicity; j++)
      {
         for (size_t k = 0; k < irreducible_factors.size(); k++)
         {
            factors.push_back(irreducible_factors[k]);
         }
      }
   }

   // Step 4. [Cleanup]
   // Group together all the identical factors found, order them by degree,
   // output the complete factorisation and terminate the algorithm
   std::sort(factors.begin(), factors.end());
   if (factors.size() > 0)
   {
      factors[0] *= l_poly_p;
   }
}

template <class INTEGER, class MODULAR_INTEGER>
Polynomial<INTEGER> lift(const Polynomial<MODULAR_INTEGER>& fp)
{
   std::vector<INTEGER> c;
   c.resize(fp.deg() + 1);
   for (int i = 0; i <= fp.deg(); i++)
   {
      MODULAR_INTEGER co = fp.coefficient(i);
      c[i] = get_integer<INTEGER, MODULAR_INTEGER>(co);
   }
   return Polynomial<INTEGER>(c);
}

template <class INTEGER, class MODULAR_INTEGER>
Polynomial<INTEGER> monic_lift(const Polynomial<MODULAR_INTEGER>& fp)
{
   std::vector<INTEGER> c;
   Polynomial<MODULAR_INTEGER> fp_tmp = fp;
   fp_tmp.make_monic();
   c.resize(fp.deg() + 1);
   for (int i = 0; i <= fp.deg(); i++)
   {
      MODULAR_INTEGER co = fp_tmp.coefficient(i);
      c[i] = get_integer<INTEGER, MODULAR_INTEGER>(co);
   }
   return Polynomial<INTEGER>(c);
}

// NOTE: Irreducible polynomials over F_p
// A polynomial A in F_p[X] of degree n is irreducible iff
// X^p^n = X mod A(X)
// and for each prime q | n
// gcd(X^p^(n/q) - X, A(X)) = 1

template <class F>
void check_root(const Polynomial<F>& f, F& root, int id)
{
   if (f.evaluate(root) != 0L)
   {
      std::cout << id << ". Bad root : " << root << std::endl;
   }
}

// Algorithm borrowed from Matthew Kudzin with thanks
// quick calculation of x^p mod f(x)
template <class INTEGER, class MODULAR_INTEGER>
Polynomial<MODULAR_INTEGER> xtopmodf(INTEGER p, const Polynomial<MODULAR_INTEGER>& f0)
{
   const MODULAR_INTEGER one(1L);
   long int p_li = get_long(p);
   if (p_li < (long)f0.deg())
   {
      std::vector<MODULAR_INTEGER> c;
      c.resize(p_li+1);
      c[p_li] = one;
      return Polynomial<MODULAR_INTEGER>(c);
   }
   Polynomial<MODULAR_INTEGER> f = f0;

   // xi[0] (X) = X
   Polynomial<MODULAR_INTEGER> xi[2];
   std::vector<MODULAR_INTEGER> c;
   c.clear();
   c.resize(2);
   c[1] = one;
   xi[0] = Polynomial<MODULAR_INTEGER>(c);


   // power[0] (X) = 1
   Polynomial<MODULAR_INTEGER> power[2];
   power[0] = Polynomial<MODULAR_INTEGER>(1L);


   int from = 0;
   int to = 1;
   int powfrom = 0;
   int powto = 1;

   int i = 0;

   for (i = 0; (1<<i) <= p_li; i++)
   {
      if (p_li & (1<<i))
      {
         //cout << "1<<i = " << (1<<i) << endl;
         power[powto] = (power[powfrom] * xi[from]) % f;
         powfrom = 1 - powfrom;
         powto = 1 - powto;
      }

      //cout << "xi[" << from << "] = " << xi[from] << endl;
      //cout << "f = " << f << endl;
      xi[to] = (xi[from] * xi[from]) % f;

      from = 1 - from;
      to = 1 - to;
   }
   return power[powfrom];
}

template <class INTEGER, class INTEGER2, class MODULAR_INTEGER>
void find_roots_mod_p_1(const Polynomial<MODULAR_INTEGER>& f,
                        const Polynomial<MODULAR_INTEGER>& g0,
                        MODULAR_INTEGER b,
                        Polynomial<MODULAR_INTEGER>& X_minus_1,
                        INTEGER2 p,
                        std::vector<MODULAR_INTEGER>& roots)
{
   Polynomial<MODULAR_INTEGER> g = g0;

   const MODULAR_INTEGER zero(0L);
   while (g.evaluate(zero) == zero)
   {
      //check_root(f,b,2);
      roots.push_back(b);
      if (g.deg() == 1) return;
      g.divide_by_X();
   }
   if (g.deg() == 1)
   {
      MODULAR_INTEGER root = g.coefficient(0) * (-1L) + b;
      check_root(f,root,10);
      roots.push_back(root);
      return;
   }

   // fast x^(p-1)/2 mod g(x)
   Polynomial<MODULAR_INTEGER> hp1 = xtopmodf<INTEGER2, MODULAR_INTEGER>((p-1L)/2L, g);
   Polynomial<MODULAR_INTEGER> hp2 = hp1 + Polynomial<MODULAR_INTEGER>(1L);
   hp1 -= Polynomial<MODULAR_INTEGER>(1L);

   Polynomial<MODULAR_INTEGER> g1;
   if (hp1 == 0L) g1 = g;
//   else if (hp1.deg() == 0) g1 = hp1;
   else g1 = gcd(g, hp1);
   g1.make_monic();
   Polynomial<MODULAR_INTEGER> g2;
   if (hp2 == 0L) g2 = g;
//   else if (hp2.deg() == 0) g2 = hp2;
   else g2 = gcd(g, hp2);
   g2.make_monic();

   // cout << "g0 = " << g0 << endl;
   // cout << "g = " << g << endl;
   // cout << "g1 = " << g1 << endl;
   // cout << "g2 = " << g2 << endl;
   // cout << "b = " << b << endl;

   if (g1.deg() == 0 || g2.deg() == 0)
   {
      Polynomial<MODULAR_INTEGER> gg = g.evaluate(X_minus_1);
      find_roots_mod_p_1<INTEGER, INTEGER2, MODULAR_INTEGER>(f, gg, b - 1L, X_minus_1, p, roots);
   }
   else
   {
      if (g1.deg() == 1)
      {
         // found a root
         MODULAR_INTEGER root = g1.coefficient(0) * (-1L) + b;
         //check_root(f,root,3);
         roots.push_back(root);
      }
      else
      {
         Polynomial<MODULAR_INTEGER> gg = g1.evaluate(X_minus_1);
         find_roots_mod_p_1<INTEGER, INTEGER2, MODULAR_INTEGER>(f, gg, b - 1L, X_minus_1, p, roots);
      }

      if (g2.deg() == 1)
      {
         // found a root
         MODULAR_INTEGER root = g2.coefficient(0) * (-1L) + b;
         //check_root(f, root, 4);
         roots.push_back(root);
      }
      else
      {
         Polynomial<MODULAR_INTEGER> gg = g2.evaluate(X_minus_1);
         find_roots_mod_p_1<INTEGER, INTEGER2, MODULAR_INTEGER>(f, gg, b - 1L, X_minus_1, p, roots);
      }
   }

   return;
}

template <class INTEGER, class INTEGER2, class MODULAR_INTEGER>
void find_roots_mod_p(const Polynomial<INTEGER>& poly, INTEGER2 p, std::vector<MODULAR_INTEGER>& roots)
{
   MODULAR_INTEGER::set_default_modulus(p);
   const MODULAR_INTEGER one(1L);
   //cout << "Finding roots of " << poly << " mod " << p << " ..." << endl;
   //cout << "one = " << one << endl;
   Polynomial<MODULAR_INTEGER> f = convert_to_F_p<INTEGER, INTEGER2, MODULAR_INTEGER>(poly, p);
   f.make_monic();
   //cout << f << endl;

   std::vector<MODULAR_INTEGER> fpc;
   fpc.resize(2);

   // Construct X
   fpc[1] = MODULAR_INTEGER(1L);
   Polynomial<MODULAR_INTEGER> X(fpc);

   // calculate first gcd using xtopmodf
   Polynomial<MODULAR_INTEGER> total_product_mod_f = xtopmodf<INTEGER2, MODULAR_INTEGER>(p, f); // x^p mod f(x)
   total_product_mod_f -= X; // x^p - x mod f(x)
   total_product_mod_f.make_monic();

   //cout << "total_product_mod_f = " << total_product_mod_f << endl;
   //cout << "f = " << f << endl;

   Polynomial<MODULAR_INTEGER> g;
   if (total_product_mod_f == 0L) g = f;
   else g = gcd(f, total_product_mod_f);
   g.make_monic();
   //cout << "g = " << g << endl;

   if (g.deg() == 0)
   {
      return;
   }

   if (g.deg() == 1)
   {
      // found a root
      MODULAR_INTEGER root = g.coefficient(0) * (-1L);
      roots.push_back(root);
      return;
   }

   // Construct X - b
   MODULAR_INTEGER b = 0L;
   std::vector<MODULAR_INTEGER> fpc2;
   fpc2.resize(2);
   fpc2[0] = p - one;
   fpc2[1] = one;
   Polynomial<MODULAR_INTEGER> X_minus_1(fpc2);

   find_roots_mod_p_1<INTEGER, INTEGER2, MODULAR_INTEGER>(f, g, b, X_minus_1, p, roots);
}

template <class INTEGER, class INTEGER2, class MODULAR_INTEGER>
long int count_roots_mod_p(const Polynomial<INTEGER>& poly, const INTEGER2 p)
{
   MODULAR_INTEGER::set_default_modulus(p);
   Polynomial<MODULAR_INTEGER> f = convert_to_F_p<INTEGER, INTEGER2, MODULAR_INTEGER>(poly, p);
   f.make_monic();

   //cout << "f = " << f << endl;
   std::vector<MODULAR_INTEGER > fpc;
   fpc.resize(2);

   // Construct X
   const MODULAR_INTEGER one(1L);
   fpc[1] = one;
   Polynomial<MODULAR_INTEGER> X(fpc);

   // calculate first gcd using xtopmodf
   Polynomial<MODULAR_INTEGER> total_product_mod_f = xtopmodf<INTEGER2, MODULAR_INTEGER>(p, f); // x^p mod f(x)
   total_product_mod_f -= X; // x^p - x mod f(x)
   total_product_mod_f.make_monic();
   //cout << "total_product_mod_f = " << total_product_mod_f << endl;

   Polynomial<MODULAR_INTEGER> g;
   if (total_product_mod_f == 0L) return f.deg();

   g = gcd(f, total_product_mod_f);
   //cout << "g = " << g << endl;
   return g.deg();
}

#if 0
template <class INTEGER, class INTEGER2, class MODULAR_INTEGER>
void find_roots_mod_p_k(const Polynomial<INTEGER>& f, INTEGER2& p, long int k, std::vector<INTEGER>& roots)
{
   // find roots of f(X) = 0 mod p^k
   // by first solving f(X) = 0 mod p and then
   // lifting the roots to mod p^k

   std::vector<MODULAR_INTEGER> roots_mod_p;
   find_roots_mod_p<INTEGER, INTEGER2, MODULAR_INTEGER>(f, p, roots_mod_p);
   if (roots_mod_p.size() == 0)
   {
      roots.resize(0);
      return;
   }

   roots.resize(roots_mod_p.size());

   // get derivative of f
   Polynomial<INTEGER> f1 = f.derivative();

   const MODULAR_INTEGER one(1L);
   for (size_t r = 0; r < roots_mod_p.size(); r++)
   {
      INTEGER c0 = get_integer<INTEGER, MODULAR_INTEGER>(roots_mod_p[r]);
      MODULAR_INTEGER::set_default_modulus(p);
      // calculate inverse of f'(c0) mod p
      MODULAR_INTEGER t_vlm = one / MODULAR_INTEGER(f1.evaluate(c0));
      INTEGER t = get_integer<INTEGER, MODULAR_INTEGER>(t_vlm);
      INTEGER c = c0;
      INTEGER p_k = p;
      for (int power = 2; power <= k; power++)
      {
         c += ((p - (f.evaluate(c) / p_k) * t) % p) * p_k;
         p_k *= p;
         // sanity check
         if (f.evaluate(c) % p_k != 0L)
         {
            std::cout << "Error" << std::endl;
            std::cout << "c = " << c << std::endl;
            std::cout << "f(c) mod p^k = " << f.evaluate(c) % p_k << std::endl;
         }
      }
      // check
      //cout << "c = " << c << endl;
      //cout << "f(c) mod p^k = " << f.evaluate(c) % p_k << endl;
      roots[r] = c;
   }

}
#endif

// Algorithm 3.1.1 (Euclidean Division)
template <class F>
void euclidean_division(const Polynomial<F>& A,
                        const Polynomial<F>& B,
                        Polynomial<F>& Q,
                        Polynomial<F>& R)
{
   // Step 1. [Initialize]
   R = A;
   Q = Polynomial<F>(F(0L));

   // Step 2. [Finished]
   while (R.deg() >= B.deg())
   {
      // Step 3. [Find coefficient]
      int d = R.deg() - B.deg();
      std::vector<F> c;
      c.resize(d + 1);
      c[d] = R.coefficient(R.deg()) / B.coefficient(B.deg());
      Polynomial<F> S(c);
      Q = Q + S;
      R = R - S * B;
   }
}

template <class F>
int divides(const Polynomial<F>& f, const Polynomial<F>& p)
{
   if ((p / f) * f == p) return 1;
   return 0;
}

template <class F>
int exponent(const Polynomial<F>& p, const Polynomial<F>& f)
{
   int e = 0;
   Polynomial<F> pp = p;
   while (divides(f, pp))
   {
      pp = pp / f;
      e++;
   }
   return e;
}
#endif
