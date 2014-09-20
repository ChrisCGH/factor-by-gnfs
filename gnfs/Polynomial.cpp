#include "Polynomial.h"
#include "VeryLong.h"
#include "LongModular.h"
#include "VeryLongModular.h"
#include "Polynomial.inl"
#include <algorithm>
#include <deque>
#include <ctype.h>
#include "MPFloat.h"
#if 1
template <> long int get_long(const VeryLong& i)
{
   return i.get_long();
}
template <> long int get_long(const long int& i)
{
   return i;
}
template <> long int get_long(const long long int& i)
{
   return static_cast<long int>(i);
}
template <> long int get_integer(const LongModular& i)
{
   return i.get_long();
}
template <> long long int get_integer(const LongModular& i)
{
   return i.get_long();
}
template <> VeryLong get_integer(const LongModular& i)
{
   return i.get_long();
}
template <> VeryLong get_integer(const VeryLongModular& i)
{
   return i.get_very_long();
}
template <> VeryLong get_integer2(const VeryLong& i)
{
   return i;
}
template <> long int get_integer2(const VeryLong& i)
{
   return i.get_long();
}
template <> long long int get_integer2(const VeryLong& i)
{
   return i.get_long();
}
#endif
Polynomial<VeryLong> read_skewed_polynomial(std::istream& is, VeryLong& m, long int& s)
{
   // reads a skewed polynomial from the input stream
   // assumes format
   // 7111977108918472837611244937303928 - 2523764323493368834122550403406 X - 100831972740733548080518208
   // X^2 + 7467483273731340429431 X^3 + 150036498813232980 X^4 + 2859632424000 X^5
   // m = 20651964074656838659712611715
   // s = 18428
   std::vector<VeryLong> coeff;
   coeff.resize(6);
   int done = 0;
   std::string poly_str;
   std::string str;
   while (!done)
   {
      int poly_done = 0;
      while (!poly_done)
      {
         getline(is, str);
         if (str.c_str()[0] == 'm') poly_done = 1;
         else
         {
            poly_str = str;
         }
      }
      m = str.substr(4);
      getline(is, str);
      //s = atol(&buf[4]);
      s = std::atol(str.substr(4).c_str());
      done = 1;
   }

   char tmp[1024];
   const char* c = poly_str.c_str();
   char* d = tmp;
   int coeff_index = 0;
   while (*c)
   {
      while (*c != ' ')
      {
         *d = *c;
         c++;
         d++;
      }
      *d = '\0';

      coeff[coeff_index] = tmp;
      coeff_index++;
      c++; // skip space
      if (*c == 'X')
      {
         while (*c && *c != ' ') c++; // skip to space before sign or to end of string
         if (*c) c++;
      }

      if (*c)
      {
         d = tmp;
         if (*c == '-')
         {
            *d = *c;
            d++;
         }
         c++; // skip space
         c++;
      }
   }


   return Polynomial<VeryLong>(coeff);
}

template <> Polynomial<VeryLong> Polynomial<VeryLong>::read_polynomial(std::istream& is)
{
   // reads a polynomial from the input stream
   // assumes format
   // 7111977108918472837611244937303928 - 2523764323493368834122550403406 X - 100831972740733548080518208
   // X^2 + 7467483273731340429431 X^3 + 150036498813232980 X^4 + 2859632424000 X^5
   std::string poly_str;
   getline(is, poly_str);
   if (poly_str.find("DEGREE = ") == 0)
   {
      std::string::size_type eqpos = poly_str.find('=');
      if (std::string::npos != eqpos)
      {
         if (eqpos + 2 < poly_str.size())
         {
            std::string s = poly_str.substr(eqpos + 2);
            int degree = std::atoi(s.c_str());
            std::vector<VeryLong> coeff;
            for (int i = 0; i <= degree; ++i)
            {
               getline(is, s);
               VeryLong c(s.c_str());
               coeff.push_back(c);
            }
            return Polynomial<VeryLong>(coeff);
         }
      }
   }
   return Polynomial<VeryLong>::read_polynomial(poly_str.c_str());
}

// Algorithm 3.3.1 (Sub-Resultant GCD)
Polynomial<VeryLong> sub_resultant_GCD(const Polynomial<VeryLong>& AA, const Polynomial<VeryLong>& BB)
{
   Polynomial<VeryLong> A = AA;
   Polynomial<VeryLong> B = BB;
   // Step 1. [Initializations and reductions]
   if (B.deg() > A.deg())
   {
      A = BB;
      B = AA;
   }
   const VeryLong zero(0L);
   if (B == Polynomial<VeryLong>(zero)) return A;

   VeryLong a = A.content();
   VeryLong b = B.content();
   VeryLong d = gcd<VeryLong>(a, b);
   A = A / a;
   B = B / b;
   VeryLong g(1L);
   VeryLong h(1L);

   while (1)
   {
      // Step 2. [Pseudo division]
      int delta = A.deg() - B.deg();
      Polynomial<VeryLong> Q;
      Polynomial<VeryLong> R;
      pseudo_divide(A, B, Q, R);
      if (R == Polynomial<VeryLong>(zero))
      {
         Polynomial<VeryLong> res = B / B.content();
         res = d * res;
         return res;
      }
      if (R.deg() == 0)
      {
         Polynomial<VeryLong> res(d);
         return res;
      }
      A = B;
      VeryLong e = g * pow<VeryLong, int>(h, delta);
      B = R / e;
      g = A.coefficient(A.deg());
      VeryLong f = pow<VeryLong, int>(h, delta - 1);
      h = pow<VeryLong, int>(g, delta) / f;
   }
}

//Algorithm 3.5.5 (Hensel Lift)
void hensel_lift(const VeryLong& p,
                 const VeryLong& q,
                 const Polynomial<VeryLong>& A,
                 const Polynomial<VeryLong>& B,
                 const Polynomial<VeryLong>& C,
                 const Polynomial<VeryLong>& U,
                 const Polynomial<VeryLong>& V,
                 Polynomial<VeryLong>& A1,
                 Polynomial<VeryLong>& B1)
{
   // AB = C mod q
   // AU + BV = 1 mod q
   //cout << "hensel_lift:" << endl;
   //cout << "p = " << p << ", q = " << q << endl;
   //cout << "A = " << A << ", B = " << B << ", C = " << C << endl;
   //cout << "U = " << U << ", V = " << V << endl;
   Polynomial<VeryLong> check3 = C - A * B;
   Polynomial<VeryLongModular> check4 = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(check3, q);
   if (check4 != Polynomial<VeryLongModular>(VeryLongModular(0L)))
   {
      std::cout << "Problem: C - A * B != 0 mod " << q << std::endl;
      std::cout << "check4 = " << check4 << std::endl;
   }
   //cout << "C - AB = " << C - A * B << endl;
   VeryLong r = gcd<VeryLong>(p, q);
   // Step 1. [Euclidean division]
   Polynomial<VeryLong> f1 = (C - A * B) / q;
   //cout << "f1 = " << f1 << endl;
   VeryLongModular::set_default_modulus(r);
   Polynomial<VeryLongModular> f = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(f1, r);
   Polynomial<VeryLongModular> V_ = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(V, r);
   Polynomial<VeryLongModular> A_ = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(A, r);
   //cout << "f = " << f << ", V_ = " << V_ << ", A_ = " << A_ << endl;
   Polynomial<VeryLongModular> t;
   Polynomial<VeryLongModular> R;
   euclidean_division<VeryLongModular>(V_ * f, A_, t, R);
   //cout << "t = " << t << endl;
   //cout << "R = " << R << endl;

   // check
   if (R.deg() >= A.deg())
   {
      std::cout << "Problem: euclidean_division of " << V_ * f << " and " << A_ << std::endl;
      std::cout << "gave t = " << t << ", R = " << R << std::endl;
   }
   Polynomial<VeryLongModular> check = V_ * f - A_ * t;
   if (check != R)
   {
      std::cout << "Problem: check != R, check = " << check << ", R = " << R << std::endl;
   }

   // Step 2. [Terminate]
   Polynomial<VeryLong> A0 = lift<VeryLong, VeryLongModular>(R);
   Polynomial<VeryLongModular> U_ = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(U, r);
   Polynomial<VeryLongModular> B_ = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(B, r);
   Polynomial<VeryLong> B0 = lift<VeryLong, VeryLongModular>(U_ * f + B_ * t);
   //cout << "A0 = " << A0 << endl;
   //cout << "B0 = " << B0 << endl;

   A1 = A + q * A0;
   B1 = B + q * B0;
   //cout << "A1 = " << A1 << endl;
   //cout << "B1 = " << B1 << endl;
   // should have C = A1 * B1 mod qr
   Polynomial<VeryLong> check1 = A1 * B1;
   VeryLongModular::set_default_modulus(q * r);
   Polynomial<VeryLongModular> C2 = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(C, q * r);
   Polynomial<VeryLongModular> check2 = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(check1, q * r);
   if (C2 != check2)
   {
      std::cout << "Problem in result of hensel lift : " << std::endl;
      std::cout << "qr = " << q * r << std::endl;
      std::cout << "C mod qr = " << C2 << std::endl;
      std::cout << "A1 * B1 mod qr = " << check2 << std::endl;
   }
}

void hensel_lift(const VeryLong& p,
                 const VeryLong& q,
                 const Polynomial<VeryLong>& C,
                 std::deque<Polynomial<VeryLong> >& VV,
                 std::deque<Polynomial<VeryLong> >& W)
{
   //cout << "generalised hensel lift : " << endl;
   Polynomial<VeryLong> A = VV[0];
   Polynomial<VeryLong> B = VV[1];

   for (size_t i = 2; i < VV.size(); i++) B *= VV[i];

   VeryLongModular::set_default_modulus(p);
   Polynomial<VeryLongModular> Ap_ = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(A, p);
   Polynomial<VeryLongModular> Bp_ = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(B, p);
   Polynomial<VeryLongModular> U_;
   Polynomial<VeryLongModular> V_;
   //cout << "Ap_ = " << Ap_ << ", Bp_ = " << Bp_ << ", gcd(Ap_,Bp_) = " << gcd(Ap_, Bp_) << endl;
   Polynomial<VeryLongModular> G = extended_gcd<Polynomial<VeryLongModular> >(Ap_, Bp_, U_, V_);
   //cout << "G = " << G << ", U_ = " << U_ << ", V_ = " << V_ << endl;
   VeryLongModular g = G.coefficient(0);
   U_ /= g;
   V_ /= g;
   //cout << "U_ = " << U_ << ", V_ = " << V_ << endl;
   Polynomial<VeryLongModular> check = U_ * Ap_ + V_ * Bp_;
   const VeryLongModular one(1L);
   if (check != Polynomial<VeryLongModular>(one))
   {
      std::cout << "Problem U_ Ap_ + V_ Bp_ != 1, check = " << check << std::endl;
      std::cout << "Ap_ = " << Ap_ << ", Bp_ = " << Bp_ << ", gcd(Ap_,Bp_) = " << gcd(Ap_, Bp_) << std::endl;
      std::cout << "U_ = " << U_ << ", V_ = " << V_ << std::endl;
   }
   Polynomial<VeryLong> U = lift<VeryLong, VeryLongModular>(U_);
   Polynomial<VeryLong> V = lift<VeryLong, VeryLongModular>(V_);
   Polynomial<VeryLong> A1;
   Polynomial<VeryLong> B1;

   hensel_lift(p, q, A, B, C, U, V, A1, B1);
   W.push_back(A1);
   VV.pop_front();
   if (VV.size() > 1)
   {
      hensel_lift(p, q, B1, VV, W);
   }
   else
   {
      W.push_back(B1);
   }
}

static long int combinations(int n, int j)
{
   // number of ways of choosing j from n
   if (j > n) return 0;
   // n! / (n-j)!j! = n(n-1)...(n-j+1) / j!
   // we only need to consider n up to 5 for the moment, so we can use
   // a naive calculation
   int res = 1;
   for (int i = 0; i < j; i++)
   {
      res *= n - i;
   }
   for (int i = 0; i < j; i++)
   {
      res /= i + 1;
   }
   return (long int)res;
}

static VeryLong size(const Polynomial<VeryLong>& A)
{
   int n = A.deg();
   VeryLong s(0L);
   for (int i = 0; i <= n; i++)
   {
      s += A.coefficient(i) * A.coefficient(i);
   }
   VeryLong t = s.nth_root(2);
   if (t*t < s) t += 1L;
   return t;
}

static VeryLong bound(const Polynomial<VeryLong>& A, int n)
{
   const VeryLong zero(0L);
   /* find bound on coefficients of any polynomial B that is a factor of A
   // with A having degree m, and B having degree n, n | m
   // From Theorem 3.3.1 the bound is the maximum of the bounds given by
   //            / n - 1 \         / n - 1 \
   // | b_j | <= |       | | A | + |       | | a_m |
   //            \   j   /         \ j - 1 /
   */

   VeryLong sizeA = size(A);
   VeryLong B(zero);
   if (n >= A.deg()) return B;
   VeryLong a_m = A.coefficient(A.deg());
   if (a_m < zero) a_m = -a_m;

   for (int j = 0; j <= n; j++)
   {
      VeryLong B_j = sizeA * VeryLong(combinations(n - 1, j)) + a_m * VeryLong(combinations(n - 1, j - 1));
      if (B_j > B) B = B_j;
   }

   return B;
}

void check(const Polynomial<VeryLong>& U,
           const std::deque<Polynomial<VeryLong> >& Ufactors,
           const VeryLong p_power)
{
   // check that U = prod Ufactors[i] mod p_power
   Polynomial<VeryLong> check(VeryLong(1L));
   for (size_t i = 0; i < Ufactors.size(); i++)
   {
      check *= Ufactors[i];
   }
   Polynomial<VeryLongModular> check1 = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(check, p_power);
   Polynomial<VeryLongModular> check2 = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(U, p_power);
   if (check1 != check2)
   {
      std::cout << "Problem: check1 - check2 = " << check1 - check2 << std::endl;
      std::cout << "check1 = " << check1 << ", check2 = " << check2 << std::endl;
   }
}

void display(const std::vector<std::vector<int> >& v)
{
   std::cout << "Total combinations = " << static_cast<long int>(v.size()) << std::endl;
   for (size_t i = 0; i < v.size(); i++)
   {
      std::cout << "(";
      for (size_t j = 0; j < v[i].size(); j++)
      {
         std::cout << v[i][j];
         if (j < v[i].size() - 1) std::cout << " ";
      }
      std::cout << ")" << std::endl;
   }
}

std::vector<std::vector<int> > generate_combinations(int r, int d)
{
   static int recurse = 0;
   recurse++;
//   cout << "generate_combinations(" << r << ", " << d << ")" << endl;
   std::vector<std::vector<int> > res;
   if (d == 1)
   {
      for (int i = 0; i < r; i++)
      {
         std::vector<int> v;
         v.push_back(i + 1);
         res.push_back(v);
      }
      recurse--;
      return res;
   }

   if (d == r)
   {
      std::vector<int> v;
      for (int i = 0; i < r; i++)
      {
         v.push_back(i + 1);
      }
      res.push_back(v);
      recurse--;
      return res;
   }

   if (recurse == 1 && 2 * d == r)
   {
      std::vector<std::vector<int> > res1 = generate_combinations(r - 1, d - 1);
      for (size_t k = 0; k < res1.size(); k++)
      {
         res1[k].push_back(r);
         res.push_back(res1[k]);
      }
      recurse--;
      return res;
   }

   for (int j = r; j >= d; --j)
   {
      std::vector<std::vector<int> > res1 = generate_combinations(j - 1, d - 1);
      for (size_t k = 0; k < res1.size(); k++)
      {
         res1[k].push_back(j);
         res.push_back(res1[k]);
      }
   }
   recurse--;
   return res;
}

// Algorithm 3.5.7 (Factor in Z[X])
template <> void Polynomial<VeryLong>::factor(const Polynomial<VeryLong>& AA, std::vector<Polynomial<VeryLong> >& factors, VeryLong& cont)
{
   const VeryLong zero(0L);
   const VeryLong two(2L);
   Polynomial<VeryLong> A = AA;
   // Step 1. [Reduce to squarefree and primitive]
   cont = A.content();
   A = A / cont;
   if (A.coefficient(A.deg()) < zero && cont > zero) cont = -cont;
   //cout << "A = " << A << endl;
   Polynomial<VeryLong> A1 = A.derivative();
   //cout << "A' = " << A1 << endl;

   Polynomial<VeryLong> U = A / sub_resultant_GCD(A, A1);
   //cout << "U = " << U << endl;

   if (U.deg() == 1)
   {
      if (U.coefficient(1) < zero) U = -U;
      for (int i = 0; i < A.deg(); i++) factors.push_back(U);
      return;
   }
   VeryLong l_U = U.coefficient(U.deg());

   VeryLong::generate_prime_table();
   VeryLong p = VeryLong::firstPrime();

   // Step 2. [Find a squarefree factorization mod p]
   int done = 0;
   std::vector<Polynomial<VeryLongModular> > Ufactors_;
   while (!done)
   {
      while (l_U % p == zero) p = VeryLong::nextPrime();
      //cout << "p = " << p << endl;
      VeryLongModular::set_default_modulus(p);
      const VeryLongModular one(1L);
      Polynomial<VeryLongModular> Up = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(U, p);
      Polynomial<VeryLongModular> U1p = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(U.derivative(), p);
      //cout << "U mod p = " << Up << endl;
      //cout << "U' mod p = " << U1p << endl;
      Polynomial<VeryLongModular> g = gcd<Polynomial<VeryLongModular> >(Up, U1p);
      //cout << "(U, U') mod p = " << g << endl;

      if (g.deg() == 0)
      {
         factor_over_F_p<VeryLong, VeryLong, VeryLongModular>(U, p, Ufactors_);
         //cout << "factors of " << U << " mod " << p << " are:" << endl;
         Polynomial<VeryLongModular> check(one);
         for (size_t i = 0; i < Ufactors_.size(); i++)
         {
            //cout << Ufactors_[i] << endl;
            check *= Ufactors_[i];
         }
         if (check != Up)
         {
            std::cout << "Problem: check = " << check << ", Up = " << Up << std::endl;
         }
         done = 1;
      }
      else p = VeryLong::nextPrime();
   }

   // Step 3. [Find bound]
   VeryLong B = bound(U, U.deg() / 2);
   //cout << "B = " << B << endl;
   int e = 1;
   VeryLong p_e = p;
   VeryLong B1 = two * B * U.coefficient(U.deg());
   if (B1 < zero) B1 = -B1;

   while (p_e <= B1)
   {
      p_e *= p;
      e++;
   }
   //cout << "e = " << e << endl;
   // e is the smallest integer such that p^e > 2*B*l(U)

   // Step 4. [Lift factorization]
   // Lift factors in Ufactors_ to factors mod p^e
   // Just do it simply first, by repeated calls to hensel_lift
   std::deque<Polynomial<VeryLong> > Ufactors;
   VeryLongModular::set_default_modulus(p);
   for (size_t i = 0; i < Ufactors_.size(); i++)
   {
      Ufactors.push_back(lift<VeryLong, VeryLongModular>(Ufactors_[i]));
   }
   if (Ufactors.size() == 1)
   {
      Ufactors.push_back(Polynomial<VeryLong>(VeryLong(1L)));
   }
   VeryLong p_power = p;
   std::deque<Polynomial<VeryLong> > liftedUfactors;
   for (int i = 0; i < e - 1; i++)
   {
      // at this point
      // i = 0, p_power = p
      // i = 1, p_power = p^2
      // ...
      // i = e - 2, p_power = p^(e-1), so lift to factors mod p^e
      hensel_lift(p, p_power, U, Ufactors, liftedUfactors);
      Ufactors = liftedUfactors;
      liftedUfactors.clear();
      p_power *= p;
      check(U, Ufactors, p_power);
   }
   // leaves loop with p_power = p^e
   //cout << "factors of " << U << " mod " << p << "^" << e << " = " << p_power << ", are:" << endl;

   VeryLongModular::set_default_modulus(p_power);
   Polynomial<VeryLong> check(VeryLong(1L));
   Ufactors_.clear();
   for (size_t i = 0; i < Ufactors.size(); i++)
   {
      //cout << Ufactors[i] << endl;
      Polynomial<VeryLongModular> Ui_ = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(Ufactors[i], p_power);
      //cout << Ui_ << endl;
      Ui_.make_monic();
      //cout << Ui_ << endl;
      Ufactors_.push_back(Ui_);
      Ufactors[i] = lift<VeryLong, VeryLongModular>(Ui_);
      //cout << Ufactors[i] << endl;
      check *= Ufactors[i];
   }
   check *= l_U;
   Polynomial<VeryLongModular> check1 = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(check, p_power);
   Polynomial<VeryLongModular> check2 = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(U, p_power);

   if (check1 - check2 != Polynomial<VeryLongModular>(VeryLongModular(0L)))
   {
      std::cout << "Problem in factors:" << std::endl;
      std::cout << "check1 = " << check1 << std::endl;
      std::cout << "check2 = " << check2 << std::endl;
      std::cout << "check1 - check2 = " << check1 - check2 << std::endl;
   }

   int r = static_cast<int>(Ufactors.size());
   int d = 1;

   Polynomial<VeryLongModular> U_ = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(U, p_power);
   Polynomial<Quotient<VeryLong> > Uq = convert_to_quotient_field<VeryLong>(U);
   const Polynomial<Quotient<VeryLong> > zero_qp(Quotient<VeryLong>(0L));
   done = 0;
   while (!done)
   {
      int done5 = 0;
      while (!done5)
      {
         /* Step 5. [Try combination]
         // find all combinations of factors in Ufactor of length d
         // and check whether they divide U
         // If d = r/2, take i_d = 1
         // If we have r factors and length d, there are
         // / r \
         // |   | = r! / d! (r - d)! combinations
         // \ d /
         // but in the case d = r/2 we only want to consider
         // half of them, since by induction the remainder must be
         // the other factor.
         */
         Polynomial<VeryLong> V;
         std::vector<std::vector<int> > combination_set = generate_combinations(r, d);
         //display(combination_set);
         int factorFound = 0;
         for (size_t comb = 0; !factorFound && comb < combination_set.size(); comb++)
         {
            Polynomial<VeryLongModular> V_(VeryLongModular(1L));
            for (int j = 0; j < d; j++)
            {
               V_ *= Ufactors_[combination_set[comb][j] - 1];
            }
            if (2 * V.deg() <= U.deg())
            {
               V_ *= VeryLongModular(l_U);
            }
            else
            {
               V_ = U_ / V_;
            }
            V = lift<VeryLong, VeryLongModular>(V_);
            for (int j = 0; j <= V.deg(); j++)
            {
               if (two * V.coefficient(j) >= p_power)
               {
                  V.set_coefficient(j, V.coefficient(j) - p_power);
               }
            }
            //cout << "V = " << V << endl;
            Polynomial<VeryLong> Q = (l_U * U) / V;
            Polynomial<VeryLong> UUU = Q;
            UUU *= V;
            //if (Q * V == l_U * U)
            if (UUU == l_U * U)
            {
               //cout << "V divides l(U)U!" << endl;
               //cout << "l(U)U = " << l_U * U << ", V = " << V << endl;
               //cout << "l(U)U / V = " << (l_U * U) / V << endl;
               Polynomial<VeryLong> F = V;
               F.make_primitive();
               if (F.coefficient(F.deg()) < zero) F = -F;
               //cout << "F = " << F << endl;
               U = U / F;
               int v = exponent(A, F);
               for (int k = 0; k < v; k++) factors.push_back(F);
               l_U = U.coefficient(U.deg());
               // remove Ui from Ufactors_
               std::vector<Polynomial<VeryLongModular> > newUf_;
               if (2 * d <= r)
               {
                  for (int k = 0; k < r; k++)
                  {
                     int j = 0;
                     //for (j = 0; combination_set[comb][j] - 1 != k && j < d; j++);
                     for (j = 0; j < d && combination_set[comb][j] - 1 != k; j++);
                     if (j >= d) newUf_.push_back(Ufactors_[k]);
                  }
               }
               else
               {
                  for (int j = 0; j < d; j++)
                  {
                     newUf_.push_back(Ufactors_[combination_set[comb][j] - 1]);
                  }
               }
               Ufactors_ = newUf_;
               r = static_cast<int>(Ufactors_.size());
               if (2 * d > r) done5 = 1;
               factorFound = 1;
            }
         }
         // We are here if we've found a factor among the combinations, or if
         // the combinations have been exhausted with no fact found.
         // In the latter case we want to go to step 6:
         if (!factorFound) done5 = 1;
         // Otherwise we go round again with the reduced U
      } // this loop is exited if done5 = 1

      // Step 6.
      d++;
      if (2 * d > r) done = 1;
   }

   // primitive part of U is remaining factor.
   if (U.deg() > 0)
   {
      Polynomial<VeryLong> F = U;
      F.make_primitive();
      if (F.coefficient(F.deg()) < zero) F = -F;
      //cout << "F = " << F << endl;
      int v = exponent(A, F);
      for (int k = 0; k < v; k++) factors.push_back(F);
   }
}
