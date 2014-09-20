#ifndef _GCD_H
#define _GCD_H
template <class I> I gcd(const I& a, const I& b)
{
   const I zero(0L);
   if (a == zero) return b;
   if (b == zero) return a;
   if (a == b) return a;
   if (a < zero) return gcd<I>(-a,b);
   if (b < zero) return gcd<I>(a,-b);
   //if (a > b) return gcd<I>(b,a);
   if (!(b < a)) return gcd<I>(b,a);
   I prev_r = b;
   I r = a;
   I next_r = b % a;
   while (next_r != zero)
   {
      prev_r = r;
      r = next_r;
      next_r = prev_r % next_r;
   }
   return r;
}

// Algorithm 1.3.6 (Euclid Extended)
template <class I> I extended_gcd(const I& a, const I& b, I& u, I& v)
{
   const I zero(0L);
   const I one(1L);
   // Step 1. [Initialize]
   u = one;
   I d = a;
   if (b == zero)
   {
      v = zero;
      if (d < zero)
      {
         d = -d;
         u = -u;
         v = -v;
      }
      return d;
   }
   I v1 = zero;
   I v3 = b;

   while (1)
   {
      // Step 2. [Finished?]
      if (v3 == zero)
      {
         v = d;
         v -= a * u;
         v /= b;
         if (d < zero)
         {
            d = -d;
            u = -u;
            v = -v;
         }
         return d;
      }

      // Step 3. [Euclidean step]
      I q = d / v3;
      I t3 = d % v3;
      I t1 = u;
      t1 -= q * v1;
      u = v1;
      d = v3;
      v1 = t1;
      v3 = t3;
   }
}

// Algorithm 1.3.6 (Euclid Extended)
// but with |u| and |v| minimal
template <class I> I minimal_extended_gcd(const I& a, const I& b, I& u, I& v)
{
   const I zero(0L);
   const I one(1L);
   // Step 1. [Initialize]
   u = one;
   I d = a;
   if (b == zero)
   {
      v = zero;
      return d;
   }
   I v1 = zero;
   I v3 = b;

   while (1)
   {
      // Step 2. [Finished?]
      if (v3 == zero)
      {
         v = d;
         v -= a * u;
         v /= b;
         // adjust u and v to be minimal
         const I minusone(-1L);
         I lambda = a / d;
         if (lambda < zero) lambda *= minusone;
         I mu = b / d;
         if (mu < zero) mu *= minusone;
         if (v > zero && b > zero)
         {
            I k = v / lambda + one;
            v -= k * lambda;
            u += k * mu;
         }
         else if (v < zero && b < zero)
         {
            I k = v / lambda - one;
            v += k * lambda;
            u -= k * mu;
         }
         else if (v > lambda && b < zero)
         {
            I k = v / lambda;
            v -= k * lambda;
            u += k * mu;
         }
         else if (v < -lambda && b > zero)
         {
            I k = v / lambda;
            v += k * lambda;
            u -= k * mu;
         }
         return d;
      }

      // Step 3. [Euclidean step]
      I q = d / v3;
      I t3 = d % v3;
      I t1 = u;
      t1 -= q * v1;
      u = v1;
      d = v3;
      v1 = t1;
      v3 = t3;
   }
}

// From Algorithm 1.3.6 (Euclid Extended)
template <class I> inline I inverse(I n, I p)
{
   if (n == 1L) return n;
   bool negative = false;
   if (n < 0L)
   {
      n = -n;
      negative = true;
   }

   I u = 1L;
   I q = n;
   I w = 0L;
   I nr = p;
   I nq;
   I nw;
   I r = p;

   if (q < r) q += r;
   while (nr)
   {
      if ((nq = w, nr = q - r) < r)
         {}
      else if ((nq += w, nr -= r) < r)
         {}
      else if ((nq += w, nr -= r) < r)
         {}
      else if ((nq += w, nr -= r) < r)
         {}
      else
      {
         nq = q / r;
         nr = q % r;
         nq *= w;
      }

      nw = u - nq;
      u = w;
      w = nw;

      q = r;
      r = nr;
   }
   if (negative)
   {
      u = -u;
   }
   if (u < 0) return (u + p);
   return u;
}
#endif
