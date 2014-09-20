#ifndef _POW_H
#define _POW_H

// template for generic power function
// Algorithm 1.2.1 (Right-Left Binary)
template <class T, class Z>
T pow(const T& g, const Z& n)
{
   const Z zero(0L);
   const Z one(1L);
   const Z two(2L);
   // Step 1. [Initialize]
   T y(1L);
   if (n == zero) return y;

   if (n < zero)
   {
      return T(0L);
   }

   Z N = n;
   T z = g;

   while (1)
   {
      // Step 2. [Multiply]
      if (N % two == one) y *= z;

      // Step 3. [Halve N]
      N = N / two;
      if (N == zero) return y;
      z *= z;
   }
}

// Adaptation, where n > 0 avoiding need to use identity element
template <class T, class Z>
T pospow(const T& g, const Z& n)
{
   const Z zero(0L);
   const Z one(1L);
   const Z two(2L);
   // Step 1. [Initialize]
   T y = g;
   if (n == one) return y;

   Z N = n / two;
   T z = g * g;
   if (n == two) return z;

   while (1)
   {
      // Step 2. [Multiply]
      if (N % two == one) y *= z;

      // Step 3. [Halve N]
      N = N / two;
      if (N == zero) return y;
      z *= z;
   }
}
#endif
