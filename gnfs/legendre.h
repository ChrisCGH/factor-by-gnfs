#ifndef LEGENDRE_H
#define LEGENDRE_H

// Algorith 1.4.10 (Kronecker)
// adapted to assume b is prime
template <class I> int legendre(const I& aa, const I& bb)
{
   const I zero(0L);
   const I one(1L);
   const I two(2L);
   static int tab2[8] =
   {
      0, 1, 0, -1, 0, -1, 0, 1
   };
   I a(aa);
   I b(bb);

   // 1. [Test b equal to 0]
   // ignore this step since b is assumed to be prime

   // 2. [Remove 2's from b]
   int v = 0;
   int k = 1;

   while (true)
   {
      // 3. [Finished?]
      if (a == zero)
      {
         if (b > one)
         {
            return 0;
         }
         else
         {
            return k;
         }
      }

      while (a % two == zero)
      {
         ++v;
         a /= two;
      }

      if (v & 1)
      {
         int temp = b % 8L;
         if (tab2[temp] < 0)
         {
            k = -k;
         }
      }

      // 4. [Apply reciprocity]

      long int aa = a % 4L;
      long int bb = b % 4L;
      if (aa & bb & 2)
      {
         k = -k;
      }

      I r = a;
      a = b % r;
      b = r;
   }

}
#endif
