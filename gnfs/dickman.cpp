#include <math.h>
#include <vector>

// functions to do with Dickman's rho function
// rho(u) = lim r->inf phi(r, r^(1/u))
// phi(a,b) = |{ x in Z+ : x <= a and P1(x) <= b }|
// P1(x) = largest prime factor of x.

long int P1(long int x)
{
   long int i = 2;
   long int p1 = 1;

   while (i <= sqrt((double)x))
   {
      while (x % i == 0)
      {
         x /= i;
         p1 = i;
      }
      if (i == 2) i++;
      else i += 2;
   }
   if (x > p1) p1 = x;
   return p1;
}

long int phi(long int a, long int b)
{
   int count = 0;
   long int x = 1;
   while (x <= a)
   {
      if (P1(x) <= b) count++;
      x++;
   }
   return count;
}

/*
 *  Compute  rho( x ) where  rho  is the Dickman rho function.
 *  This can be optimized by precomputing the coefficients of
 *  the polynomials and storing those coefficients for future calls.
 *  Code adpated from Scott Contini and from Brian Murphy's thesis.
 */

const int NUMBER_OF_COEFFICIENTS = 55;
class rho_coefficients
{
   public:
      double& operator[] (int i)
      {
         return c[i];
      }
   private:
      double c[NUMBER_OF_COEFFICIENTS];
};

double dickman_rho( double x)
{
   //cout << "rho(" << x << ")" << endl;
   if (x <= 1.0) return 1.0;
   if (x > 50) return 0.0;

   double k = ceil(x);
   double zeta = k - x;

   static std::vector<rho_coefficients*> saved_coefficients;
   static int first_time = 1;

   if (first_time)
   {
      // pre-calculate coefficients
      first_time = 0;
      rho_coefficients* c = new rho_coefficients;
      // calculate the coefficients of rho (2) (u) defined on u in [1,2]
      (*c)[0] = 1.0 - log(2.0);

      double p = 2.0;
      for (int i = 1; i < NUMBER_OF_COEFFICIENTS; i++)
      {
         (*c)[i] = 1.0 /  (i * p);
         p *= 2;
      }
      saved_coefficients.push_back(c);
      double prev_c[NUMBER_OF_COEFFICIENTS];
      double sum = 0.0;

      long int l = 2;
      while (l < 50)
      {
         int i = 0;
         for (i = 0; i < NUMBER_OF_COEFFICIENTS; i++)
         {
            prev_c[i] = (*c)[i];
         }

         c = new rho_coefficients;
         for (i = 1; i < NUMBER_OF_COEFFICIENTS; i++)
         {
            sum = 0.0;

            for (int j = 0; j < i; j++)
            {
               sum += prev_c[j]/(i * pow((l + 1.0),i-j));
            }
            (*c)[i] = sum;
         }

         sum = 0.0;
         for (int j = 1; j < NUMBER_OF_COEFFICIENTS; j++)
         {
            sum += (*c)[j]/(j + 1.0);
         }
         (*c)[0] = sum/l;
         saved_coefficients.push_back(c);
         l++;
      }
   }

   rho_coefficients* c = saved_coefficients[(int)(k-2)];

   // we now have the coefficients for pho (l) (u) defined on u in [l-1,l]
   double rho_k = 0.0;
   double z = 1.0;
   for (int i = 0; i < NUMBER_OF_COEFFICIENTS; i++)
   {
      //rho_k += (*c)[i] * pow(zeta,i);
      rho_k += (*c)[i] * z;
      z *= zeta;
   }

   return rho_k;
}
