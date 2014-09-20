#include <vector>
#ifdef MS
#include <xutility>
#endif
#include "Matrix.h"
#include "VeryLong.h"

// lll.cpp - implementation of integral LLL lattice basis reduction
// based on description in Cohen (Algorithm 2.6.7, p94).

void myswap(std::vector<VeryLong>& v1, std::vector<VeryLong>& v2)
{
   std::vector<VeryLong> tmp = v1;
   v1 = v2;
   v2 = tmp;
}

static VeryLong dot(const std::vector<VeryLong>& a, const std::vector<VeryLong>& b)
{
   VeryLong result = 0L;
   for (size_t i = 0; i < a.size(); i++)
   {
      result += a[i] * b[i];
   }
   return result;
}

void SWAPI(int kk, int k_max,
           std::vector<std::vector<VeryLong > >& basis,
           std::vector<std::vector<VeryLong > >& H,
           std::vector<std::vector<VeryLong > >& lambda,
           std::vector<VeryLong>& d)
{
   //cout << "SWAPI(" << kk << "," << k_max << ")" << endl;
   myswap(H[kk-1], H[kk-2]);
   myswap(basis[kk-1], basis[kk-2]);
   if (kk > 2)
   {
      for (int j = 1; j <= kk - 2; j++)
      {
         VeryLong tmp = lambda[kk-1][j-1];
         lambda[kk-1][j-1] = lambda[kk-2][j-1];
         lambda[kk-2][j-1] = tmp;
      }
   }
   VeryLong lambda1 = lambda[kk-1][kk-2];
   VeryLong B = (d[kk-2] * d[kk] + lambda1 * lambda1) / d[kk-1];
   for (int i = kk + 1; i <= k_max; i++)
   {
      VeryLong t = lambda[i-1][kk-1];
      lambda[i-1][kk-1] = (d[kk] * lambda[i-1][kk-2] - lambda1 * t) / d[kk-1];
      lambda[i-1][kk-2] = (B * t + lambda1 * lambda[i-1][kk-1]) / d[kk];
   }
   d[kk-1] = B;
   //cout << "d[" << kk - 1 << "] <- B = " << B << endl;
}

void REDI(int kk, int ll,
          std::vector<std::vector<VeryLong > >& basis,
          std::vector<std::vector<VeryLong > >& H,
          std::vector<std::vector<VeryLong > >& lambda,
          std::vector<VeryLong>& d)
{
   //cout << "REDI(" << kk << "," << ll << ")" << endl;
   int n = basis.size();
   VeryLong a = lambda[kk-1][ll-1];
   VeryLong b = d[ll];
   VeryLong tmp = a * 2L;
   if (tmp < 0L) tmp = -tmp;
   VeryLong tmp1 = b;
   if (tmp1 < 0L) tmp1 = -tmp1;
   if (tmp > b)
   {
      VeryLong q = (2L * a + b) / (2L * b);

      //cout << "q = " << q << endl;
      //cout << "lambda[" << kk << "][" << ll << "] = " << a << endl;
      //cout << "d[" << ll << "] = " << b << endl;

      for (int jj = 0; jj < n; jj++)
      {
         H[kk-1][jj] -= q * H[ll-1][jj];
      }
      for (size_t jjj = 0; jjj < basis[0].size(); jjj++)
      {
         basis[kk-1][jjj] -= q * basis[ll-1][jjj];
      }
      lambda[kk-1][ll-1] -= q * d[ll];
      for (int ii = 1; ii <= ll - 1; ii++)
      {
         lambda[kk-1][ii-1] -= q * lambda[ll-1][ii-1];
      }
   }
}

int LLL_reduce(std::vector<std::vector<VeryLong > >& basis,
               std::vector<std::vector<VeryLong > >& H)
{
   // step 1. [Initialize]
   //cout << "LLL step 1. [Initialize]" << endl;
   int k = 2;
   int k_max = 1;
   int n = basis.size();
   std::vector<VeryLong> d;
   d.resize(n+1);
   d[0] = 1L;
   d[1] = dot(basis[0], basis[0]);

   std::vector<std::vector<VeryLong > > lambda;
   lambda.resize(n);
   int ii = 0;
   for (ii = 0; ii < n; ii++) lambda[ii].resize(n);

   // H = I_n
   for (ii = 0; ii < n; ii++)
   {
      for (int jj = 0; jj < n; jj++)
      {
         if (ii == jj) H[ii][jj] = 1L;
         else H[ii][jj] = 0L;
      }
   }

   while (k <= n)
   {
      // step 2. [Incremental Gram-Schmidt]
      //cout << "step 2. [Incremental Gram-Schmidt] : " << k << endl;
      if (k > k_max)
      {
         k_max = k;
         for (int j = 1; j <= k; j++)
         {
            VeryLong u = dot(basis[k-1], basis[j-1]);
            for (int i = 1; i <= j - 1; i++)
            {
               u = (d[i] * u - lambda[k-1][i-1] * lambda[j-1][i-1])/d[i-1];
            }
            if (j < k)
            {
               lambda[k-1][j-1] = u;
            }
            else if (j == k)
            {
               d[k] = u;
               //cout << "d[" << k << "] <- u = " << u << endl;
               if (d[k] == 0L)
               {
                  std::cout << "LLL_reduce: ERROR: basis is not a basis" << std::endl;
                  return -1;
               }
            }
         }
      }
      // step 3. [Test LLL condition]
      //cout << "step 3. [Test LLL condition]" << endl;
      // sub-algorithm REDI in line
      int step_3_done = 0;
      while (!step_3_done)
      {
         REDI(k, k-1, basis, H, lambda, d);
         if (4L * d[k] * d[k-2] < 3L * d[k-1] * d[k-1] - 4L * lambda[k-1][k-2] * lambda[k-1][k-2])
         {
            SWAPI(k, k_max, basis, H, lambda, d);
            k = std::max(2, k - 1);
         }
         else
         {
            for (int l = k - 2; l >= 1; l--)
            {
               REDI(k, l, basis, H, lambda, d);
            }
            k++;
            step_3_done = 1;
         }
      }
   }
   // step 4. [Finished]
   //cout << "step 4. [Finished]" << endl;
   return 0;
}

void transform_row(std::vector<VeryLong>& c1, std::vector<VeryLong>& c2,
                   const VeryLong& x, const VeryLong& y,
                   const VeryLong& u, const VeryLong& v)
{
   int m = c1.size();
   VeryLong t1;
   VeryLong t2;
   VeryLong t3;
   VeryLong t4;
   for (int i = 0; i < m; i++)
   {
      t1 = x * c1[i];
      t2 = y * c2[i];
      t1 += t2;
      t3 = u * c1[i];
      t4 = v * c2[i];
      t3 += t4;
      c1[i] = t1;
      c2[i] = t3;
   }
}

void transform_row(VeryLong& c1, VeryLong& c2,
                   const VeryLong& x, const VeryLong& y,
                   const VeryLong& u, const VeryLong& v)
{
   VeryLong t1;
   VeryLong t2;
   VeryLong t3;
   VeryLong t4;

   t1 = x * c1;
   t2 = y * c2;
   t1 += t2;
   t3 = u * c1;
   t4 = v * c2;
   t3 += t4;
   c1 = t1;
   c2 = t3;

}


void SWAPI_1(int kk, int k_max,
             Matrix<VeryLong>& b,
             std::vector<int>& p,
             Matrix<VeryLong>& H,
             Matrix<VeryLong>& lambda,
             std::vector<VeryLong>& d)
{
//   cout << "SWAPI_1(" << kk << "," << k_max << ")" << endl;
   //int n = b.size();
   int i = 0;
   int j = 0;
   if (p[kk-1] != 0)
   {
      //myswap(H[kk-1], H[kk-2]);
      H.swap(kk-1, kk-2);
      //myswap(b[kk-1], b[kk-2]);
      b.swap(kk-1, kk-2);

      for (j = 1; j <= kk - 2; j++)
      {
         if (p[j-1] != 0)
         {
            VeryLong tmp = lambda(kk-1, p[j-1]-1);
            lambda(kk-1, p[j-1]-1) = lambda(kk-2, p[j-1]-1);
            lambda(kk-2, p[j-1]-1) = tmp;
            //cout << "lambda[" << kk << "][" << p[j-1] << "] = " << lambda[kk-1][p[j-1]-1] << endl;
            //cout << "lambda[" << kk-1 << "][" << p[j-1] << "] = " << lambda[kk-2][p[j-1]-1] << endl;
         }
      }

      for (i = kk + 1; i <= k_max; i++)
      {
         VeryLong t1 = (lambda(i-1, p[kk-1]-2) * lambda(kk-1, p[kk-1]-2) +
                        lambda(i-1, p[kk-1]-1) * d[p[kk-1]-2])/d[p[kk-1]-1];
         VeryLong t2 = (lambda(i-1, p[kk-1]-2) * d[p[kk-1]] -
                        lambda(i-1, p[kk-1]-1) * lambda(kk-1, p[kk-1]-2))/d[p[kk-1]-1];
         lambda(i-1, p[kk-1]-2) = t1;
         lambda(i-1, p[kk-1]-1) = t2;
         //cout << "lambda[" << i << "][" << p[kk-1] << "] = " << t2 << endl;
         //cout << "lambda[" << i << "][" << p[kk-1]-1 << "] = " << t1 << endl;
      }
      d[p[kk-1]-1] = (d[p[kk-1]] * d[p[kk-1]-2] +
                      lambda(kk-1, p[kk-1]-2) * lambda(kk-1, p[kk-1]-2))/d[p[kk-1]-1];
      //cout << "d[p[kk-1]-1 = " << p[kk-1]-1 << "] <- " << d[p[kk-1]-1] << endl;

   }
   else if (lambda(kk-1, p[kk-2]-1) != 0L)
   {
      // swap case 2:
      VeryLong x;
      VeryLong y;
      VeryLong e = extended_gcd(lambda(kk-1, p[kk-2]-1), d[p[kk-2]], x, y);
      //cout << "check xgcd : " << lambda[kk-1][p[kk-2]-1] % e << endl;
      //cout << "           : " << d[p[kk-2]] % e << endl;

      VeryLong t1 = lambda(kk-1, p[kk-2]-1) / e;
      VeryLong t2 = d[p[kk-2]] / e;
      //cout << "e = " << e << endl;
      //cout << "x = " << x << endl;
      //cout << "y = " << y << endl;
      //cout << "t1 = " << t1 << endl;
      //cout << "t2 = " << t2 << endl;

      VeryLong t3 = t2;
      t2 = t2 * -1L;
      //transform_row(b[kk-2], b[kk-1], t1, t2, y, x);
      b.transform_row(kk-2, kk-1, t1, t2, y, x);
      //transform_row(H[kk-2], H[kk-1], t1, t2, y ,x);
      H.transform_row(kk-2, kk-1, t1, t2, y ,x);
      for (int j = 1; j <= kk-2; j++)
      {
         if (p[j-1] != 0)
         {
            //cout << "lambda[" << kk-1 << "][" << p[j-1] << "] = " << lambda[kk-2][p[j-1]-1] << endl;
            //cout << "lambda[" << kk << "][" << p[j-1] << "] = " << lambda[kk-1][p[j-1]-1] << endl;
            transform_row(lambda(kk-2, p[j-1]-1), lambda(kk-1, p[j-1]-1), t1, t2, y, x);
            //cout << "lambda[" << kk-1 << "][" << p[j-1] << "] = " << lambda[kk-2][p[j-1]-1] << endl;
            //cout << "lambda[" << kk << "][" << p[j-1] << "] = " << lambda[kk-1][p[j-1]-1] << endl;
         }
      }
      t2 = t2 * t2;
      d[p[kk-2]] /= t2;
      //cout << "d[p[kk-2] = " << p[kk-2] << "] <- " << d[p[kk-2]] << endl;

      for (i = kk+1; i <= k_max; i++)
      {
         if (p[i-1] != 0)
         {
            d[p[i-1]] /= t2;
            //cout << "d[p[i-1] = " << p[i-1] << "] <- " << d[p[i-1]] << endl;
            for (int j = i+1; j <= k_max; j++)
            {
               lambda(j-1, p[i-1]-1) /= t2;
               //cout << "lambda[" << j << "][" << p[i-1] << "] = " << lambda[j-1][p[i-1]-1] << endl;
            }
         }
      }

      for (i = kk+1; i <= k_max; i++)
      {
         lambda(i-1, p[kk-2]-1) /= t3;
         //cout << "lambda[" << i << "][" << p[kk-2] << "] = " << lambda[i-1][p[kk-2]-1] << endl;
      }

      int tmp = p[kk-2];
      p[kk-2] = p[kk-1];
      p[kk-1] = tmp;
   }
   else
   {
      // swap case 3:
      //myswap(H[kk-1], H[kk-2]);
      H.swap(kk-1, kk-2);
      //myswap(b[kk-1], b[kk-2]);
      b.swap(kk-1, kk-2);

      for (int j = 1; j <= kk - 2; j++)
      {
         if (p[j-1] != 0)
         {
            VeryLong tmp = lambda(kk-1, p[j-1]-1);
            lambda(kk-1, p[j-1]-1) = lambda(kk-2, p[j-1]-1);
            lambda(kk-2, p[j-1]-1) = tmp;
            //cout << "lambda[" << kk << "][" << p[j-1] << "] = " << lambda[kk-1][p[j-1]-1] << endl;
            //cout << "lambda[" << kk-1 << "][" << p[j-1] << "] = " << lambda[kk-2][p[j-1]-1] << endl;

         }
      }

      int tmp = p[kk-2];
      p[kk-2] = p[kk-1];
      p[kk-1] = tmp;
   }

//   cout << "H : " << endl << H.transpose();
}

void REDI_1(int kk, int ll,
            Matrix<VeryLong>& b,
            std::vector<int>& p,
            Matrix<VeryLong>& H,
            Matrix<VeryLong>& lambda,
            std::vector<VeryLong>& d)
{
//   cout << "REDI_1(" << kk << "," << ll << ")" << endl;
   if (p[ll-1] == 0) return;

   int n = b.rows();
   VeryLong a = lambda(kk-1, p[ll-1]-1);
   VeryLong c = d[p[ll-1]];
   VeryLong tmp = a * 2L;
   if (tmp < 0L) tmp = -tmp;
   if (tmp > c)
   {
      VeryLong q = a / c;
      VeryLong r = a % c;
      r = r + r;
      if (r > c || (r == c && q < 0L)) q += 1L;


//	  cout << "q = " << q << endl;
      //cout << "lambda[" << kk << "][" << p[ll-1] << "] = " << a << endl;
      //cout << "d[" << p[ll-1] << "] = " << c << endl;

      for (int jj = 0; jj < n; jj++)
      {
         H(kk-1, jj) -= q * H(ll-1, jj);
      }
      //for (size_t jjj = 0; jjj < b[0].size(); jjj++)
      for (size_t jjj = 0; jjj < b.columns(); jjj++)
      {
         b(kk-1, jjj) -= q * b(ll-1, jjj);
      }
      lambda(kk-1, p[ll-1]-1) -= q * d[p[ll-1]];
      //cout << "lambda[" << kk << "][" << p[ll-1] << "] = " << lambda[kk-1][p[ll-1]-1] << endl;
      for (int ii = 1; ii <= ll - 1; ii++)
      {
         if (p[ii-1] != 0)
         {
            lambda(kk-1, p[ii-1]-1) -= q * lambda(ll-1, p[ii-1]-1);
            //cout << "lambda[" << kk << "][" << p[ii-1] << "] = " << lambda[kk-1][p[ii-1]-1] << endl;
         }
      }
   }
//   cout << "H : " << endl << H.transpose();
}

int LLL_reduce_1(Matrix<VeryLong>& b,
                 Matrix<VeryLong>& H)
{
   // NOTE: this does LLL reduction on rows of b, not
   // columns as in Cohen
   // generalisation to a set of not necessarily independent
   // vectors b
   // acknowledgements to NTL
   // step 1. [Initialize]
   //cout << "LLL step 1. [Initialize]" << endl;
   int k = 1;
   int k_max = 0;
   int n = b.rows();

   // corresponds to P in NTL LLL.cpp
   std::vector<int> p;
   p.resize(n);
   int s = 0;

   std::vector<VeryLong> d;
   d.resize(n+1);
   d[0] = 1L;
   //cout << "d[0] <- " << d[0] << endl;

   Matrix<VeryLong> lambda(n, VeryLong(0L), 0);

   H = Matrix<VeryLong>(n, VeryLong(1L), 0);

   while (k <= n)
   {
      // step 2. [Incremental Gram-Schmidt]
      //cout << "step 2. [Incremental Gram-Schmidt] : " << k << endl;
      if (k > k_max)
      {
         k_max = k;
         VeryLong u;
         for (int j = 1; j <= k-1; j++)
         {
            int posj = p[j-1];
            if (posj == 0) continue;
            //u = dot(b[k-1], b[j-1]);
            u = b.dot(k-1, j-1);
            for (int i = 1; i <= posj-1; i++)
            {
               u = (d[i] * u - lambda(k-1, i-1) * lambda(j-1, i-1))/d[i-1];
            }
            lambda(k-1, posj-1) = u;
            //cout << "lambda[" << k << "][" << posj << "] = " << u << endl;
         }

         //u = dot(b[k-1], b[k-1]);
         u = b.dot(k-1, k-1);

         //cout << "u = " << u << endl;

         for (int i = 1; i <= s; i++)
         {
            u = (d[i] * u - lambda(k-1, i-1) * lambda(k-1, i-1))/d[i-1];
            if (u < 0L)
            {
               //cout << "u = " << u << endl;
            }
         }
         if (u == 0L)
         {
            p[k-1] = 0;
         }
         else
         {
            s++;
            p[k-1] = s;
            //cout << "d[s = " << s << "] <- " << u << endl;
            d[s] = u;
         }
      }
      if (k == 1)
      {
         k++;
//         cout << "(1) k changed up : " << k << endl;
         continue;
      }

      // step 3. [Test LLL condition]
      //cout << "step 3. [Test LLL condition]" << endl;
      // sub-algorithm REDI in line

      REDI_1(k, k-1, b, p, H, lambda, d);
      if (p[k-2] != 0 &&
            (p[k-1] == 0 ||
             (4L * d[p[k-1]] * d[p[k-1]-2] < 3L * d[p[k-1]-1] * d[p[k-1]-1] - 4L * lambda(k-1, p[k-1]-2) * lambda(k-1, p[k-1]-2))))
      {
         SWAPI_1(k, k_max, b, p, H, lambda, d);
         k--;
//         cout << "(1) k changed down : " << k << endl;
      }
      else
      {
         for (int l = k - 2; l >= 1; l--)
         {
            REDI_1(k, l, b, p, H, lambda, d);
         }
         k++;
//         cout << "(1) k changed up : " << k << endl;
      }
   }
   // step 4. [Finished]
   //cout << "step 4. [Finished] : s = " << s << endl;
   return s;
}

int LLL_reduce_1_on_columns(Matrix<VeryLong>& b,
                            Matrix<VeryLong>& H)
{
   Matrix<VeryLong> bt = b.transpose();
   Matrix<VeryLong> HT = H.transpose();
   int s = LLL_reduce_1(bt, HT);
   b = bt.transpose();
   H = HT.transpose();
   return s;
}

void SWAPI_2(int kk, int k_max,
             Matrix<VeryLong>& b,
             std::vector<int>& p,
             Matrix<VeryLong>& lambda,
             std::vector<VeryLong>& d)
{
//   cout << "SWAPI_1(" << kk << "," << k_max << ")" << endl;
   //int n = b.size();
   int i = 0;
   int j = 0;
   if (p[kk-1] != 0)
   {
      //myswap(b[kk-1], b[kk-2]);
      b.swap(kk-1, kk-2);

      for (j = 1; j <= kk - 2; j++)
      {
         if (p[j-1] != 0)
         {
            VeryLong tmp = lambda(kk-1, p[j-1]-1);
            lambda(kk-1, p[j-1]-1) = lambda(kk-2, p[j-1]-1);
            lambda(kk-2, p[j-1]-1) = tmp;
            //cout << "lambda[" << kk << "][" << p[j-1] << "] = " << lambda[kk-1][p[j-1]-1] << endl;
            //cout << "lambda[" << kk-1 << "][" << p[j-1] << "] = " << lambda[kk-2][p[j-1]-1] << endl;
         }
      }

      for (i = kk + 1; i <= k_max; i++)
      {
         VeryLong t1 = (lambda(i-1, p[kk-1]-2) * lambda(kk-1, p[kk-1]-2) +
                        lambda(i-1, p[kk-1]-1) * d[p[kk-1]-2])/d[p[kk-1]-1];
         VeryLong t2 = (lambda(i-1, p[kk-1]-2) * d[p[kk-1]] -
                        lambda(i-1, p[kk-1]-1) * lambda(kk-1, p[kk-1]-2))/d[p[kk-1]-1];
         lambda(i-1, p[kk-1]-2) = t1;
         lambda(i-1, p[kk-1]-1) = t2;
         //cout << "lambda[" << i << "][" << p[kk-1] << "] = " << t2 << endl;
         //cout << "lambda[" << i << "][" << p[kk-1]-1 << "] = " << t1 << endl;
      }
      d[p[kk-1]-1] = (d[p[kk-1]] * d[p[kk-1]-2] +
                      lambda(kk-1, p[kk-1]-2) * lambda(kk-1, p[kk-1]-2))/d[p[kk-1]-1];
      //cout << "d[p[kk-1]-1 = " << p[kk-1]-1 << "] <- " << d[p[kk-1]-1] << endl;

   }
   else if (lambda(kk-1, p[kk-2]-1) != 0L)
   {
      // swap case 2:
      VeryLong x;
      VeryLong y;
      VeryLong e = extended_gcd(lambda(kk-1, p[kk-2]-1), d[p[kk-2]], x, y);
      //cout << "check xgcd : " << lambda[kk-1][p[kk-2]-1] % e << endl;
      //cout << "           : " << d[p[kk-2]] % e << endl;

      VeryLong t1 = lambda(kk-1, p[kk-2]-1) / e;
      VeryLong t2 = d[p[kk-2]] / e;
      //cout << "e = " << e << endl;
      //cout << "x = " << x << endl;
      //cout << "y = " << y << endl;
      //cout << "t1 = " << t1 << endl;
      //cout << "t2 = " << t2 << endl;

      VeryLong t3 = t2;
      t2 = t2 * -1L;
      //transform_row(b[kk-2], b[kk-1], t1, t2, y, x);
      b.transform_row(kk-2, kk-1, t1, t2, y, x);
      for (int j = 1; j <= kk-2; j++)
      {
         if (p[j-1] != 0)
         {
            //cout << "lambda[" << kk-1 << "][" << p[j-1] << "] = " << lambda[kk-2][p[j-1]-1] << endl;
            //cout << "lambda[" << kk << "][" << p[j-1] << "] = " << lambda[kk-1][p[j-1]-1] << endl;
            transform_row(lambda(kk-2, p[j-1]-1), lambda(kk-1, p[j-1]-1), t1, t2, y, x);
            //cout << "lambda[" << kk-1 << "][" << p[j-1] << "] = " << lambda[kk-2][p[j-1]-1] << endl;
            //cout << "lambda[" << kk << "][" << p[j-1] << "] = " << lambda[kk-1][p[j-1]-1] << endl;
         }
      }
      t2 = t2 * t2;
      d[p[kk-2]] /= t2;
      //cout << "d[p[kk-2] = " << p[kk-2] << "] <- " << d[p[kk-2]] << endl;

      for (i = kk+1; i <= k_max; i++)
      {
         if (p[i-1] != 0)
         {
            d[p[i-1]] /= t2;
            //cout << "d[p[i-1] = " << p[i-1] << "] <- " << d[p[i-1]] << endl;
            for (int j = i+1; j <= k_max; j++)
            {
               lambda(j-1, p[i-1]-1) /= t2;
               //cout << "lambda[" << j << "][" << p[i-1] << "] = " << lambda[j-1][p[i-1]-1] << endl;
            }
         }
      }

      for (i = kk+1; i <= k_max; i++)
      {
         lambda(i-1, p[kk-2]-1) /= t3;
         //cout << "lambda[" << i << "][" << p[kk-2] << "] = " << lambda[i-1][p[kk-2]-1] << endl;
      }

      int tmp = p[kk-2];
      p[kk-2] = p[kk-1];
      p[kk-1] = tmp;
   }
   else
   {
      // swap case 3:
      //myswap(b[kk-1], b[kk-2]);
      b.swap(kk-1, kk-2);

      for (int j = 1; j <= kk - 2; j++)
      {
         if (p[j-1] != 0)
         {
            VeryLong tmp = lambda(kk-1, p[j-1]-1);
            lambda(kk-1, p[j-1]-1) = lambda(kk-2, p[j-1]-1);
            lambda(kk-2, p[j-1]-1) = tmp;
            //cout << "lambda[" << kk << "][" << p[j-1] << "] = " << lambda[kk-1][p[j-1]-1] << endl;
            //cout << "lambda[" << kk-1 << "][" << p[j-1] << "] = " << lambda[kk-2][p[j-1]-1] << endl;

         }
      }

      int tmp = p[kk-2];
      p[kk-2] = p[kk-1];
      p[kk-1] = tmp;
   }

}

void SWAPI_2_on_columns(int kk, int k_max,
             Matrix<VeryLong>& b,
             std::vector<int>& p,
             Matrix<VeryLong>& lambda,
             std::vector<VeryLong>& d)
{
//   cout << "SWAPI_1(" << kk << "," << k_max << ")" << endl;
   //int n = b.size();
   int i = 0;
   int j = 0;
   if (p[kk-1] != 0)
   {
      //myswap(b[kk-1], b[kk-2]);
      b.swap_columns(kk-1, kk-2);

      for (j = 1; j <= kk - 2; j++)
      {
         if (p[j-1] != 0)
         {
            VeryLong tmp = lambda(kk-1, p[j-1]-1);
            lambda(kk-1, p[j-1]-1) = lambda(kk-2, p[j-1]-1);
            lambda(kk-2, p[j-1]-1) = tmp;
            //cout << "lambda[" << kk << "][" << p[j-1] << "] = " << lambda[kk-1][p[j-1]-1] << endl;
            //cout << "lambda[" << kk-1 << "][" << p[j-1] << "] = " << lambda[kk-2][p[j-1]-1] << endl;
         }
      }

      for (i = kk + 1; i <= k_max; i++)
      {
         VeryLong t1 = (lambda(i-1, p[kk-1]-2) * lambda(kk-1, p[kk-1]-2) +
                        lambda(i-1, p[kk-1]-1) * d[p[kk-1]-2])/d[p[kk-1]-1];
         VeryLong t2 = (lambda(i-1, p[kk-1]-2) * d[p[kk-1]] -
                        lambda(i-1, p[kk-1]-1) * lambda(kk-1, p[kk-1]-2))/d[p[kk-1]-1];
         lambda(i-1, p[kk-1]-2) = t1;
         lambda(i-1, p[kk-1]-1) = t2;
         //cout << "lambda[" << i << "][" << p[kk-1] << "] = " << t2 << endl;
         //cout << "lambda[" << i << "][" << p[kk-1]-1 << "] = " << t1 << endl;
      }
      d[p[kk-1]-1] = (d[p[kk-1]] * d[p[kk-1]-2] +
                      lambda(kk-1, p[kk-1]-2) * lambda(kk-1, p[kk-1]-2))/d[p[kk-1]-1];
      //cout << "d[p[kk-1]-1 = " << p[kk-1]-1 << "] <- " << d[p[kk-1]-1] << endl;

   }
   else if (lambda(kk-1, p[kk-2]-1) != 0L)
   {
      // swap case 2:
      VeryLong x;
      VeryLong y;
      VeryLong e = extended_gcd(lambda(kk-1, p[kk-2]-1), d[p[kk-2]], x, y);
      //cout << "check xgcd : " << lambda[kk-1][p[kk-2]-1] % e << endl;
      //cout << "           : " << d[p[kk-2]] % e << endl;

      VeryLong t1 = lambda(kk-1, p[kk-2]-1) / e;
      VeryLong t2 = d[p[kk-2]] / e;
      //cout << "e = " << e << endl;
      //cout << "x = " << x << endl;
      //cout << "y = " << y << endl;
      //cout << "t1 = " << t1 << endl;
      //cout << "t2 = " << t2 << endl;

      VeryLong t3 = t2;
      t2 = t2 * -1L;
      //transform_row(b[kk-2], b[kk-1], t1, t2, y, x);
      b.transform_column(kk-2, kk-1, t1, t2, y, x);
      for (int j = 1; j <= kk-2; j++)
      {
         if (p[j-1] != 0)
         {
            //cout << "lambda[" << kk-1 << "][" << p[j-1] << "] = " << lambda[kk-2][p[j-1]-1] << endl;
            //cout << "lambda[" << kk << "][" << p[j-1] << "] = " << lambda[kk-1][p[j-1]-1] << endl;
            transform_row(lambda(kk-2, p[j-1]-1), lambda(kk-1, p[j-1]-1), t1, t2, y, x);
            //cout << "lambda[" << kk-1 << "][" << p[j-1] << "] = " << lambda[kk-2][p[j-1]-1] << endl;
            //cout << "lambda[" << kk << "][" << p[j-1] << "] = " << lambda[kk-1][p[j-1]-1] << endl;
         }
      }
      t2 = t2 * t2;
      d[p[kk-2]] /= t2;
      //cout << "d[p[kk-2] = " << p[kk-2] << "] <- " << d[p[kk-2]] << endl;

      for (i = kk+1; i <= k_max; i++)
      {
         if (p[i-1] != 0)
         {
            d[p[i-1]] /= t2;
            //cout << "d[p[i-1] = " << p[i-1] << "] <- " << d[p[i-1]] << endl;
            for (int j = i+1; j <= k_max; j++)
            {
               lambda(j-1, p[i-1]-1) /= t2;
               //cout << "lambda[" << j << "][" << p[i-1] << "] = " << lambda[j-1][p[i-1]-1] << endl;
            }
         }
      }

      for (i = kk+1; i <= k_max; i++)
      {
         lambda(i-1, p[kk-2]-1) /= t3;
         //cout << "lambda[" << i << "][" << p[kk-2] << "] = " << lambda[i-1][p[kk-2]-1] << endl;
      }

      int tmp = p[kk-2];
      p[kk-2] = p[kk-1];
      p[kk-1] = tmp;
   }
   else
   {
      // swap case 3:
      //myswap(b[kk-1], b[kk-2]);
      b.swap_columns(kk-1, kk-2);

      for (int j = 1; j <= kk - 2; j++)
      {
         if (p[j-1] != 0)
         {
            VeryLong tmp = lambda(kk-1, p[j-1]-1);
            lambda(kk-1, p[j-1]-1) = lambda(kk-2, p[j-1]-1);
            lambda(kk-2, p[j-1]-1) = tmp;
            //cout << "lambda[" << kk << "][" << p[j-1] << "] = " << lambda[kk-1][p[j-1]-1] << endl;
            //cout << "lambda[" << kk-1 << "][" << p[j-1] << "] = " << lambda[kk-2][p[j-1]-1] << endl;

         }
      }

      int tmp = p[kk-2];
      p[kk-2] = p[kk-1];
      p[kk-1] = tmp;
   }

}

void REDI_2(int kk, int ll,
            Matrix<VeryLong>& b,
            std::vector<int>& p,
            Matrix<VeryLong>& lambda,
            std::vector<VeryLong>& d)
{
//   cout << "REDI_1(" << kk << "," << ll << ")" << endl;
   if (p[ll-1] == 0) return;

   //int n = b.rows();
   VeryLong a = lambda(kk-1, p[ll-1]-1);
   VeryLong c = d[p[ll-1]];
   VeryLong tmp = a * 2L;
   if (tmp < 0L) tmp = -tmp;
   if (tmp > c)
   {
      VeryLong q = a / c;
      VeryLong r = a % c;
      r = r + r;
      if (r > c || (r == c && q < 0L)) q += 1L;


//	  cout << "q = " << q << endl;
      //cout << "lambda[" << kk << "][" << p[ll-1] << "] = " << a << endl;
      //cout << "d[" << p[ll-1] << "] = " << c << endl;

      //for (size_t jjj = 0; jjj < b[0].size(); jjj++)
      for (size_t jjj = 0; jjj < b.columns(); jjj++)
      {
         b(kk-1, jjj) -= q * b(ll-1, jjj);
      }
      lambda(kk-1, p[ll-1]-1) -= q * d[p[ll-1]];
      //cout << "lambda[" << kk << "][" << p[ll-1] << "] = " << lambda[kk-1][p[ll-1]-1] << endl;
      for (int ii = 1; ii <= ll - 1; ii++)
      {
         if (p[ii-1] != 0)
         {
            lambda(kk-1, p[ii-1]-1) -= q * lambda(ll-1, p[ii-1]-1);
            //cout << "lambda[" << kk << "][" << p[ii-1] << "] = " << lambda[kk-1][p[ii-1]-1] << endl;
         }
      }
   }
}

void REDI_2_on_columns(int kk, int ll,
            Matrix<VeryLong>& b,
            const std::vector<int>& p,
            Matrix<VeryLong>& lambda,
            const VeryLong& c)
{
   if (p[ll-1] == 0) return;

   //int n = b.columns();
   VeryLong a = lambda(kk-1, p[ll-1]-1);
   VeryLong tmp = a * 2L;
   if (tmp < 0L) tmp = -tmp;
   if (tmp > c)
   {
      VeryLong q = a / c;
      VeryLong r = a % c;
      r = r + r;
      if (r > c || (r == c && q < 0L)) q += 1L;

      for (size_t jjj = 0; jjj < b.rows(); jjj++)
      {
         b(jjj, kk-1) -= q * b(jjj, ll-1);
      }
      lambda(kk-1, p[ll-1]-1) -= q * c;
      for (int ii = 0; ii < ll - 1; ii++)
      {
         if (p[ii] != 0)
         {
            lambda(kk-1, p[ii]-1) -= q * lambda(ll-1, p[ii]-1);
         }
      }
   }
}

int LLL_reduce_2(Matrix<VeryLong>& b)
{
   // NOTE: this does LLL reduction on rows of b, not
   // columns as in Cohen
   // generalisation to a set of not necessarily independent
   // vectors b
   // acknowledgements to NTL
   // step 1. [Initialize]
   //cout << "LLL step 1. [Initialize]" << endl;
   int k = 1;
   int k_max = 0;
   int n = b.rows();

   // corresponds to P in NTL LLL.cpp
   std::vector<int> p;
   p.resize(n);
   int s = 0;

   std::vector<VeryLong> d;
   d.resize(n+1);
   d[0] = 1L;
   //cout << "d[0] <- " << d[0] << endl;

   Matrix<VeryLong> lambda(n, VeryLong(0L), 0);

   while (k <= n)
   {
      // step 2. [Incremental Gram-Schmidt]
      //cout << "step 2. [Incremental Gram-Schmidt] : " << k << endl;
      if (k > k_max)
      {
         k_max = k;
         VeryLong u;
         for (int j = 1; j <= k-1; j++)
         {
            int posj = p[j-1];
            if (posj == 0) continue;
            //u = dot(b[k-1], b[j-1]);
            u = b.dot(k-1, j-1);
            for (int i = 1; i <= posj-1; i++)
            {
               u = (d[i] * u - lambda(k-1, i-1) * lambda(j-1, i-1))/d[i-1];
            }
            lambda(k-1, posj-1) = u;
            //cout << "lambda[" << k << "][" << posj << "] = " << u << endl;
         }

         //u = dot(b[k-1], b[k-1]);
         u = b.dot(k-1, k-1);

         //cout << "u = " << u << endl;

         for (int i = 1; i <= s; i++)
         {
            u = (d[i] * u - lambda(k-1, i-1) * lambda(k-1, i-1))/d[i-1];
            if (u < 0L)
            {
               //cout << "u = " << u << endl;
            }
         }
         if (u == 0L)
         {
            p[k-1] = 0;
         }
         else
         {
            s++;
            p[k-1] = s;
            //cout << "d[s = " << s << "] <- " << u << endl;
            d[s] = u;
         }
      }
      if (k == 1)
      {
         k++;
//         cout << "(1) k changed up : " << k << endl;
         continue;
      }

      // step 3. [Test LLL condition]
      //cout << "step 3. [Test LLL condition]" << endl;
      // sub-algorithm REDI in line

      REDI_2(k, k-1, b, p, lambda, d);
      if (p[k-2] != 0 &&
            (p[k-1] == 0 ||
             (4L * d[p[k-1]] * d[p[k-1]-2] < 3L * d[p[k-1]-1] * d[p[k-1]-1] - 4L * lambda(k-1, p[k-1]-2) * lambda(k-1, p[k-1]-2))))
      {
         SWAPI_2(k, k_max, b, p, lambda, d);
         k--;
//         cout << "(1) k changed down : " << k << endl;
      }
      else
      {
         for (int l = k - 2; l >= 1; l--)
         {
            REDI_2(k, l, b, p, lambda, d);
         }
         k++;
//         cout << "(1) k changed up : " << k << endl;
      }
   }
   // step 4. [Finished]
   //cout << "step 4. [Finished] : s = " << s << endl;
   return s;
}

int LLL_reduce_2_on_columns(Matrix<VeryLong>& b)
{
   Matrix<VeryLong> bt = b.transpose();
   int s = LLL_reduce_2(bt);
   b = bt.transpose();
   return s;
}

int LLL_reduce_3_on_columns(Matrix<VeryLong>& b)
{
   // NOTE: this does LLL reduction on columns of b
   // generalisation to a set of not necessarily independent
   // vectors b
   // acknowledgements to NTL
   // step 1. [Initialize]
   //cout << "LLL step 1. [Initialize]" << endl;
   int k = 1;
   int k_max = 0;
   int n = b.columns();

   // corresponds to P in NTL LLL.cpp
   std::vector<int> p;
   p.resize(n);
   int s = 0;

   std::vector<VeryLong> d;
   d.resize(n+1);
   d[0] = 1L;
   //cout << "d[0] <- " << d[0] << endl;

   Matrix<VeryLong> lambda(n, VeryLong(0L), 0);

   while (k <= n)
   {
      // step 2. [Incremental Gram-Schmidt]
      //cout << "step 2. [Incremental Gram-Schmidt] : " << k << endl;
      if (k > k_max)
      {
         k_max = k;
         VeryLong u;
         for (int j = 1; j <= k-1; j++)
         {
            int posj = p[j-1];
            if (posj == 0) continue;
            //u = dot(b[k-1], b[j-1]);
            u = b.dot_on_columns(k-1, j-1);
            for (int i = 1; i <= posj-1; i++)
            {
               u = (d[i] * u - lambda(k-1, i-1) * lambda(j-1, i-1))/d[i-1];
            }
            lambda(k-1, posj-1) = u;
            //cout << "lambda[" << k << "][" << posj << "] = " << u << endl;
         }

         //u = dot(b[k-1], b[k-1]);
         u = b.dot_on_columns(k-1, k-1);

         //cout << "u = " << u << endl;

         for (int i = 1; i <= s; i++)
         {
            u = (d[i] * u - lambda(k-1, i-1) * lambda(k-1, i-1))/d[i-1];
            if (u < 0L)
            {
               //cout << "u = " << u << endl;
            }
         }
         if (u == 0L)
         {
            p[k-1] = 0;
         }
         else
         {
            s++;
            p[k-1] = s;
            //cout << "d[s = " << s << "] <- " << u << endl;
            d[s] = u;
         }
      }
      if (k == 1)
      {
         k++;
//         cout << "(1) k changed up : " << k << endl;
         continue;
      }

      // step 3. [Test LLL condition]
      //cout << "step 3. [Test LLL condition]" << endl;
      // sub-algorithm REDI in line

      REDI_2_on_columns(k, k-1, b, p, lambda, d[p[k-2]]);
      if (p[k-2] != 0 &&
            (p[k-1] == 0 ||
             (4L * d[p[k-1]] * d[p[k-1]-2] < 3L * d[p[k-1]-1] * d[p[k-1]-1] - 4L * lambda(k-1, p[k-1]-2) * lambda(k-1, p[k-1]-2))))
      {
         SWAPI_2_on_columns(k, k_max, b, p, lambda, d);
         k--;
//         cout << "(1) k changed down : " << k << endl;
      }
      else
      {
         for (int l = k - 2; l >= 1; l--)
         {
            REDI_2_on_columns(k, l, b, p, lambda, d[p[l-1]]);
         }
         k++;
//         cout << "(1) k changed up : " << k << endl;
      }
   }
   // step 4. [Finished]
   //cout << "step 4. [Finished] : s = " << s << endl;
   return s;
}

static VeryLong dot(const Matrix<VeryLong>& b, int col1, int col2)
{
   VeryLong result = 0L;
//   cout << "b[" << col1 << "] . b[" << col2 << "]" << endl;
   //cout << "b = " << endl << b;
   for (size_t i = 0; i < b.rows(); i++)
   {
      result += b(i, col1-1) * b(i, col2-1);
   }
   //cout << "result = " << result << endl;
   return result;
}

void REDI(int kk, int ll,
          Matrix<VeryLong>& basis,
          Matrix<VeryLong>& H,
          Matrix<VeryLong>& lambda,
          std::vector<VeryLong>& d)
{
//   cout << "REDI(" << kk << "," << ll << ")" << endl;
   VeryLong a = lambda(kk-1, ll-1);
   VeryLong b = d[ll];
   VeryLong tmp = a * 2L;
   if (tmp < 0L) tmp = -tmp;
   VeryLong tmp1 = b;
   if (tmp1 < 0L) tmp1 = -tmp1;
   if (tmp > b)
   {
      VeryLong q = (2L * a + b) / (2L * b);
      VeryLong q1 = a / b;
      VeryLong r = a % b;
      r = r + r;
      if (r > b || (r == b && q1 < 0L)) q1 += 1L;

      q = q1;
      for (size_t jj = 0; jj < H.rows(); jj++)
      {
         H(jj, kk-1) -= q * H(jj, ll-1);
      }
      for (size_t jjj = 0; jjj < basis.rows(); jjj++)
      {
         basis(jjj, kk-1) -= q * basis(jjj, ll-1);
      }
      lambda(kk-1, ll-1) -= q * d[ll];
      for (int ii = 1; ii <= ll - 1; ii++)
      {
         lambda(kk-1, ii-1) -= q * lambda(ll-1, ii-1);
      }
   }
}

void myswap(Matrix<VeryLong>& m, int col1, int col2)
{
   for (size_t row = 0; row < m.rows(); row++)
   {
      VeryLong tmp = m(row, col1 - 1);
      m(row, col1 - 1) = m(row, col2 - 1);
      m(row, col2 - 1) = tmp;
   }
}

void SWAPK(int k, int k_max, Matrix<VeryLong>& b, Matrix<VeryLong>& H, Matrix<VeryLong>& lambda, std::vector<VeryLong>& d, std::vector<int>& f)
{
//   cout << "SWAPK(" << k << "," << k_max << ")" << endl;
   myswap(H, k, k - 1);
   myswap(b, k, k - 1);
   if (k > 2)
   {
      for (int j = 1; j <= k - 2; j++)
      {
         VeryLong tmp = lambda(k-1, j-1);
         lambda(k-1, j-1) = lambda(k-2, j-1);
         lambda(k-2, j-1) = tmp;
      }
   }
   VeryLong lamb = lambda(k-1, k-2);
   if (lamb == 0L)
   {
      d[k-1] = d[k-2];
      f[k-2] = 0;
      f[k-1] = 1;
      lambda(k-1, k-2) = 0L;
      for (int i = k + 1; i <= k_max; i++)
      {
         lambda(i-1, k-1) = lambda(i-1, k-2);
         lambda(i-1, k-2) = 0L;
      }
   }
   else
   {
      for (int i = k + 1; i <= k_max; i++)
      {
         lambda(i-1, k-2) = (lamb * lambda(i-1, k-2)) / d[k-1];
      }
      VeryLong t = d[k];
      d[k-1] = (lamb * lamb) / d[k-1];
      d[k] = d[k-1];
      for (int j = k + 1; j <= k_max - 1; j++)
      {
         for (int i = j + 1; i <= k_max; i++)
         {
            lambda(i-1, j-1) = (lambda(i-1, j-1) * d[k-1]) / t;
         }
      }
      for (int j = k + 1; j <= k_max; j++)
      {
         d[j] = (d[j]*d[k-1])/t;
      }
   }
}

// Algorithm 2.7.2 (Kernel over Z Using LLL).
Matrix<VeryLong> kernel_over_Z_using_LLL(const Matrix<VeryLong>& AA)
{
   // NOTE: this does LLL reduction over columns of Z
   // Given an m x n matrix A with integral entries, this algorithm finds an LLL-reduced
   // Z-basis for the kernel of A. We use an auxiliary n x n integral matrix H.
   Matrix<VeryLong> A = AA;
   int n = A.columns();
   // step 1. [Initialize]
   int k = 2;
   int k_max = 1;
   std::vector<VeryLong> d;
   d.resize(n+1);
   d[0] = 1L;
   VeryLong t = dot(A, 1, 1);
   Matrix<VeryLong> H(n, VeryLong(1L), 0);
   std::vector<int> f;
   f.resize(n);
   Matrix<VeryLong> lambda(n, VeryLong(0L), 0);
   d[1] = 1L;
   f[0] = 0;
   if (t != 0L)
   {
      d[1] = t;
      f[0] = 1;
   }

   // step 2. [Incremental Gram-Schmidt]
   while (k <= n)
   {
      if (k > k_max)
      {
         k_max = k;
         for (int j = 1; j <= k; j++)
         {
            if (f[j-1] == 0 && j < k)
            {
               lambda(k-1, j-1) = 0L;
            }
            else
            {
               VeryLong u = dot(A, k, j);
               for (int i = 1; i <= j - 1; i++)
               {
                  if (f[i-1] != 0)
                  {
                     u = d[i] * u - lambda(k-1, i-1) * lambda(j-1, i-1);
                     if (u % d[i-1] != 0L)
                     {
                        std::cout << "Problem: u = " << u << ", d[i-1] = " << d[i-1] << std::endl;
                     }
                     u = u / d[i-1];
                  }
               }
               if (j < k)
               {
                  lambda(k-1, j-1) = u;
               }
               else if (j == k)
               {
                  if (u != 0L)
                  {
                     d[k] = u;
                     f[k-1] = 1;
                  }
                  else
                  {
                     d[k] = d[k-1];
                     f[k-1] = 0;
                  }
               }
            }
         }
      }

      // step 3. [Test f[k] == 0 and f[k-1] != 0]
      int done = 0;
      while (!done)
      {
         if (f[k-2] != 0)
         {
            REDI(k, k-1, A, H, lambda, d);
         }
         if (f[k-2] != 0 && f[k-1] == 0)
         {
            SWAPK(k, k_max, A, H, lambda, d, f);
            if (k > 2) k--;
            //cout << "k changed down: " << k << endl;
         }
         else done = 1;
      }
      for (int l = k - 2; l >= 1; --l)
      {
         if (f[l-1] != 0)
         {
            REDI(k, l, A, H, lambda, d);
         }
      }
      k++;
   }

   // step 4. [Finished?]
   int i = 1;
   while (i <= n && f[i-1] == 0) i++;
   int r = i - 1;

   Matrix<VeryLong> kernel(n, r);
   for (int row = 0; row < n; row++)
   {
      for (int col = 0; col < r; col++)
      {
         kernel(row, col) = H(row, col);
      }
   }
   Matrix<VeryLong> H1(r, r);

   LLL_reduce_1_on_columns(kernel,H1);

   return kernel;
}
