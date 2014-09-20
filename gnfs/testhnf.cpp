#include <string>
#include <fstream>
#include "Matrix.h"
#include "VeryLong.h"
#include "timings.h"
Timing timing("testhnf.tim", false);

std::istream& operator>> (std::istream& is, Matrix<VeryLong>& m)
{
   std::string str;
   int i = 0;
   int j = 0;
   while (std::getline(is, str))
   {
      std::string::size_type pos = str.find(" ");
      while (pos != std::string::npos)
      {
         std::string s = str.substr(0, pos);
         str = str.substr(pos + 1);
         m(i,j) = VeryLong(s);
         pos = str.find(" ");
         ++j;
      }
      j = 0;
      ++i;
   }
   return is;
}

int main(int argc, char* argv[])
{
   const long int d = 4;
   Matrix<VeryLong> N(d, d*d);
   //Matrix<VeryLong> N(3, 4);
   std::fstream fin(argv[1], std::ios::in);
   fin >> N;
   std::cout << N << std::endl;
   //timing.start("HNF");
   Matrix<VeryLong> N1 = HNF(N);
   //timing.stop();
   Matrix<VeryLong> N2 = HNF1(N);
   VeryLong det = determinant_in_integral_domain(N1);
   Matrix<VeryLong> N3 = HNF_mod_D(N, det);
   //if (N1 != N2)
   {
      //std::cout << "N1, N2 differ!" << std::endl;
      std::cout << N1 << std::endl;
      std::cout << determinant_in_integral_domain(N1) << std::endl;
      std::cout << N2 << std::endl;
      std::cout << determinant_in_integral_domain(N2) << std::endl;
      std::cout << N3 << std::endl;
      std::cout << determinant_in_integral_domain(N3) << std::endl;
   }
//   else
//   {
//      std::cout << "N1, N2 identical!" << std::endl;
//   }
   for (int i = 0; i < 10000; ++i)
   {
      timing.start("HNF");
      Matrix<VeryLong> N2 = HNF(N);
      timing.stop();
      if (N1 != N2)
      {
         std::cout << "N1, N2 differ!" << std::endl;
      }
   }
   for (int i = 0; i < 10000; ++i)
   {
      timing.start("HNF1");
      Matrix<VeryLong> N2 = HNF1(N);
      timing.stop();
      if (N1 != N2)
      {
         std::cout << "N1, N2 differ!" << std::endl;
      }
   }
   det = determinant_in_integral_domain(N1);
   std::cout << "determinant(N1) = " << det << std::endl;
   for (int i = 0; i < 10000; ++i)
   {
      timing.start("HNF_mod_D");
      Matrix<VeryLong> N2 = HNF_mod_D(N, det);
      timing.stop();
      if (N1 != N2)
      {
         std::cout << "N1, N2 differ!" << std::endl;
      }
   }

   timing.summary();
   return 0;
}
