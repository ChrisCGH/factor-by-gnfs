#include "squfof.h"
#include <iostream>

int main()
{
   long long int N = 2757413891990958491LL;
   long int factor;
   if (SQUFOF(N, factor))
   {
      long int other_factor = N / factor;
      std::cout << N << " = " << factor << " * " << other_factor << std::endl;
   }

   return 0;
}
