#include <time.h>
#include "Siever.h"

int main(int argc, char** argv)
{
   sgenrand(time(0));
   Siever siever;
   long int min_A = 1L;
   long int max_A = 30000L;
   if (argc > 2 && strcmp(argv[1], "-check") == 0)
   {
      siever.checkRelations(argv[2]);
      return 0;
   }
   if (argc > 1) min_A = std::atol(argv[1]);
   if (argc > 2) max_A = std::atol(argv[2]);
   if (max_A <= min_A)
   {
      if (min_A < 0) max_A = -min_A;
      else max_A = 2L * min_A;
   }
   long int B_span = 1L;
   long int B_first = 1L;
   long int B_last = B_first + 1L;
   if (argc > 3) B_first = std::atol(argv[3]);
   if (B_first < 1) B_first = 1L;
   if (argc > 4) B_last = std::atol(argv[4]);
   if (B_last <= B_first) B_last = B_first + 1L;
   if (argc > 5) B_span = std::atol(argv[5]);
   if (B_span < 1) B_span = 1;
   for (int B = B_first; B < B_last; B += B_span)
   {
      siever.sieve(min_A, max_A, B);
   }

   return 0;
}
