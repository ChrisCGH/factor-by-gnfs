#include <time.h>
#include "LatticeSiever.h"
#include <algorithm>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <sstream>

namespace
{
void get_saved_values(long int& min_q, long int& max_q)
{
   std::fstream f("savedq.txt", std::ios::in);
   std::string str;
   if (getline(f, str))
   {
      std::stringstream ss(str);
      ss >> min_q;
      ss >> max_q;
   }
}
void put_saved_values(long int min_q, long int max_q)
{
   std::fstream f("savedq.txt", std::ios::out);
   f << min_q << " " << max_q << std::endl;
}
LatticeSiever siever;
};

int main(int argc, char** argv)
{
   sgenrand(static_cast<unsigned long int>(time(0)));
   long int min_q = -1L;
   long int max_q = -1L;
   int arg = 1;
   bool sample = false;
   while (arg < argc)
   {
      if (strcmp(argv[arg], "-s") == 0)
      {
         sample = true;
         ++arg;
      }
      else if (isdigit(argv[arg][0]))
      {
         if (min_q < 0) min_q = std::atol(argv[arg]);
         else if (max_q < 0) max_q = std::atol(argv[arg]);
         ++arg;
      }
   }
   if (!sample)
   {
      if (min_q == -1L && max_q == -1L)
      {
         get_saved_values(min_q, max_q);
         if (min_q == -1L && max_q == -1L)
         {
            std::cerr << "no values for min_q/max_q specified" << std::endl;
            return 0;
         }
      }
      if (max_q <= min_q)
      {
         std::swap(min_q, max_q);
      }
      for (long int q = min_q; q < max_q; q++)
      {
         VeryLong q_vl(q);
         if (q_vl.is_probable_prime(10))
         {
            put_saved_values(q + 1, max_q);
            if (!siever.sieve(q))
            {
                //return 0;
            }
         }
      }
   }
   else
   {
      const int sample_count = 100;
      int sample_span = (max_q - min_q) / sample_count;
      long int q = min_q;
      while (q < max_q)
      {
         VeryLong q_vl(q);
         while (!q_vl.is_probable_prime(10)) q_vl += 1L;
         q = q_vl.get_long();
         std::cout << q << std::endl;
         siever.sieve(q);
         q += sample_span;
      }
   }

   return 0;
}
