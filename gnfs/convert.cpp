#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <string>
#include <algorithm>
#include <vector>
#include <iomanip>
#include "MemoryMappedFile.h"
#include "mod.h"
#include "gcd.h"
// relation converter
// converts relations between different formats
namespace
{
char* buf = 0;
std::string::size_type buflen = 0;
void copy_str(const std::string& str)
{
   if (str.empty()) return;
   if (str.size() > buflen)
   {
      buflen = str.size();
      delete [] buf;
      buf = new char [ buflen + 1 ];
   }
   strcpy(buf, str.c_str());
}

long long strtoll(const char* str)
{
   double x = std::atof(str);
   return (long long)x;
}
};

namespace Convert
{
char* parse_FBGNFS_relation(const std::string& str, long long int& a, long long int&b)
{
   copy_str(str);
   char* c = buf;
   char* d = c;
   while (d && *d && *d != ' ') ++d;
   if (d && *d) *d = '\0';
   ++d;
   // str -> aaaaaa bbbbbbb : alg : rat
   //        ^      ^
   //        c      d
   a = strtoll(c);
   c = d;
   while (d && *d && *d != ':') ++d;
   if (d && *d) *(d - 1) = '\0';
   ++d;
   ++d;
   // str -> aaaaaa bbbbbbb : alg : rat
   //               ^         ^
   //               c         d
   b = strtoll(c);
   return d;
}

void parse_FBGNFS(const std::string& str, long long int& a, long long int& b, char*& alg_str, char*& rat_str)
{
   // str -> aaaaaa bbbbbbb : alg : rat :
   char* d = parse_FBGNFS_relation(str, a, b);
   char* c = d;
   while (d && *d && *d != ':') ++d;
   if (d && *d) *(d - 1) = '\0';
   ++d;
   ++d;
   // str -> aaaaaa bbbbbbb : alg : rat :
   //                         ^     ^
   //                         c     d

   alg_str = c;
   c = d;
   while (d && *d && *d != ':') ++d;
   if (d && *d) *(d - 1) = '\0';
   rat_str = c;
}

void extract_FBGNFS_primes(char* str, std::vector<long int>& primes)
{
   primes.clear();
   // str is a set of strings separated by spaces
   // s1 s2 s3
   char* c = str;
   char* d = c;
   while (*c)
   {
      while (d && *d && *d != ' ') ++d;
      if (*d)
      {
         *d = '\0';
         ++d;
      }
      long int primeNum = std::atol(c);
      primes.push_back(primeNum);
      c = d;
   }

   std::sort(primes.begin(), primes.end());
}

};
