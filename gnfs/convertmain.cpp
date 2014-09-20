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
#include "convert.h"
// relation converter
// converts relations between different formats
namespace
{
MemoryMappedFile* mm_input_file = 0;
std::ostream* output_file = 0;
std::istream* input_file = 0;

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

enum format { FBGNFS = 0, GGNFS, NOT_DEFINED };
const char* format_string[] = { "FBGNFS", "GGNFS", 0 };

format string_to_format(const char* s)
{
   int i = 0;

   while (format_string[i] && strcmp(format_string[i], s) != 0) ++i;
   return static_cast<format>(i);
}

const char* format_to_string(format f)
{
   return format_string[static_cast<int>(f)];
}

format from_format(NOT_DEFINED);
format to_format(NOT_DEFINED);

void usage()
{
   std::cerr << "Usage: convert -from format -to format [-i[m] input] [-o output]" << std::endl;
   std::cerr << " where format is GGNFS or FBGNFS" << std::endl;
   std::cerr << " If -im is specified then the input is read as a memory-mapped file" << std::endl;
   std::exit(-1);
}

void init(int argc, char* argv[])
{
   if (argc != 5 && argc != 7 && argc != 9) usage();
   int arg = 1;
   while (arg < argc)
   {
      if (strcmp(argv[arg], "-from") == 0)
      {
         ++arg;
         if (arg >= argc) usage();
         from_format = string_to_format(argv[arg]);
         ++arg;
      }
      else if (strcmp(argv[arg], "-to") == 0)
      {
         ++arg;
         if (arg >= argc) usage();
         to_format = string_to_format(argv[arg]);
         ++arg;
      }
      else if (strcmp(argv[arg], "-im") == 0)
      {
         ++arg;
         if (arg >= argc) usage();
         mm_input_file = new MemoryMappedFile(argv[arg]);
         ++arg;
      }
      else if (strcmp(argv[arg], "-i") == 0)
      {
         ++arg;
         if (arg >= argc) usage();
         input_file = new std::fstream(argv[arg], std::ios::in);
         ++arg;
      }
      else if (strcmp(argv[arg], "-o") == 0)
      {
         ++arg;
         if (arg >= argc) usage();
         output_file = new std::fstream(argv[arg], std::ios::out);
         ++arg;
      }
   }
   if (from_format == NOT_DEFINED ||
         to_format == NOT_DEFINED) usage();
   if (!output_file) output_file = &std::cout;
   if (!mm_input_file && !input_file) input_file = &std::cin;

   std::cerr << "converting from " << format_to_string(from_format) << " to " << format_to_string(to_format) << std::endl;
}

long long strtoll(const char* str)
{
   double x = std::atof(str);
   return (long long)x;
}

void convert_from_FBGNFS_to_GGNFS(const std::string& input_string, std::ostream& os)
{
   // FBGNFS format for relation:
   // a b : p1/r1 p2/r2 .... pn/rn : q1 q2 ... qm :
   //
   // GGNFS format for relation:
   // a,b:q1,q2,...,qm:p1,p2,...,pn
   //
   long long a;
   long long b;
   char* alg_str;
   char* rat_str;
   Convert::parse_FBGNFS(input_string, a, b, alg_str, rat_str);

   std::vector<long int> alg_primes;
   std::vector<long int> rat_primes;
   Convert::extract_FBGNFS_primes(alg_str, alg_primes);
   Convert::extract_FBGNFS_primes(rat_str, rat_primes);

   static char buf[2048];
   char tmp[132];
#ifdef WIN32
   std::sprintf(buf, "%I64d,%I64d:", a, b);
#else
   std::sprintf(buf, "%lld,%lld:", a, b);
#endif

   for (size_t i = 0; i < rat_primes.size(); ++i)
   {
      std::sprintf(tmp, "%lx", static_cast<size_t>(rat_primes[i]));
      strcat(buf, tmp);
      if (i < rat_primes.size() - 1) strcat(buf, ",");
   }
   strcat(buf, ":");
   for (size_t i = 0; i < alg_primes.size(); ++i)
   {
      std::sprintf(tmp, "%lx", static_cast<size_t>(alg_primes[i]));
      strcat(buf, tmp);
      if (i < alg_primes.size() - 1) strcat(buf, ",");
   }
   os << buf << std::endl;
}

char* parse_GGNFS_relation(const std::string& str, long long int& a, long long int& b)
{
   copy_str(str);
   char* c = buf;
   char* d = c;
   while (d && *d && *d != ',') ++d;
   if (d && *d) *d = '\0';
   ++d;
   // str -> aaaaaa,bbbbbbb:alg:rat
   //        ^      ^
   //        c      d
   a = strtoll(c);
   c = d;
   while (d && *d && *d != ':') ++d;
   if (d && *d) *d = '\0';
   ++d;
   // str -> aaaaaa,bbbbbbb:alg:rat
   //               ^       ^
   //               c       d
   b = strtoll(c);
   return d;
}

void parse_GGNFS(const std::string& str, long long int& a, long long int& b, char*& alg_str, char*& rat_str)
{
   // str -> aaaaaa,bbbbbbb:alg:rat
   char* d = parse_GGNFS_relation(str, a, b);
   char* c = d;
   while (d && *d && *d != ':') ++d;
   if (d && *d) *d = '\0';
   ++d;
   // str -> aaaaaa,bbbbbbb:alg:rat
   //                       ^   ^
   //                       c   d

   rat_str = c;
   alg_str = d;
}

void extract_GGNFS_primes(char* str, std::vector<long int>& primes)
{
   primes.clear();
   // str is a set of strings separated by commas
   // s1,s2,s3
   char* c = str;
   char* d = c;
   while (*c)
   {
      while (d && *d && *d != ',') ++d;
      if (*d)
      {
         *d = '\0';
         ++d;
      }
      long int primeNum = std::strtol(c, 0, 16);
      primes.push_back(primeNum);
      c = d;
   }

   std::sort(primes.begin(), primes.end());
}

void convert_from_GGNFS_to_FBGNFS(const std::string& input_string, std::ostream& os)
{
   long long a;
   long long b;
   char* alg_str;
   char* rat_str;
   parse_GGNFS(input_string, a, b, alg_str, rat_str);
   std::vector<long int> alg_primes;
   std::vector<long int> rat_primes;
   extract_GGNFS_primes(alg_str, alg_primes);
   extract_GGNFS_primes(rat_str, rat_primes);

   static char buf[2048];
   char tmp[132];
#ifdef WIN32
   std::sprintf(buf, "%I64d %I64d : ", a, b);
#else
   std::sprintf(buf, "%lld %lld : ", a, b);
#endif

   for (size_t i = 0; i < alg_primes.size(); ++i)
   {
      long int p = alg_primes[i];
      long int b_inv = inverse<long int>((long int)b, p);
      long int r = 0;
      mulmodasm2(a, b_inv, p, r);
      std::sprintf(tmp, "%ld/%ld ", p, r);
      strcat(buf, tmp);
   }
   strcat(buf, ": ");
   for (size_t i = 0; i < rat_primes.size(); ++i)
   {
      std::sprintf(tmp, "%ld ", rat_primes[i]);
      strcat(buf, tmp);
   }
   strcat(buf, ":");
   os << buf << std::endl;
}

void cleanup()
{
   delete mm_input_file;
   if (output_file != &std::cout) delete output_file;
   if (input_file && input_file != &std::cin) delete input_file;
}
};

int main(int argc, char* argv[])
{
   init(argc, argv);

   std::string input_line;
   while ((mm_input_file && getline(*mm_input_file, input_line)) ||
          (input_file && std::getline(*input_file, input_line)))
   {
      if (from_format == FBGNFS && to_format == GGNFS)
      {
         convert_from_FBGNFS_to_GGNFS(input_line, *output_file);
      }
      else
      {
         convert_from_GGNFS_to_FBGNFS(input_line, *output_file);
      }
   }

   cleanup();
   return 0;
}
