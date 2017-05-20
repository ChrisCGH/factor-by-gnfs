#include "NumberField.h"
#include "Ideal.h"
#include "pow.h"
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include "MPFloat.h"
#include <fstream>
#include "root.h"
#include "RootConfig.h"
#include <string>
#include <unordered_map>

namespace
{
long long strtoll(const char* str)
{
   double x = std::atof(str);
   return (long long)x;
}

int readRatRelations(const std::string& filename, RelationList& ratRelations)
{
   std::fstream infile(filename.c_str(), std::ios::in);
   int done = 0;
   int line = 0;
   std::string str;
   static char* buf = 0;
   static std::string::size_type buflen = 0;
   while (!done)
   {
      getline(infile, str);
      if (str.empty()) continue;
      if (str.size() > buflen)
      {
         delete [] buf;
         buflen = str.size();
         buf = new char [ buflen + 1 ];
      }
      strcpy(buf, str.c_str());
      if (buf[0] == '!') done = 1;
      else
      {
         line++;
         // find space
         char* c = buf;
         while (c && *c && *c != ' ') c++;
         if (!c || !*c || !*(c+1))
         {
            std::cerr << "Problem: bad format in line " << line << ":" << std::endl;
            std::cerr << "[" << buf << "]" << std::endl;
            return 0;
         }
         *c = '\0';
         c++;
         //VeryLong aa = VeryLong(buf);
         long long int aa = strtoll(buf);
         long int b = std::atol(c);
         //Relation* rel = new Relation(aa.get_long_long(), b);
         Relation* rel = new Relation(aa, b);
         std::unordered_map<long int, int> primes;
         //std::cout << "(aa,b) = (" << aa << "," << b << ")" << std::endl;
         // read factors from file, if any
         int done1 = 0;
         while (!done1)
         {
            // find next space or end of line
            while (c && *c && *c != ' ') c++;
            if (!c || !*c || !*(c+1)) done1 = 1;
            else
            {
               // we're at space before next factor
               c++;
               long int p = std::atol(c);
               //std::cout << "p = " << p << std::endl;
               primes[p]++;
            }
         }
         for (auto& pr: primes)
         {
            rel->primes_.push_back(PrimeValuation(pr.first, pr.second));
         }
         ratRelations.push_back(rel);
      }
   }
   return 1;
}
};

int main(int argc, char** argv)
{
   RootConfig config("root.cfg");
   config.display();
   std::string relfile_str = config.RELATION_FILE();
   std::string ratRelFile(relfile_str);
   const char* relfile = relfile_str.c_str();
   if (argc > 1) relfile = argv[1];

   Polynomial<VeryLong> f1 = config.f1();
   std::cout << "f1 = " << f1 << std::endl;

   for (int i = 0; i < config.EXTRA_PRIMES(); i++)
   {
      VeryLong::addPrime(config.EXTRA_PRIME(i));
   }
   char fbFile[132];
   strcpy(fbFile, config.ROOT_ID().c_str());
   strcat(fbFile, ".fb.dat");
   NumberField nf(f1, fbFile);

   AlgebraicNumber::setNumberField(nf);

   std::vector<AlgebraicNumber*> delta;
   std::vector<int> s;
   AlgebraicNumber root_gamma_L;
   RelationList relationNumer;
   RelationList relationDenom;
   bool debug = config.DEBUG();
   std::string dump_file = config.DUMP_FILE();
   squareRoot(relfile, delta, s, root_gamma_L, relationNumer, relationDenom, debug, dump_file);

   VeryLong N = config.N();
   VeryLong m = config.m();
   VeryLongModular::set_default_modulus(N);

   VeryLongModular phi_gamma = root_gamma_L.Phi(N, m);

   for (size_t l = 0; l < delta.size(); l++)
   {
      if (s[l] > 0)
      {
         phi_gamma *= delta[l]->Phi(N, m);
      }
      else
      {
         phi_gamma /= delta[l]->Phi(N, m);
      }
   }

   std::cout << "Phi(gamma) = " << phi_gamma << std::endl;

   // Now calculate Phi for the rational side
   // from relationNumer and relationDenom returned by squareRoot

   // read relations from relfile."rat"
   ratRelFile += ".rat";
   RelationList relationRat;
   readRatRelations(ratRelFile, relationRat);

   Polynomial<VeryLong> f2 = config.f2();
   std::vector<VeryLong> factors;
   std::unordered_map<long int, int> primes;
   auto numerIter = relationNumer.begin();
   auto denomIter = relationDenom.begin();
   auto ratIter = relationRat.begin();
   // ratRelations is in same order as relationNumer and relationDenom interspersed, e.g.
   //
   // 1 2 3 4 5 6 7 8 ...
   // 1     4 5   7   ...
   //   2 3     6   8 ...
   //
   // so we can do a kind of merge operation to iterate through

   while (ratIter != relationRat.end())
   {
      VeryLong a = (*ratIter)->a;
      VeryLong b = (*ratIter)->b;
      //std::cout << "(a,b) = (" << a << "," << b << ")" << std::endl;
      VeryLong f = abs(f2.evaluate_homogeneous(a, b));
      if (numerIter != relationNumer.end() &&
            a == VeryLong((*numerIter)->a) && b == (*numerIter)->b)
      {
         for (auto& p1: (*ratIter)->primes_)
         {
            long int p = p1.first;
            int v = p1.second;
            primes[p] += v;
            while (f % p == 0L)
            {
               f /= p;
               v--;
            }
            if (v != 0)
            {
               std::cout << "Problem: incorrect valuation : (a,b) = (" << a << "," << b << "), p = " << p << ", v = " << v << ", f = " << f << std::endl;
            }
         }
         if (f != 1L)
         {
            std::cout << "Problem: f(" << a << "," << b << ") not smooth, remaining quotient = " << f << std::endl;
         }
         ++numerIter;
      }
      else if (denomIter != relationDenom.end() &&
               a == VeryLong((*denomIter)->a) && b == (*denomIter)->b)
      {
         for (auto& p1: (*ratIter)->primes_)
         {
            long int p = p1.first;
            int v = p1.second;
            primes[p] -= v;
            while (f % p == 0L)
            {
               f /= p;
               v--;
            }
            if (v != 0)
            {
               std::cout << "Problem: incorrect valuation : (a,b) = (" << a << "," << b << "), p = " << p << ", v = " << v << ", f = " << f << std::endl;
            }
         }
         if (f != 1L)
         {
            std::cout << "Problem: f(" << a << "," << b << ") not smooth, remaining quotient = " << f << std::endl;
         }
         ++denomIter;
      }
      ++ratIter;
   }

   VeryLongModular phi_g(1L);
   VeryLongModular c_d(f2.coefficient(f2.deg()));
   long int S = relationNumer.size() - relationDenom.size();

   for (auto& pr: primes)
   {
      VeryLong p = pr.first;
      int v = pr.second;
      //std::cout << "(p,v) = (" << p << "," << v << ")" << std::endl;
      if (v % 2 != 0)
      {
         std::cerr << "Problem: rational side not a square, contains " << p << "^" << v << std::endl;
      }
      v /= 2;
      if (v > 0)
      {
         phi_g *= pow<VeryLongModular, int>(p, v);
      }
      else
      {
         phi_g /= pow<VeryLongModular, int>(p, -v);
      }
   }
   if (S > 0)
   {
      phi_g /= pow<VeryLongModular, long int>(c_d, S / 2L);
   }
   else
   {
      phi_g *= pow<VeryLongModular, long int>(c_d, -S / 2L);
   }
   std::cout << "Phi(g) = " << phi_g << std::endl;

   VeryLong u = phi_gamma.get_very_long();
   VeryLong v = phi_g.get_very_long();
   VeryLong w = u - v;
   if (w < 0L) w = -w;
   VeryLong f = ::gcd(w, N);
   std::cout << "f = " << f << std::endl;
   if (f > 1L && f < N)
   {
      std::cout << "Non-trivial factor found : " << f << std::endl;
      std::cout << "N % f = " << N % f << std::endl;
      VeryLong g = N / f;
      std::cout << "g = N / f = " << g << std::endl;
      std::cout << "N - fg = " << N - f * g << std::endl;
      return 0;
   }

   return 1;

}
