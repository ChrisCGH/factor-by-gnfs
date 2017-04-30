#pragma GCC diagnostic ignored "-Wundefined-var-template"
#include "NumberField.h"
#include "AlgebraicNumber.h"
#include "AlgebraicNumber_in_O_pO.h"
#include "Ideal.h"
#include "LongModular.h"
#include "pow.h"
#include "crt.h"
#include "lll.h"
#include <vector>
#include <deque>
#include <map>
#include <limits.h>
#include <time.h>
#include "timings.h"
#include "root.h"
#include "ExceptionalPrimes.h"
#include <iomanip>
//#define CHECKNORM 1
#ifdef CHECKNORM
#include "RootConfig.h"
#endif

PrimeList::iterator PrimeList::next_ = 0;
PrimeList::iterator PrimeList::last_ = 0;
std::vector<PrimeList::iterator> PrimeList::chunks_;

//#define CHECKCALC 1
#define CHECKCALC1 1
namespace
{
#ifdef CHECKCALC
double ln_norm_H;
double ln_norm_G;
#endif
#ifdef CHECKNORM
std::map<VeryLong, long int> check_norm_prime_decomposition;
std::map<long int, VeryLong> unfactored_norm_delta;
#endif
long long strtoll(const char* str)
{
   double x = std::atof(str);
   return (long long)x;
}

bool Debug = false;
Timing* timing = 0;
std::string Dump_file("");
ExceptionalPrimes* SpecialPrimes = 0;
std::map<long int, int> ProjectivePrimes;

typedef std::map<std::pair<long int, long int>, PrimeIdealRep*> NormalPrimesType;
//std::map<std::pair<long int, long int>, PrimeIdealRep*> NormalPrimes;
NormalPrimesType NormalPrimes;

typedef std::map<PrimeIdealRep*, int> PrimeIdealDecomposition;
const NumberField* nf;

void set_nf()
{
   if (!nf)
   {
      nf = &AlgebraicNumber::nf();
   }
}

void initialiseSpecialPrimes()
{
   set_nf();
   SpecialPrimes = new ExceptionalPrimes(nf);
}

typedef std::map<PrimeIdeal*, int> valuation_map_type;
valuation_map_type LeadingCoefficientValuationMap;

void findLeadingCoefficientValuations()
{
   // Calculate the valuations of <c_d> at each special prime
   std::vector<Quotient<VeryLong> > c;
   c.resize(nf->degree());
   VeryLong c_d = nf->c_d();
   c[0] = c_d;
   AlgebraicNumber an(c);
   Ideal I_c_d(an);
   VeryLong lcm(0L);
   Matrix<VeryLong> AA(1,1);
   Ideal::integralPart(I_c_d, lcm, AA);

   std::vector<long int> sp_list = SpecialPrimes->getListOfExceptionalPrimes();
   for (auto& p: sp_list)
   {
      for (std::vector<PrimeIdealRep*>::const_iterator piit = SpecialPrimes->begin(p);
            piit != SpecialPrimes->end(p);
            ++piit)
      {
         PrimeIdeal* pi = (*piit)->getPrimeIdeal();
         int val = PrimeIdeal::padicValuation(p, *pi, lcm, AA);
         LeadingCoefficientValuationMap[pi] = val;
         if (Debug)
         {
            std::cout << "Valuation of <c_d> is " << val << " at prime ideal " << *pi << std::endl;
         }
      }
   }
}

void initialiseProjectivePrimes()
{
   set_nf();
   VeryLong c_d = nf->c_d();
   std::vector<VeryLong> factors;
   c_d.factorise(&factors);
   VeryLong p;
   const VeryLong LONG_MAX_VL(LONG_MAX);
   for (size_t i = 0; i < factors.size(); i++)
   {
      p = factors[i];
      if (p < LONG_MAX_VL)
      {
         long int pp = p.get_long();
         ProjectivePrimes[pp]++;
      }
   }
}

int isSpecial(long int p)
{
   return SpecialPrimes->isExceptional(p);
}

void addSpecialPrimeFactorisation(const VeryLong& a,
                                  const VeryLong& b,
                                  long int p,
                                  int sign,
                                  PrimeIdealDecomposition& primeIdealProduct)
{
   VeryLong c_d = nf->c_d();
   for (std::vector<PrimeIdealRep*>::const_iterator specialPrimeIter = SpecialPrimes->begin(p);
         specialPrimeIter != SpecialPrimes->end(p);
         ++specialPrimeIter)
   {
      PrimeIdealRep* pir = *specialPrimeIter;
      PrimeIdeal* pi = pir->getPrimeIdeal();
      int val = SpecialPrimes->padicValuation(a * c_d, b * c_d, pi);
      if (val < 0)
      {
         // must use slow method
         std::vector<Quotient<VeryLong> > c;
         c.resize(nf->degree());
         c[0] = VeryLong(a);
         c[1] = -VeryLong(b);
         AlgebraicNumber an(c);
         Ideal I(an);
         VeryLong lcm(0L);
         Matrix<VeryLong> AA(1,1);
         Ideal::integralPart(I, lcm, AA);
         val = PrimeIdeal::padicValuation(p, *pi, lcm, AA);
      }
      else
      {
         val -= LeadingCoefficientValuationMap[pi];
      }
      val *= sign;
      if (val)
      {
         primeIdealProduct[*specialPrimeIter] += val;
         //std::cout << "(" << a << "," << b << ") : " << *pir << " : " << val << std::endl;
      }
   }
}

void addPrimeFactorisation(long int p, long int r, int e_p_r,
                           PrimeIdealDecomposition& primeIdealProduct)
{
   //std::cerr << "addPrimeFactorisation: p = " << p << ", r = " << r << ", e_p_r = " << e_p_r << std::endl;
   PrimeIdealRep* pi = 0;
   NormalPrimesType::iterator it = NormalPrimes.find(std::pair<long int, long int>(p, r));
   if (it == NormalPrimes.end())
   {
//      cout << "not found (" << p << "," << r << ")" << endl;
      pi = new PrimeIdealRep(p, r);
      NormalPrimes.insert(NormalPrimesType::value_type(std::make_pair(p, r), pi));
      //NormalPrimes[std::pair<long int, long int>(p, r)] = pi;
   }
   else
   {
      //pi = NormalPrimes[std::pair<long int, long int>(p, r)];
      pi = it->second;
//      cout << "found (" << p << "," << r << ")" << endl;
//      cout << "*pi" << endl;
   }
   primeIdealProduct[pi] += e_p_r;
}

void addNormalPrimeFactorisation(const VeryLong& a,
                                 long int b,
                                 long int p,
                                 int isProjective,
                                 int e_p_r,
                                 PrimeIdealDecomposition& primeIdealProduct)
{
   set_nf();
   LongModular::set_default_modulus(p);
   FactorBase& fb = nf->factorBase();
   //fb.add_extra(p);
   for (FactorBase::a_const_root_iterator rootIter = fb.begin(p);
         rootIter != fb.end(p);
         ++rootIter)
   {
      long int r = *rootIter;
      LongModular x(b % p);
      x *= LongModular(r);
      LongModular y((a % VeryLong(p)).get_long());
      if ((r == p && !isProjective) || x == y)
      {
         addPrimeFactorisation(p, r, e_p_r, primeIdealProduct);
      }
   }
}

void addRelationFactorisation(const Relation& rel, int sign,
                              PrimeIdealDecomposition& primeIdealProduct)
{
   VeryLong a = rel.a;
   long int b = rel.b;
   VeryLong c_d = nf->c_d();
   const VeryLong zero(0L);

   if (Debug)
   {
      std::cout << "addRelationFactorisation: (a,b) = (" << a << "," << b << ")" << std::endl;
   }
   for (auto& pr: rel.primes_)
   {
      long int p = pr.first;
      VeryLong p_vl(p);
      int e_p_r = pr.second;
//      cout << "p = " << p << ", e_p_r = " << e_p_r << endl;

      if (isSpecial(p))
      {
         addSpecialPrimeFactorisation(a, b, p, sign, primeIdealProduct);
      }
      else
      {
         int isProjective = (c_d % p_vl == zero);
         if (e_p_r != -1)
         {
            addNormalPrimeFactorisation(a, b, p, isProjective, e_p_r * sign, primeIdealProduct);
         }
         if (isProjective)
         {
            int v = -ProjectivePrimes[p];
            addPrimeFactorisation(p, p, v * sign, primeIdealProduct);
         }
      }
   }
}

void printPrimeIdealProduct(const PrimeIdealDecomposition& primeIdealProduct)
{
   //if (Debug)
   {
      std::fstream file("primeIdealProduct.dmp", std::ios::out);
      for (auto& pi: primeIdealProduct)
      {
         PrimeIdealRep* pir = pi.first;
         int v = pi.second;
         file << "v = <" << v << ">, ";
		 file << *pir << std::endl;
      }
   }
}

//----------------------------------------------------------------------------
// decompose a quotient of products of relations
void producePrimeDecomposition(const RelationList& relationNumer,
                               const RelationList& relationDenom,
                               PrimeIdealDecomposition& primeIdealProduct)
{
   std::cout << "producing prime decomposition ..." << std::endl;
   //int firstTime = 1;
   for (auto& rn: relationNumer)
   {
      Relation* rel = rn;
      //cout << "Adding (" << rel->a << "," << rel->b << ")" << endl;
      addRelationFactorisation(*rel, 1, primeIdealProduct);
   }
   if (relationDenom.size() > 0)
   {
      for (auto& rd: relationDenom)
      {
         Relation* rel = rd;
         addRelationFactorisation(*rel, -1, primeIdealProduct);
      }
   }
   NormalPrimes.clear();
   ProjectivePrimes.clear();
   delete SpecialPrimes;
   LeadingCoefficientValuationMap.clear();
}

//----------------------------------------------------------------------------
// valuation of F(a,b) at a prime p

int e_p_r(Relation& rel, long int p, long int v, long int r)
{
   //long int a = (rel.a % VeryLong(p)).get_long();
   long int a = rel.a % p;
   long int b = rel.b % p;
   long long int br = (long long int)b * r;

   //if (r != p && LongModular(rel.b) * LongModular(r) == LongModular(a))
   if (r != p && br % p == a)
      return v;

   //if (r == p && rel.b % p == 0)
   if (r == p && b == 0)
      return v;

   return 0;
}

std::map<std::pair<long int, long int>, int> S;
//static std::map<long int, int> V_p;

//----------------------------------------------------------------------------
// Complexity of a given quotient of products of relations
// In this case, the complete set of relations is given
// and a corresponding vector of signs indicates whether a
// particular relation is in the numerator or the denominator

double complexity(const RelationList& relations,
                  const std::vector<char>& e)
{
   timing->start("complexity full");
   S.clear();
   const VeryLong zero(0L);
   FactorBase& fb = nf->factorBase();
   int i = 0;
   for (auto& rel: relations)
   {
      for (auto& p1: rel->primes_)
      {
         long int p = p1.first;
         long int v = p1.second;
         LongModular::set_default_modulus(p);
         VeryLong pp(p);
         for (FactorBase::a_const_root_iterator rootIter = fb.begin(p);
               rootIter != fb.end(p);
               ++rootIter)
         {
            long int r = *rootIter;
            S[std::pair<long int, long int>(p, r)] += e[i] * e_p_r(*rel, p, v, r);
         }
         if (nf->c_d() % pp == zero)
         {
            S[std::pair<long int, long int>(p, p)] -= e[i] * ProjectivePrimes[p];
         }
      }
      i++;
   }

   double C = 0.0;
   for (auto& s1: S)
   {
      long int p = (s1.first).first;
      C += ln(p) * fabs((double)s1.second);
   }
   timing->stop();
   return C;
}

//----------------------------------------------------------------------------
// incremental complexity when e[j] -> -e[j] and resetting -e[k] -> e[k]
// If k < 0, then only change affect of e[j]
// If j >= n, then only reset -e[k] -> e[k]

double complexity(const RelationList& relations,
                  const std::vector<char>& e,
                  int j, int k, double& C)
{
//   double C = 0.0;
   int tmp = 0;
   const VeryLong zero(0L);
   FactorBase& fb = nf->factorBase();
   std::pair<long int, long int> a_pair(0L, 0L);
   if (j < (int)e.size())
   {
      Relation* rel = relations[j];
      for (auto& p1: rel->primes_)
      {
         long int p = p1.first;
         long int v = p1.second;
         a_pair.first = p;
         LongModular::set_default_modulus(p);
         //VeryLong pp(p);
         for (FactorBase::a_const_root_iterator rootIter = fb.begin(p);
               rootIter != fb.end(p);
               ++rootIter)
         {
            long int r = *rootIter;
            a_pair.second = r;
            //S[std::pair<long int, long int>(p, r)] -= 2 * e[j] * e_p_r(*rel, p, r);
            tmp = e_p_r(*rel, p, v, r) << 1;
            if (e[j] < 0) tmp = -tmp;
            C -= ln(p) * fabs((double)S[a_pair]);
            S[a_pair] -= tmp;
            C += ln(p) * fabs((double)S[a_pair]);
         }
         if (nf->c_d() % p == 0L)
         {
            a_pair.second = p;
            //S[std::pair<long int, long int>(p, p)] += 2 * e[j] * ProjectivePrimes[p];
            tmp = ProjectivePrimes[p] << 1;
            if (e[j] < 0) tmp = -tmp;
            C -= ln(p) * fabs((double)S[a_pair]);
            S[a_pair] += tmp;
            C += ln(p) * fabs((double)S[a_pair]);
         }
      }
   }
   if (k >= 0)
   {
      Relation* rel = relations[k];
      for (auto& p1: rel->primes_)
      {
         long int p = p1.first;
         long int v = p1.second;
         a_pair.first = p;
         LongModular::set_default_modulus(p);
         //VeryLong pp(p);
         for (FactorBase::a_const_root_iterator rootIter = fb.begin(p);
               rootIter != fb.end(p);
               ++rootIter)
         {
            long int r = *rootIter;
            a_pair.second = r;
            //S[std::pair<long int, long int>(p, r)] += 2 * e[k] * e_p_r(*rel, p, r);
            tmp = e_p_r(*rel, p, v, r) << 1;
            if (e[k] < 0) tmp = -tmp;
            C -= ln(p) * fabs((double)S[a_pair]);
            S[a_pair] += tmp;
            C += ln(p) * fabs((double)S[a_pair]);
         }
         if (nf->c_d() % p == 0L)
         {
            //S[std::pair<long int, long int>(p, p)] -= 2 * e[k] * ProjectivePrimes[p];
            a_pair.second = p;
            tmp = ProjectivePrimes[p] << 1;
            if (e[k] < 0) tmp = -tmp;
            C -= ln(p) * fabs((double)S[a_pair]);
            S[a_pair] -= tmp;
            C += ln(p) * fabs((double)S[a_pair]);
         }
      }
   }

   return C;
}

//----------------------------------------------------------------------------
// print relations
int printRelations(const RelationList& relations)
{
   for (size_t i = 0; i < relations.size(); i++)
   {
      Relation* rel = relations[i];
      VeryLong a = rel->a;
      long int b = rel->b;
      if (Debug)
      {
         std::cout << "(" << a << "," << b << ") : ";
         for (auto& p1: rel->primes_)
         {
            long int p = p1.first;
            long int v = p1.second;
            if (v > 1)
               std::cout << p << "^" << v << " ";
            else if (v == 1)
               std::cout << p << " ";
            else
               std::cout << p << "^(-1)" << " ";
         }
         std::cout << std::endl;
      }
   }
   return 0;
}

// read relations from a file
bool readRelations(const char* filename, RelationList& numerRelations, RelationList& denomRelations)
{
   const VeryLong LONG_MAX_VL(LONG_MAX);
   set_nf();
   Polynomial<VeryLong> min_poly = nf->min_poly();
   FactorBase& fb = nf->factorBase();

   std::fstream infile(filename, std::ios::in);
   // file format should be
   // a b
   std::string str;
   int done = 0;
   int line = 0;
   bool numeratorRelation = true;
   const std::string NumeratorStr("Numerator");
   const std::string DenominatorStr("Denominator");
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
         if (NumeratorStr == buf) continue;
         if (DenominatorStr == buf)
         {
            if (numeratorRelation)
            {
               numeratorRelation = false;
            }
            continue;
         }
         // find space
         char* c = buf;
         while (c && *c && *c != ' ') c++;
         if (!c || !*c || !*(c+1))
         {
            std::cerr << "Problem: bad format in line " << line << ":" << std::endl;
            std::cerr << "[" << buf << "]" << std::endl;
            return false;
         }
         *c = '\0';
         c++;
         long long int aa = strtoll(buf);
         long int b = std::atol(c);
         Relation* rel = new Relation(aa, b);
         VeryLong bb(b);
         VeryLong f = min_poly.evaluate_homogeneous(aa, bb);
         if (Debug)
         {
            std::cout << "F(" << aa << ", " << bb << ") = " << f << std::endl;
         }
         if (f < VeryLong(0L)) f = -f;

         // read factors from file, if any
         std::map<long int, int> primes;
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
               primes[p]++;
               fb.queue_extra(p);
               VeryLong pp(p);
               f /= pp;
            }
         }

         bool smooth = true;
         std::vector<VeryLong> factors;
         if (f > VeryLong(1L))
         {
            //std::cerr << "Factorising remainder of F(" << aa << "," << bb << ") = " << f << std::endl;
            f.factorise(&factors);
         }
         for (auto& pp: factors)
         {
            if (pp > LONG_MAX_VL)
            {
               smooth = false;
            }
            f /= pp;
         }
         if (f > VeryLong(1L))
         {
            std::cerr << "Problem: F(" << aa << "," << bb << ") has unfactored remainder : " << f << std::endl;
         }
         if (smooth)
         {
            for (auto& pp: factors)
            {
               long int p = pp.get_long();
               primes[p]++;
               fb.queue_extra(p);
            }
            // now look at the primes which divide c_d
            for (auto& pp: ProjectivePrimes)
            {
               long int p = pp.first;
               // for projective primes which aren't already in rel->primes set
               // value to -1
               if (primes.find(p) == primes.end())
               {
                  primes[p] = -1;
               }
            }
            // copy primes into Relation pointed to by rel
            rel->primes_.reserve(primes.size());
            for (auto& p1: primes)
            {
               rel->primes_.push_back(PrimeValuation(p1.first, p1.second));
            }
            if (numeratorRelation)
            {
               numerRelations.push_back(rel);
            }
            else
            {
               denomRelations.push_back(rel);
            }
         }
         else delete rel;
      }
   }
   fb.add_extra_queue();
   return true;
}

//----------------------------------------------------------------------------
// take list of relations in numer (denom is initially empty) and
// optimize split into numer and denom so that complexity (as defined
// by complexity() above) is minimized
void optimizeRelationsQuotient(RelationList& numer, RelationList& denom)
{
   if (std::getenv("DONT_OPTIMIZE")) return;
   std::cout << "Optimizing relations quotient ..." << std::endl;
   // initialise vector e to all ones => everything is in numer
   std::vector<char> e;
   int n = numer.size();
   e.resize(n, char(1));
   std::vector<char> best_e = e;
   double C = 0.0;
   double minC = 0.0;
   double best_local_minC = 1e200;
   int min_i = -1;
   int x = 0;
   const double DEFAULT_CUTOFF = 137.0;
   double CUTOFF = DEFAULT_CUTOFF;
   if (std::getenv("COMPLEXITY_CUTOFF"))
   {
      CUTOFF = std::atof(std::getenv("COMPLEXITY_CUTOFF"));
   }

   int done = 0;
   int goes = 1;
   while (!done && goes)
   {
      // choose random e
      for (int i = 0; i < n; i++)
      {
         x = genrand() % 2;
         if (!x) x--;
         e[i] = x;
      }
      // ecalculate complexity of e from scratch
      C = complexity(numer, e);

      // search for local minimum complexity near to e
      int local_min_found = 0;
      double local_minC = 1e200;
      int retries = 0;
      const int MAX_RETRIES = 2;
      while (!local_min_found)
      {
         minC = C;
         timing->start("complexity loop");
         char* tmp = std::getenv("NUMBER_TO_TRY");
         int NUMBER_TO_TRY = 0;
         if (tmp) NUMBER_TO_TRY = std::atoi(tmp);
         if (NUMBER_TO_TRY == 0) NUMBER_TO_TRY = 100;
         if (NUMBER_TO_TRY > n) NUMBER_TO_TRY = n;
         int i = -1;
         int prev_i = -1;
         int reduced = 0;
         int j = 0;
         for (j = 0; !reduced && j < NUMBER_TO_TRY; j++)
         {
            while (i == prev_i) i = genrand() % n;
            C = complexity(numer, e, i, prev_i, C);
            if (C < minC)
            {
               minC = C;
               min_i = i;
               reduced = 1;
            }
            prev_i = i;
         }
         timing->stop();
         if (minC < local_minC)
         {
            local_minC = minC;
            if (Debug)
            {
               std::cout << "minC = " << local_minC << ", j = " << j << std::endl;
            }
            if (min_i != prev_i) C = complexity(numer, e, min_i, prev_i, C);
            e[min_i] = -e[min_i]; // e and S are now in line
         }
         else
         {
            // e and S are out of line, since S reflects e with e[n-1] changed,
            // so change back:
            complexity(numer, e, n, prev_i, C);
            retries++;
            if (retries >= MAX_RETRIES)
            {
               local_min_found = 1;
               if (Debug)
               {
                  std::cout << "local_minC = " << local_minC << std::endl;
               }
            }
         }
      } // end local_min_found
      if (local_minC < best_local_minC)
      {
         best_local_minC = local_minC;
         // when we get here, e reflects local minimum,
         // and S is now in line with it
         best_e = e;
         if (Debug)
         {
            std::cout << "best_local_minC = " << best_local_minC << std::endl;
         }
      }

      if (best_local_minC < CUTOFF)
      {
         done = 1;
         if (Debug)
         {
            std::cout << "finished optimization" << std::endl;
            std::cout << "best_local_minC = " << best_local_minC << std::endl;
         }
         C = complexity(numer, best_e);
         if (Debug)
         {
            std::cout << "check: " << C << std::endl;
            ;
         }
      }
      goes--;
   }

   RelationList tmp;
   for (int i = 0; i < n; i++)
   {
      if (best_e[i] == -1)
      {
         denom.push_back(numer[i]);
      }
      else
      {
         tmp.push_back(numer[i]);
      }
   }
   numer.clear();
   numer = tmp;
   //cout << "norm = " << norm(numer) / norm(denom) << endl;
}

VeryLong LLL_max("100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");

Ideal selectIdeal(int s_l, PrimeIdealDecomposition& G, std::deque<PrimeIdealRep*>& G_index,
                  std::deque<PrimeIdealRep*>::reverse_iterator& pos_iter,
                  std::deque<PrimeIdealRep*>::iterator& neg_iter,
                  const Ideal& H,
                  const Quotient<VeryLong>& H_norm,
                  std::vector<PrimeIdealRep*>& contributingPrimes)
{
//   cout << "selectIdeal: s_l = " << s_l << ", H = " << H << endl;
   Ideal I(H);
   VeryLong n = H_norm.numerator();
//   cout << "N(H) = " << n << endl;
   //const VeryLong one(1L);
   if (H_norm.denominator() != 1L)
   {
      std::cerr << "Problem: N(I) not integer: " << I.norm() << std::endl;
   }

   if (s_l == 1)
   {
      // pick primes from the numerator of G,
      // i.e. from the end of G_index
      //std::deque<PrimeIdealRep*>::reverse_iterator iter = G_index.rbegin();
      std::deque<PrimeIdealRep*>::reverse_iterator& iter = pos_iter;
      int done = 0;
      while (iter != G_index.rend() && !done)
      {
         PrimeIdealRep* pir = *iter;
         int v = G[pir];
         if (v > 0)
         {
//	    cout << "pir = " << *pir << ", v = " << v << endl;
            VeryLong nm = pir->norm();
            int vv = 0;
            while (n * nm < LLL_max && vv < v)
            {
               n *= nm;
               vv++;
               //I *= *(pir->getPrimeIdeal());
               I *= *pir;
               G[pir]--;
               contributingPrimes.push_back(pir);
//	       cout << "n = " << n << ", LLL_max = " << LLL_max << endl;
            }
            if (vv < v) done = 1;
            else
            {
               pir->clearPrimeIdeal();
               ++iter;
            }
         }
         else done = 1;
      }
   }
   else // s_l == -1
   {
      // pick primes from the denominator of G,
      // i.e. from the start of G_index
      //std::deque<PrimeIdealRep*>::iterator iter = G_index.begin();
      std::deque<PrimeIdealRep*>::iterator& iter = neg_iter;
      int done = 0;
      while (iter != G_index.end() && !done)
      {
         PrimeIdealRep* pir = *iter;
         int v = G[pir];
         if (v < 0)
         {
//	    cout << "pir = " << *pir << ", v = " << v << endl;
            v = -v;
            VeryLong nm = pir->norm();
            int vv = 0;
            while (n * nm < LLL_max && vv < v)
            {
               n *= nm;
               vv++;
               //I *= *(pir->getPrimeIdeal());
               I *= *pir;
               G[pir]++;
               contributingPrimes.push_back(pir);
//	       cout << "n = " << n << ", LLL_max = " << LLL_max << endl;
            }
            if (vv < v) done = 1;
            else
            {
               pir->clearPrimeIdeal();
               ++iter;
            }
         }
         else done = 1;
      }
   }

   return I;
}

VeryLong careful_exp(long double x)
{
   const long double OVERFLOW_LIMIT = 400.0;
   if (x < OVERFLOW_LIMIT)
   {
      x = exp(x);
      // Just in case GMP has problems converting a very small number less than 1.0
      if (x < 1.0) x = 0.0;
      return VeryLong(x);
   }
   long int z = (long int)(x / OVERFLOW_LIMIT) + 1L;
   long double y = x / (long double)z;
   y = exp(y) + 1;
   VeryLong tmp(y);
   tmp = pow<VeryLong, long int>(tmp, z);
   tmp += z;
   return tmp;
}

AlgebraicNumber* selectDelta(const Ideal& I, long double ln_norm,
                             const std::vector<long double>& sigma, int s_l)
{
//   cout << "selectDelta : s_l = " << s_l << endl;
   int degree = nf->degree();
   Matrix<Quotient<VeryLong> > V(nf->w());

   const VeryLong one(1L);
   const VeryLong zero(0L);
   Matrix<VeryLong> rb = I.reducedBasisOmega();

   // V gives vectors of reduced basis in terms of alpha
   V *= rb;
//   cout << "V = " << endl << V;

   // Calculate c
   VeryLong I_norm = I.norm().numerator();

//   cout << "N(I) = " << I_norm << endl;
//   cout << "ln(LLL_max) = " << ln(LLL_max) << endl;
//   cout << "ln(N(I)) = " << ln(I_norm) << endl;
//   cout << "ln(N(gamma_l)) = " << ln_norm << endl;
   VeryLong fd = abs(nf->fieldDiscriminant());

//   cout << "ln(disc(K)) = " << ln(fd) << endl;

   long double ln_c = ln(LLL_max) - ln(I_norm) +
                      0.5 * (s_l * ln_norm - ln(fd));
   ln_c /= (long double)degree;
//   cout << "ln(c) = " << ln_c << endl;

   // Calculate lambda
   std::vector<long double> ln_lambda;
   ln_lambda.resize(degree);
   //std::vector<long double> lamf;
   //lamf.resize(degree);
   for (int j = 0; j < degree; j++)
   {
      long double ln_lamf = ln_c - 0.5 * s_l * sigma[j];
      ln_lambda[j] = ln_lamf;
//      cout << "ln_lamf[" << j << "] = " << ln_lamf << endl;
      //lamf[j] = exp((long double)ln_lamf);
      //cout << "lamf[" << j << "] = " << lamf[j] << endl;
   }

   Matrix<VeryLong> M(degree*2, degree);
   // first degree rows are filled by columns of reduced basis
   for (int row = 0; row < degree; row++)
   {
      for (int col = 0; col < degree; col++)
      {
         M(row,col) = rb(row,col);
      }
   }

   // remaining rows are of form lambda[row + degree] * v(alpha[row])
   // but if some alpha are complex (not pure real) then we must replace v() by sqrt(2)Re(v) or sqrt(2)Im(v)
   for (int col = 0; col < degree; col++)
   {
      int row = 0;
      while (row < degree)
      {

         // calculate v[col](alpha(row))
         std::vector<Quotient<VeryLong> > co;
         co.resize(degree);
         for (int i = 0; i < degree; i++)
         {
            co[i] = V(i,col);
         }
         AlgebraicNumber an(co);
         long double ln_re;
         long int re_sign;
         long double ln_im;
         long int im_sign;
         an.ln_sigma(row, ln_re, re_sign, ln_im, im_sign);
         complex<long double> alp = nf->conjugate(row);
//	 cout << row << "th conjugate is " << alp << endl;
         if (alp.imag() != (long double)0.0)
         {
            const long double ln_sqrt2 = log((long double)sqrt(2.0));
            long double tmp = ln_sqrt2 + ln_re + ln_lambda[row];
//	    cout << "ln(lambda_" << row << " * Re(sigma_" << row << "(v_" << col << ")) * sqrt(2)) = " << tmp << endl;
            VeryLong tmp1 = careful_exp(tmp);
            M(row + degree, col) = tmp1*re_sign;
            row++;
            tmp = ln_sqrt2 + ln_im + ln_lambda[row];
            tmp1 = careful_exp(tmp);
//	    cout << "ln(lambda_" << row << " * Im(sigma_" << row << "(v_" << col << ")) * sqrt(2)) = " << tmp << endl;
            M(row + degree, col) = tmp1*im_sign;
         }
         else
         {
            long double tmp = ln_re + ln_lambda[row];
//	    cout << "ln(lambda_" << row << " * (sigma_" << row << "(v_" << col << "))) = " << tmp << endl;
            VeryLong tmp1 = careful_exp(tmp);
            M(row + degree, col) = tmp1*re_sign;
         }
         row++;
      }
   }

   LLL_reduce_3_on_columns(M);

   // M gives basis in terms of integral basis
   const std::vector<AlgebraicNumber>& omega = AlgebraicNumber::integralBasis();
   AlgebraicNumber* delta = new AlgebraicNumber(0L);
   for (int i = 0; i < degree; i++)
   {
      *delta += omega[i] * M(i,0);
   }
   return delta;
}

//----------------------------------------------------------------------------
// main loop that approximates the square root

void processApproximation(const RelationList& relationNumer,
                          const RelationList& relationDenom,
                          const std::vector<long int>& good_primes,
                          const std::vector<AlgebraicNumber*>& delta,
                          const std::vector<int>& s,
                          const std::vector<VeryLong>& H_norms,
                          const std::vector<double>& ln_contributing_norms,
                          AlgebraicNumber& root_gamma_L)
{
#ifdef CHECKCALC
   double ln_norm_gamma_L = 0.0;
   double ln_norm_G = 0.0;
#endif
   int degree = nf->degree();
   std::fstream* dumpfile = 0;
   if (!Dump_file.empty())
   {
      dumpfile = new std::fstream(Dump_file.c_str(), std::ios::out);
   }

   if (dumpfile)
   {
      *dumpfile << "GOOD_PRIMES" << std::endl;
      for (size_t i = 0; i < good_primes.size(); ++i)
      {
         *dumpfile << good_primes[i] << std::endl;
      }
   }

   std::cout << good_primes.size() << " good primes chosen" << std::endl;

   std::cout << "Processing good primes ..." << std::endl;
   std::vector<AlgebraicNumber_in_O_pO_1> reducedGamma;

   for (size_t i = 0; i < good_primes.size(); i++)
   {
      //VeryLong p(good_primes[i]);
      long int p(good_primes[i]);
      std::cout << "good_primes[" << i << "] = " << p << std::endl;
      AlgebraicNumber_in_O_pO_1::set_basis(p);
      AlgebraicNumber_in_O_pO::set_basis(p);
      // find product of algebraic numbers in relationNumer, modulo p
      AlgebraicNumber_in_O_pO_1 numerProduct(1L);
      AlgebraicNumber_in_O_pO_1 numerDeltaProduct(1L);
#ifdef CHECKCALC
      AlgebraicNumber an_checkNumer(1L);
      AlgebraicNumber an_checkDenom(1L);
#endif
      //const int max_relations = 100;
      //int relation_count = 0;

      if (i == 0)
         if (dumpfile) *dumpfile << "NUMERATOR_RELATIONS" << std::endl;

      for (auto& rel: relationNumer)
      {
         AlgebraicNumber_in_O_pO_1 a(rel->a, rel->b);
         if (i == 0)
         {
            if (dumpfile) *dumpfile << rel->a << " " << rel->b << std::endl;
#ifdef CHECKCALC
            AlgebraicNumber aa(rel->a, rel->b);
            Quotient<VeryLong> norm_a = aa.norm();
            VeryLong num = norm_a.numerator();
            VeryLong den = norm_a.denominator();
            if (num < 0L) num = -num;
            if (den < 0L) den = -den;
            ln_norm_gamma_L += ln(num) - ln(den);
            an_checkNumer *= aa;
#endif
         }

         numerProduct *= a;
#ifdef CHECKCALC
         if (i == 0)
         {
            AlgebraicNumber_in_O_pO_1 check(an_checkNumer);
            if (check != numerProduct)
            {
               std::cout << "Problem!" << std::endl;
               std::cout << "   check = " << check << std::endl;
               std::cout << "   numerProduct = " << numerProduct << std::endl;
            }
         }
#endif
      }
      // find product of algebraic numbers in relationDenom, modulo p
      AlgebraicNumber_in_O_pO_1 denomProduct(1L);
      AlgebraicNumber_in_O_pO_1 denomDeltaProduct(1L);
      //int relation_count = 0;

      if (i == 0)
         if (dumpfile) *dumpfile << "DENOMINATOR_RELATIONS" << std::endl;

      for (auto& rel: relationDenom)
      {
         AlgebraicNumber_in_O_pO_1 a(rel->a, rel->b);
         if (i == 0)
         {
            if (dumpfile) *dumpfile << rel->a << " " << rel->b << std::endl;
#ifdef CHECKCALC
            AlgebraicNumber aa(rel->a, rel->b);
            Quotient<VeryLong> norm_a = aa.norm();
            VeryLong num = norm_a.numerator();
            VeryLong den = norm_a.denominator();
            if (num < 0L) num = -num;
            if (den < 0L) den = -den;
            ln_norm_gamma_L -= ln(num) - ln(den);
            an_checkDenom *= aa;
#endif
         }

         denomProduct *= a;
#ifdef CHECKCALC
         if (i == 0)
         {
            AlgebraicNumber_in_O_pO_1 check(an_checkDenom);
            if (check != denomProduct)
            {
               std::cout << "Problem!" << std::endl;
               std::cout << "   check = " << check << std::endl;
               std::cout << "   denomProduct = " << denomProduct << std::endl;
            }
         }
#endif
      }
#ifdef CHECKCALC
      if (i == 0)
      {
         ln_norm_G = ln_norm_gamma_L / 2;
         std::cerr << "ln(N(gamma ) = " << ln_norm_gamma_L << std::endl;
         std::cerr << "          0" << std::endl;
         std::cerr << "ln(N(G ) = " << ln_norm_G << std::endl;
         std::cerr << "          0" << std::endl;
      }
#endif

      // now go through the list of deltas (remember to square them)
      // and according to the sign of s[l], multiply them into numerProduct
      // or denomProduct:
      // if s[l] == 1, delta[i] goes into denomProduct
      // if s[l] == -1, delta[i] goes into numerProduct
      if (i == 0 && dumpfile)
         *dumpfile << "DELTA_LIST" << std::endl;

      for (size_t l = 0; l < delta.size(); l++)
      {
         if (i == 0 && dumpfile)
         {
            *dumpfile << l << "!" << s[l] << "!" << H_norms[l] << "!" << std::setprecision(16) << ln_contributing_norms[l] << std::endl << *delta[l] << std::endl;
         }

         if (s[l] > 0)
         {
            denomDeltaProduct *= *delta[l];
#ifdef CHECKCALC
            if (i == 0)
            {
               an_checkDenom *= *delta[l];
               an_checkDenom *= *delta[l];
               AlgebraicNumber_in_O_pO_1 check(an_checkDenom);
               if (check != denomProduct)
               {
                  std::cout << "Problem!" << std::endl;
                  std::cout << "   check = " << check << std::endl;
                  std::cout << "   denomProduct = " << denomProduct << std::endl;
               }
               Quotient<VeryLong> norm_delta = delta[l]->norm();
               VeryLong num = norm_delta.numerator();
               VeryLong den = norm_delta.denominator();
               if (den != 1L)
               {
                  std::cerr << "Problem! : delta[" << l << "] is not an algebraic integer" << std::endl;
                  std::cerr << "delta[" << l << "] = " << *delta[l] << std::endl;
                  std::cerr << "N(delta[" << l << "]) = " << norm_delta << std::endl;
               }
               if (num < 0L) num = -num;
               ln_norm_gamma_L -= ln(num);
               ln_norm_gamma_L -= ln(num);
               ln_norm_G -= ln_contributing_norms[l];
               if (l + 1 < s.size())
               {
                  std::cerr << "         s   = " << s[l + 1] << std::endl;
                  std::cerr << "          " << l + 1 << std::endl;
                  std::cerr << "ln(N(gamma )) = " << ln_norm_gamma_L << std::endl;
                  std::cerr << "          " << l + 1 << std::endl;
                  VeryLong H_norm = H_norms[l + 1];
                  if (H_norm < 0L) H_norm = -H_norm;
                  std::cerr << "    ln(N(H )) = " << ln(H_norm) << std::endl;
                  std::cerr << "          " << l + 1 << std::endl;
                  std::cerr << "    ln(N(G )) = " << ln_norm_G << std::endl;
                  std::cerr << "          " << l + 1 << std::endl;
                  double check = ln_norm_gamma_L / 2.0 - ln_norm_G;
                  if (s[l + 1] > 0)
                  {
                     check -= ln(H_norm);
                  }
                  else
                  {
                     check += ln(H_norm);
                  }
                  if (check < 0) check = -check;
                  const double epsilon = 1e-6;
                  if (check > epsilon)
                  {
                     std::cerr << "Problem! : check = " << check << std::endl;
                  }
               }
            }
#endif
         }
         else
         {
            numerDeltaProduct *= *delta[l];
#ifdef CHECKCALC
            if (i == 0)
            {
               an_checkNumer *= *delta[l];
               an_checkNumer *= *delta[l];
               AlgebraicNumber_in_O_pO_1 check(an_checkNumer);
               if (check != numerProduct)
               {
                  std::cout << "Problem!" << std::endl;
                  std::cout << "   check = " << check << std::endl;
                  std::cout << "   numerProduct = " << numerProduct << std::endl;
               }
               Quotient<VeryLong> norm_delta = delta[l]->norm();
               VeryLong num = norm_delta.numerator();
               VeryLong den = norm_delta.denominator();
               if (den != 1L)
               {
                  std::cerr << "Problem! : delta[" << l << "] is not an algebraic integer" << std::endl;
                  std::cerr << "delta[" << l << "] = " << *delta[l] << std::endl;
                  std::cerr << "N(delta[" << l << "]) = " << norm_delta << std::endl;
               }
               if (num < 0L) num = -num;
               ln_norm_gamma_L += ln(num);
               ln_norm_gamma_L += ln(num);
               ln_norm_G += ln_contributing_norms[l];
               if (l + 1 < s.size())
               {
                  std::cerr << "         s   = " << s[l + 1] << std::endl;
                  std::cerr << "          " << l + 1 << std::endl;
                  std::cerr << "ln(N(gamma )) = " << ln_norm_gamma_L << std::endl;
                  std::cerr << "          " << l + 1 << std::endl;
                  VeryLong H_norm = H_norms[l + 1];
                  if (H_norm < 0L) H_norm = -H_norm;
                  std::cerr << "    ln(N(H )) = " << ln(H_norm) << std::endl;
                  std::cerr << "          " << l + 1 << std::endl;
                  std::cerr << "    ln(N(G )) = " << ln_norm_G << std::endl;
                  std::cerr << "          " << l + 1 << std::endl;
                  double check = ln_norm_gamma_L / 2.0 - ln_norm_G;
                  if (s[l + 1] > 0)
                  {
                     check -= ln(H_norm);
                  }
                  else
                  {
                     check += ln(H_norm);
                  }
                  if (check < 0) check = -check;
                  const double epsilon = 1e-6;
                  if (check > epsilon)
                  {
                     std::cerr << "Problem! : check = " << check << std::endl;
                  }
               }
            }
#endif
         }
      }
      denomProduct *= denomDeltaProduct;
      denomProduct *= denomDeltaProduct;
      numerProduct *= numerDeltaProduct;
      numerProduct *= numerDeltaProduct;
      // we now have denomProduct and numerProduct, which are Algebraic numbers
      // expressed modulo p
      // we need to find the quotient of these
      AlgebraicNumber_in_O_pO_1 gamma_mod_p = numerProduct / denomProduct;
      std::cout << "gamma mod " << p << " = " << gamma_mod_p << std::endl;
      reducedGamma.push_back(gamma_mod_p);
   }

   if (dumpfile) *dumpfile << "FINAL_H" << std::endl;

   if (dumpfile) *dumpfile << H_norms[H_norms.size() - 1] << std::endl;

   if (dumpfile) *dumpfile << "END" << std::endl;
   if (dumpfile) delete dumpfile;

   std::cout << "Good primes processed" << std::endl;

#ifdef CHECKCALC
   std::cerr << "ln(N(H_L)) = " << ln_norm_H << std::endl;
   std::cerr << "ln(N(gamma_L)) = " << ln_norm_gamma_L << std::endl;
   double checknorm = ln_norm_H + ln_norm_H - ln_norm_gamma_L;
   if (checknorm < 0) checknorm = -checknorm;
   const double epsilon = 1e-6;
   if (checknorm > epsilon)
   {
      std::cerr << "Problem! : <gamma_L> != H_L^2" << std::endl;
   }

#endif
   const std::vector<AlgebraicNumber>& omega = AlgebraicNumber::integralBasis();

   std::vector<VeryLong> mm;
   VeryLong primeProduct(1L);
   for (size_t i = 0; i < good_primes.size(); i++)
   {
      VeryLong p(good_primes[i]);
      mm.push_back(p);
      primeProduct *= p;
   }
   AlgebraicNumber gamma_L(0L);
   std::vector<VeryLong> lifted_coefficients;
   size_t combinations = 1L;
   for (int i = 0; i < degree; i++)
   {
      std::vector<VeryLong> xx;
      for (size_t j = 0; j < reducedGamma.size(); j++)
      {
         LongModular v = reducedGamma[j].coefficient(i);
         xx.push_back(v.get_long());
      }
      VeryLong x = crt<VeryLong>(xx, mm);
      if (Debug)
      {
         std::cout << "Coefficient of omega_" << i + 1 << " = " << x << std::endl;
      }
      gamma_L += x * omega[i];
      lifted_coefficients.push_back(x);
      combinations <<= 1;
   }
#if 1
   // try each combination of +ve and -ve versions of lifted coefficients,
   // since we only know each lifted coefficient modulo primeProduct.
   // The correct combination should have a norm which is a perfect square.
   // There are 2^degree possible combinations to try ... (use Grey code to
   // iterate through them?)
   // mask[i] is bit i set, all others zero
   static size_t mask[10] =
   {
      0x0001, 0x0002, 0x0004, 0x0008, 0x0010, 0x0020, 0x0040, 0x0080, 0x0100, 0x0200
   };
   bool square_found = false;
   AlgebraicNumber best_gamma_L(0L);
   VeryLong best_norm(0L);
   for (size_t combination = 0; combination < combinations; ++combination)
   {
      gamma_L = AlgebraicNumber(0L);
      for (int i = 0; i < degree; ++i)
      {
         if (mask[i] & combination)
         {
            gamma_L += lifted_coefficients[i] * omega[i];
         }
         else
         {
            VeryLong x = lifted_coefficients[i];
            if (x < 0L) x += primeProduct;
            else x -= primeProduct;
            gamma_L += x * omega[i];
         }
      }
      Quotient<VeryLong> ng = gamma_L.norm();
      VeryLong n = ng.numerator();
      if (n < 0L) n = -n;
      if (n.is_square())
      {
         square_found = true;
         std::cout << "Square? gamma_L = " << gamma_L << std::endl;
         std::cout << "N(gamma_L) = " << ng << std::endl;
         if (best_norm == 0L || best_norm > n)
         {
            best_norm = n;
            best_gamma_L = gamma_L;
         }
      }
   }

   if (!square_found)
   {
      std::cout << "Problem! : couldn't find combination" << std::endl;
   }
   else
   {
      gamma_L = best_gamma_L;
   }
#endif
   // gamma_coeff gives coefficients of gamma_L in terms of omega
   Quotient<VeryLong> ng = gamma_L.norm();
//   cout << "N(gamma_L) = " << ng << endl;

   if (ng < Quotient<VeryLong>(0L)) ng = -ng;

   std::cout << "Calculating square root of gamma_L" << std::endl;
   std::cout << "ln(N(gamma_L)) = " << ln(ng.numerator()) - ln(ng.denominator()) << std::endl;

   root_gamma_L = gamma_L.sqrt();

   std::cout << "square root of gamma_L = " << root_gamma_L << std::endl;
   if (gamma_L != root_gamma_L * root_gamma_L)
   {
      std::cout << "Problem : root_gamma_L is *not* square root!" << std::endl;
      std::cout << "gamma_L = " << gamma_L << std::endl;
      std::cout << "root_gamma_L^2 = " << root_gamma_L * root_gamma_L << std::endl;
   }
}

void approximateSquareRoot(const RelationList& relationNumer,
                           const RelationList& relationDenom,
                           PrimeIdealDecomposition& primeIdealProduct,
                           std::vector<AlgebraicNumber*>& delta,
                           std::vector<int>& s,
                           AlgebraicNumber& root_gamma_L)
{
   std::cout << "approximating square root ..." << std::endl;
   int degree = nf->degree();
   LLL_max = nf->idealBound();
   LLL_max *= 10000L;
   s.clear();
   PrimeIdealDecomposition G;
   // G is a map between PrimeIdealRep* and the corresponding valuations.
   // We need to keep it as a map so we can update the valuations
   // easily, but we also want to be able to find the elements with largest
   // prime power norms, so we define a list of PrimeIdealRep* sorted by
   // norm, where the norm is given a -ve sign if valuation is negative
   std::deque<PrimeIdealRep*> G_index;

   // create H as <1>
   const VeryLong one(1L);
   const AlgebraicNumber an_one(one);
   Ideal H(an_one);
   delta.clear();
   delta.reserve(50000);
   std::vector<VeryLong> H_norms;
   Quotient<VeryLong> H_norm(1L);
   H_norms.push_back(one);
   std::vector<double> ln_contributing_norms;

   Ideal I;
   // Initialise variables
   // G is set to the square root of primeIdealProduct, i.e. divide each valuation by 2:
   for (auto& pip: primeIdealProduct)
   {
      PrimeIdealRep* pir = pip.first;
      int val = pip.second;
      if (val % 2 != 0)
      {
         std::cerr << "Problem: prime valuation " << val << " not even : " << *pir << std::endl;
      }
      G[pir] = val / 2;
      G_index.push_back(pir);
   }

   primeIdealProduct.clear();

   // calculate approximations to |sig_j(gamma)| where gamma is the algebraic number
   // given by relationNumer / relationDenom and sig_j is the jth embedding into C
   // given by the jth conjugate of the roots of alpha
   // These then let us approximate the norm of gamma, since N(gamma) = prod(sig_j(gamma))
   // or, taking logs, ln(N(gamma)) = sum(ln(sig_j(gamma)))
   std::vector<long double> sigma;
   sigma.resize(degree);
   for (int j = 0; j < degree; j++)
   {
      sigma[j] = 0.0;
   }
   for (auto& rel: relationNumer)
   {
      VeryLong a = rel->a;
      VeryLong b = rel->b;
      for (int j = 0; j < degree; j++)
      {
         sigma[j] += nf->ln_sigma(j, a, b);
      }
   }
   for (auto& rel: relationDenom)
   {
      VeryLong a = rel->a;
      VeryLong b = rel->b;
      for (int j = 0; j < degree; j++)
      {
         sigma[j] -= nf->ln_sigma(j, a, b);
      }
   }
   long double ln_norm = 0.0;
   for (int j = 0; j < degree; j++)
   {
	   std::cerr << "sigma[" << j << "] = " << sigma[j] << std::endl;
      ln_norm += sigma[j];
   }

   int s_l = 1;
   if (ln_norm < 0.0) s_l = -1;
#ifdef CHECKNORM
   std::cerr << "s_0 = " << s_l << std::endl;
   std::cerr << "ln_norm = " << ln_norm << std::endl;
#endif

   // Main approximation loop
   int done = 0;
   int retries = 0;
   int l = 0;

   std::sort (G_index.begin(), G_index.end(), 
              [&G](PrimeIdealRep* pir1, PrimeIdealRep* pir2)
               {
                   VeryLong n1 = pir1->norm();
                   int v1 = G[pir1];
                   VeryLong n2 = pir2->norm();
                   int v2 = G[pir2];
                   if (v1 < 0 && v2 > 0) return true;
                   if (v2 < 0 && v1 > 0) return false;
                   n1 *= (long int)v1;
                   n2 *= (long int)v2;
                   return (n1 < n2);
               }
             );
   // G_index is now sorted from large negative valuations to large positive valuations.
   // Maintain two iterators on G_index, pos_iter from the positive valuation end
   // and neg_iter from the negative valuation end
   std::deque<PrimeIdealRep*>::reverse_iterator pos_iter = G_index.rbegin();
   std::deque<PrimeIdealRep*>::iterator neg_iter = G_index.begin();

   while (!done)
   {
      if (Debug)
      {
         std::cout << "==================================================================================" << std::endl;
         std::cout << "approximation loop started for l = " << l << ", s_" << l << " = " << s_l << std::endl;
      }

      // select I_l and update G
      timing->start("selectIdeal");
      std::vector<PrimeIdealRep*> contributingPrimes;
      Ideal I = selectIdeal(s_l, G, G_index, pos_iter, neg_iter, H, H_norm, contributingPrimes);
#ifdef CHECKNORM
	  std::cerr << "l = " << l << " [";
      for (auto& cp: contributingPrimes)
	  {
		  std::cerr << *cp << ",";
	  }
	  std::cerr << "], N(I) = " << I.norm() << std::endl;
#endif

      timing->stop();

      // select delta from I
      timing->start("selectDelta");
      AlgebraicNumber* delta_l_ptr = selectDelta(I, ln_norm, sigma, s_l);
      AlgebraicNumber& delta_l = *delta_l_ptr;
      timing->stop();

      delta.push_back(delta_l_ptr);
      s.push_back(s_l);
#ifdef CHECKNORM
	  {
        VeryLong::setDebug();
		  VeryLong norm_delta = delta_l.norm().numerator();
          std::cerr << "delta_" << l << " = " << delta_l << std::endl;
          std::cerr << "N(delta_" << l << ") = " << norm_delta << std::endl;
		  if (norm_delta < 0L) norm_delta = -norm_delta;
		  for (size_t i = 0; i < contributingPrimes.size(); ++i)
		  {
			   VeryLong pp = contributingPrimes[i]->norm();
			   if (norm_delta % pp != 0L)
			   {
				   std::cerr << "Problem!: " << pp << " does not divide N(delta_" << l << ") : " << norm_delta << std::endl;
			   }
			   else
			   {
				   norm_delta /= pp;
			   }
		  }

		  std::vector<long int> factors;
		  VeryLong unfactored_part;
		  unfactored_norm_delta[l] = norm_delta;
	  }
#endif
      if (Debug)
      {
         std::cout << "delta_" << l << " = " << delta_l << std::endl;
         std::cout << "N(delta_" << l << ") = " << delta_l.norm() << std::endl;
      }

      // update sigma
      for (int j = 0; j < degree; j++)
      {
         sigma[j] -= 2.0 * s_l * delta_l.ln_sigma(j);
      }

      ln_norm = 0.0;
      for (int j = 0; j < degree; j++)
      {
         ln_norm += sigma[j];
      }

      // update H
      Ideal delta_l_ideal(delta_l);

#ifdef CHECKCALC
      /*
       * check that I / H is the product of the prime ideals we just removed from G in selectIdeal
       */
      Ideal check1 = I / H;
      Ideal check2(an_one);
      double ln_norm = 0.0;
      for (auto& pir: contributingPrimes)
      {
         check2 *= *pir;
         VeryLong norm = pir->norm();
         if (norm < 0L) norm = -norm;
         ln_norm += ln(norm);
         pir->clearPrimeIdeal();
      }
      ln_contributing_norms.push_back(ln_norm);
      if (check1 != check2)
      {
         std::cout << "Problem! : check1 != check2 :" << std::endl;
         std::cout << "I = " << I << std::endl;
         std::cout << "H = " << H << std::endl;
         std::cout << "I / H = " << check1 << std::endl;
         std::cout << "prod(contributingPrimes) = " << check2 << std::endl;
      }
#else
      if (!Dump_file.empty())
      {
         double ln_norm = 0.0;
         for (auto& pir: contributingPrimes)
         {
            VeryLong norm = pir->norm();
            if (norm < 0L) norm = -norm;
            ln_norm += ln(norm);
            pir->clearPrimeIdeal();
         }
         ln_contributing_norms.push_back(ln_norm);
      }
#endif

      H = delta_l_ideal / I;
#ifdef CHECKCALC
      /*
       * check that H * I = delta_l_ideal
       */
      Ideal check3 = H * I;
      if (check3 != delta_l_ideal)
      {
         std::cout << "Problem! : H * I != <delta> :" << std::endl;
         std::cout << "I = " << I << std::endl;
         std::cout << "H = " << H << std::endl;
         std::cout << "H * I = " << check3 << std::endl;
         std::cout << "delta_l_ideal = " << check3 << std::endl;
      }
#endif

      H_norm = H.norm();
      H_norms.push_back(H_norm.numerator());

      // update s_l
      s_l = -s_l;

      // finished when G = <1>
      // can check this by looking at first and last
      // elements of G_index
      PrimeIdealRep* first = *neg_iter;
      PrimeIdealRep* last = *pos_iter;
      if (G[first] == 0 && G[last] == 0)
      {
#ifdef CHECKCALC1
         /*
         * double check that G = <1>
         */
         for (auto& pid: G)
         {
            PrimeIdealRep* pir = pid.first;
            int v = pid.second;
            if (v != 0)
            {
               std::cerr << "Problem: G[pir] is non-zero (" << v << "), for pir = " << *pir << std::endl;
            }
         }
#endif
#ifdef CHECKCALC
         /*
          * check that <gamma_L> = H_L^2
         * or at least that ln(norm(gamma_L)) = 2 ln(norm(H_L)
          */
         Quotient<VeryLong> norm_H = H.norm();
         if (norm_H.denominator() != 1L)
         {
            std::cerr << "Problem: H_L is not an integral ideal : " << std::endl;
            std::cerr << "H_L = " << H << std::endl;
            std::cerr << "N(H_L) = " << norm_H << std::endl;
         }
         VeryLong norm_H_vl = norm_H.numerator();
         if (norm_H_vl < 0L) norm_H_vl = -norm_H_vl;
         ln_norm_H = ln(norm_H_vl);
#endif

		 if (retries <= 0 && s_l == 1)
		 {
			 done = 1;
             std::cout << "Done approximation loop" << std::endl;
		 }
		 else
		 {
			 ++l;
			 --retries;
		 }
      }
      else
      {
         ++l;
      }
   }

   // we now have a sequence of delta_l and s_l that define
   // a product formula for gamma_l.
   // we need to find the square root of this gamma_l by
   // expressing it as a polynomial with coefficents modulo
   // a number of 'good' primes.
   // Then use the CRT to reconstruct the formula for gamma_l.
   // A 'good' prime is a prime p such that f(X) is irreducible
   // mod p, and p does not divide N(H_i) for any I.
   // We need enough good primes, so that their product bounds
   // the coefficients of gamma_l (when written in terms of
   // powers of theta = c_d * alpha).
   FactorBase& fb = nf->factorBase();
   std::vector<int32_t>::const_reverse_iterator ip_iter = fb.rbegin_inert();
   std::vector<long int> good_primes;
   const int GOOD_PRIMES_NEEDED = 10;
   VeryLong product(1L);
   long int inert_p = 0L;
   std::cout << "Choosing good primes ..." << std::endl;

   while ((int)good_primes.size() < GOOD_PRIMES_NEEDED &&
          ip_iter != fb.rend_inert())
   {
      inert_p = *ip_iter;
      bool good = true;
      for (size_t i = 0; good && i < H_norms.size(); i++)
      {
         if (H_norms[i] % inert_p == 0L) good = false;
      }
      if (good)
      {
         good_primes.push_back(inert_p);
         product *= inert_p;
      }
      ++ip_iter;
   }
#ifdef CHECKNORM
      {
         if (!unfactored_norm_delta.empty())
         {
            Quotient<VeryLong> unfactored_quotient(1L);
	         std::cerr << "Contents of unfactored_norm_delta : " << std::endl;
             for (auto& und: unfactored_norm_delta)
	         {
               if (und.first % 2L == 0L)
               {
                  unfactored_quotient *= und.second;
               }
               else
               {
                  unfactored_quotient /= und.second;
               }
		         std::cerr << "l = " << und.first << ", unfactored part of N(delta_" << und.first << ") = " << und.second << std::endl;
	         }
            std::cerr << "unfactored_quotient = " << unfactored_quotient << std::endl;
         }
		   else
		   {
			   std::cerr << "unfactored_norm_delta is empty" << std::endl;
		   }
         // check that the primes in check_norm_prime_decomposition have positive even valuations
		   std::cerr << "contents of check_norm_prime_decomposition" << std::endl;
         for (auto& pd: check_norm_prime_decomposition)
         {
            VeryLong p = pd.first;
            long int v = pd.second;
            if (v % 2)
            {
               std::cerr << "Problem!: valuation of prime " << p << " in norm(gamma_L) is not even : " << v << std::endl;
            }
            if (v < 0)
            {
               std::cerr << "Problem!: valuation of prime " << p << " in norm(gamma_L) is negative : " << v << std::endl;
            }
			std::cerr << "p = " << p << ", v = " << v << std::endl;
         }
      }
#endif

   processApproximation(relationNumer, relationDenom, good_primes, delta, s, H_norms, ln_contributing_norms, root_gamma_L);
}

void readDump(const char* filename,
              std::vector<AlgebraicNumber*>& delta,
              std::vector<int>& s,
              AlgebraicNumber& root_gamma_L,
              RelationList& relationNumer, RelationList& relationDenom)
{
   std::vector<long int> good_primes;
   std::vector<VeryLong> H_norms;
   std::vector<double> ln_contributing_norms;
   {
      std::fstream dumpfile(filename, std::ios::in);

      std::string str;
      char tmp[1024];
      std::getline(dumpfile, str);

      if (str != "GOOD_PRIMES") return;

      std::getline(dumpfile, str);
      while (str != "NUMERATOR_RELATIONS")
      {
         // good primes
         long int p = std::atol(str.c_str());
         good_primes.push_back(p);
         std::getline(dumpfile, str);
      }

      std::getline(dumpfile, str);
      while (str != "DENOMINATOR_RELATIONS")
      {
         // numerator relations
         const char* c = str.c_str();
         char* d = tmp;
         while (*c != ' ')
         {
            *d = *c;
            c++;
            d++;
         }
         *d = '\0';
         long long int a = strtoll(tmp);

         d = tmp;
         while (*c)
         {
            *d = *c;
            c++;
            d++;
         }
         *d = '\0';
         long int b = std::atol(tmp);
         Relation* rel = new Relation(a, b);
         relationNumer.push_back(rel);
         std::getline(dumpfile, str);
      }

      std::getline(dumpfile, str);
      while (str != "DELTA_LIST")
      {
         // denominator relations
         const char* c = str.c_str();
         char* d = tmp;
         while (*c != ' ')
         {
            *d = *c;
            c++;
            d++;
         }
         *d = '\0';
         c++;
         long long int a = strtoll(tmp);

         d = tmp;
         while (*c)
         {
            *d = *c;
            c++;
            d++;
         }
         *d = '\0';
         long int b = std::atol(tmp);
         Relation* rel = new Relation(a, b);
         relationDenom.push_back(rel);
         std::getline(dumpfile, str);
      }

      std::getline(dumpfile, str);
      while (str != "FINAL_H")
      {
         // delta list
         const char* c = str.c_str();
         char* d = tmp;
         while (*c != '!')
         {
            *d = *c;
            c++;
            d++;
         }
         *d = '\0';
         c++;
         d = tmp;
         while (*c && *c != '!')
         {
            *d = *c;
            c++;
            d++;
         }
         *d = '\0';
         c++;
         int ss = std::atoi(tmp);
         s.push_back(ss);

         d = tmp;
         while (*c && *c != '!')
         {
            *d = *c;
            c++;
            d++;
         }
         *d = '\0';
         c++;
         VeryLong hn(tmp);
         H_norms.push_back(hn);

         d = tmp;
         while (*c)
         {
            *d = *c;
            c++;
            d++;
         }
         *d = '\0';
         double ln_cn = std::atof(tmp);
         ln_contributing_norms.push_back(ln_cn);

         AlgebraicNumber* del = new AlgebraicNumber;
         dumpfile >> *del;

         delta.push_back(del);

         std::getline(dumpfile, str);
      }

      std::getline(dumpfile, str);
      VeryLong H_norm(str);
#ifdef CHECKCALC
      if (H_norm < 0L) H_norm = -H_norm;
      ln_norm_H = ln(H_norm);
#endif
   }
   processApproximation(relationNumer, relationDenom, good_primes, delta, s, H_norms, ln_contributing_norms, root_gamma_L);
}
}

void squareRoot(const char* filename,
                std::vector<AlgebraicNumber*>& delta,
                std::vector<int>& s,
                AlgebraicNumber& root_gamma_L,
                RelationList& relationNumer, RelationList& relationDenom,
                bool debug,
                const std::string& dump_file)
{
   sgenrand(time(0));
   Debug = debug;
   Dump_file = dump_file;
   timing = new Timing("root.tim", false);
   initialiseSpecialPrimes();
   initialiseProjectivePrimes();
   findLeadingCoefficientValuations();
   PrimeIdealDecomposition primeIdealProduct;
   std::string str(filename);
   if (str.find(".dmp") == str.length() - 4)
   {
      readDump(filename, delta, s, root_gamma_L, relationNumer, relationDenom);
   }
   else
   {
      if (readRelations(filename, relationNumer, relationDenom))
      {
         std::cout << relationNumer.size() + relationDenom.size() << " relations read, #numerator relations = " << relationNumer.size() << ", #denominator relations = " << relationDenom.size() << std::endl;
         if (relationDenom.size() == 0)
         {
            optimizeRelationsQuotient(relationNumer, relationDenom);
         }
         if (Debug)
         {
            nf->factorBase().dump("factorbase.dmp");
            std::cout << "Numerator relations:" << std::endl;
            printRelations(relationNumer);
            std::cout << "Denominator relations:" << std::endl;
            printRelations(relationDenom);
         }
         producePrimeDecomposition(relationNumer, relationDenom, primeIdealProduct);
         nf->factorBase().clear();
         printPrimeIdealProduct(primeIdealProduct);
         approximateSquareRoot(relationNumer, relationDenom, primeIdealProduct, delta, s, root_gamma_L);
      }
   }
   timing->summary();

   delete timing;
}
