#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <vector>
#include <list>
#include <string>
#include <memory>
#include <math.h>
#include <stdlib.h>
#include "graph.h"
#include "PriorityQueue.h"
#include "RelationManager.h"
#include "MemoryMappedFile.h"
#include "convert.h"
#include "Logger.h"

namespace
{
char* MemoryReserve = 0;
const size_t MemoryReserveSize = 1024*1024;
//const size_t MemoryReserveSize = 1024*1024*200;
bool NewHandlerCalled = false;
std::new_handler Old_newhandler;
void newhandler()
{
   if (!NewHandlerCalled)
   {
      delete [] MemoryReserve;
   }
   else
      {}
   NewHandlerCalled = true;
   std::set_new_handler(Old_newhandler);
}

// table of pi(x) to allow estimation of ExcessMin
struct pi_x_item
{
   long int p;
   long int pi_x;
};
const int pi_x_table_size = 8;
pi_x_item pi_x_table[pi_x_table_size] =
{
   {100, 25}, {1000, 168}, {5000, 669}, {10000, 1229},
   {50000, 5133}, {100000, 9592}, {500000, 41538}, {1000000, 78498}
};

long int estimateRequiredExcess(long int filtmin)
{
   int i = 0;
   while (i < pi_x_table_size && filtmin > pi_x_table[i].p) ++i;
   if (i < pi_x_table_size) return pi_x_table[i].pi_x * 3L;
   return 100000L;
}

std::string RelationsFile("");
std::ostream* Reloutfile = 0;
std::string Reloutfilename("");
std::ostream* Relsetoutfile = 0;
std::ostream* Relsetprimeoutfile = 0;
MemoryMappedFile* Relsetinmmfile = 0;
size_t RelationSetCount = 0;
MemoryMappedFile* Relsetprimeinmmfile = 0;
size_t MergeLevel = 1;
size_t MinMergeLevel = 3;
int FiltMin = 0;
size_t ExcessMin = 0;
int MaxPass = 0;
bool MergeOnly = false;
int MaxDiscard = 0;
long int RelationsMergedInThisPass = 0;
const long int CutoffForRelationsMergedInThisPass = 600000;

RelationTable* relationTable = 0;
struct Prime
{
   uint32_t p;
   uint32_t r;
   bool operator<(const Prime& prime) const
   {
      if (p < prime.p) return true;
      if (p == prime.p && r < prime.r) return true;
      return false;
   }
   bool operator!=(const Prime& prime) const
   {
      if (p != prime.p) return true;
      if (p == prime.p && r != prime.r) return true;
      return false;
   }
   bool operator==(const Prime& prime) const
   {
      if (p != prime.p) return false;
      if (p == prime.p && r != prime.r) return false;
      return true;
   }
   Prime(const char* c)
   {
      p = atol(c);
      r = 0;
      while (*c && *c != '/') ++c;
      if (*c)
      {
         r = atol(c+1);
      }
   }
   Prime() : p(0), r(0)
   {}
   struct Hasher
   {
      static uint32_t hash(const Prime& p)
      {
         return p.p;
      }
      std::size_t operator()(const Prime& p) const
      {
         return p.p;
      }
   };
};
typedef std::unordered_map<Prime, int, Prime::Hasher> prime_map_type;
int Excess = 0;

PrimeFrequencyTable* frequencyTable = 0;

typedef PriorityQueue<int, 50000> clique_queue_type;

unsigned long long int hash_relation(long long int a, long long int b)
{
   const long long int PI = 314159265358979323LL;
   const long long int E =  271828182845904523LL;

   unsigned long long int h = a * PI + b * E;
   return h;
}

struct RelationHasher
{
   static unsigned long int hash(const std::pair<long long int, long int>& r)
   {
      return hash_relation(r.first, r.second);
   };
   std::size_t operator()(const std::pair<long long int, long int>& r) const
   {
      return hash_relation(r.first, r.second);
   }
};

void usage()
{
   std::cerr << "usage: filter [-r[elations] relfile]" << std::endl;
   std::cerr << "              [-m[ergelevel] level]" << std::endl;
   std::cerr << "              [-f[iltmin] filtmin]" << std::endl;
   std::cerr << "              [-e[xcessmin] excessmin]" << std::endl;
   std::cerr << "              [-[mp|maxpass] maxpass]" << std::endl;
   std::cerr << "              [-[ro|relations_output] relations_output_file]" << std::endl;
   std::cerr << "              [-[rso|relation_sets_output] relation_sets_output_file]" << std::endl;
   std::cerr << "              [-[rspo|relation_sets_prime_output] relation_sets_prime_output_file]" << std::endl;
   std::cerr << "              [-[rsi|relation_sets_input] relation_sets_input_file]" << std::endl;
   std::cerr << "              [-[rspi|relation_sets_prime_input] relation_sets_prime_input_file]" << std::endl;
   std::cerr << "              [-merge[_only]]" << std::endl;
   exit(-1);
}

void init(int argc, char* argv[])
{
   try
   {
      int arg = 1;
      while (arg < argc)
      {
         if (strcmp(argv[arg], "-relations") == 0 ||
               strcmp(argv[arg], "-r") == 0)
         {
            arg++;
            RelationsFile = argv[arg];
         }
         else if (strcmp(argv[arg], "-mergelevel") == 0 ||
                  strcmp(argv[arg], "-m") == 0)
         {
            arg++;
            MergeLevel = atoi(argv[arg]);
         }
         else if (strcmp(argv[arg], "-minmergelevel") == 0 ||
                  strcmp(argv[arg], "-mm") == 0)
         {
            arg++;
            MinMergeLevel = atoi(argv[arg]);
         }
         else if (strcmp(argv[arg], "-filtmin") == 0 ||
                  strcmp(argv[arg], "-f") == 0)
         {
            arg++;
            FiltMin = atoi(argv[arg]);
         }
         else if (strcmp(argv[arg], "-excessmin") == 0 ||
                  strcmp(argv[arg], "-e") == 0)
         {
            arg++;
            ExcessMin = atoi(argv[arg]);
         }
         else if (strcmp(argv[arg], "-maxpass") == 0 ||
                  strcmp(argv[arg], "-mp") == 0)
         {
            arg++;
            MaxPass = atoi(argv[arg]);
         }
         else if (strcmp(argv[arg], "-merge_only") == 0 ||
                  strcmp(argv[arg], "-merge") == 0)
         {
            MergeOnly = true;
         }
         else if (strcmp(argv[arg], "-maxdiscard") == 0 ||
                  strcmp(argv[arg], "-md") == 0)
         {
            arg++;
            MaxDiscard = atoi(argv[arg]);
         }
         else if (strcmp(argv[arg], "-ro") == 0 ||
                  strcmp(argv[arg], "-relations_output_file") == 0)
         {
            arg++;
            Reloutfilename = argv[arg];
            Reloutfile = new std::fstream(argv[arg], std::ios::out);
         }
         else if (strcmp(argv[arg], "-rso") == 0 ||
                  strcmp(argv[arg], "-relation_sets_output_file") == 0)
         {
            arg++;
            Relsetoutfile = new std::fstream(argv[arg], std::ios::out);
         }
         else if (strcmp(argv[arg], "-rspo") == 0 ||
                  strcmp(argv[arg], "-relation_sets_prime_output_file") == 0)
         {
            arg++;
            Relsetprimeoutfile = new std::fstream(argv[arg], std::ios::out);
         }
         else if (strcmp(argv[arg], "-rsi") == 0 ||
                  strcmp(argv[arg], "-relation_sets_input_file") == 0)
         {
            arg++;
            Relsetinmmfile = new MemoryMappedFile(argv[arg]);
            std::string str;
            getline(*Relsetinmmfile, str);
            RelationSetCount = ::atol(str.c_str());
            delete Relsetinmmfile;
            Relsetinmmfile = new MemoryMappedFile(argv[arg]);
         }
         else if (strcmp(argv[arg], "-rspi") == 0 ||
                  strcmp(argv[arg], "-relation_sets_primes_input_file") == 0)
         {
            arg++;
            Relsetprimeinmmfile = new MemoryMappedFile(argv[arg]);
         }
         else if (strcmp(argv[arg], "-help") == 0 ||
                  strcmp(argv[arg], "-h") == 0)
         {
            usage();
         }
         arg++;
      }
      if (!Relsetoutfile) Relsetoutfile = &std::cout;
      if (FiltMin > 0)
      {
         size_t excessMinEst = estimateRequiredExcess(FiltMin);
         if (excessMinEst > ExcessMin)
         {
            std::cerr << "ExcessMin = " << ExcessMin << " is too small, resetting to " << excessMinEst << std::endl;
            ExcessMin = excessMinEst;
         }
      }
      if (!MergeOnly && (Relsetinmmfile || Relsetprimeinmmfile))
      {
         std::cerr << "-rsi and -rspi can only be specified when merge_only is specified" << std::endl;
         usage();
      }

      if (MergeOnly)
      {
         MemoryReserve = new char [ MemoryReserveSize ];
         Old_newhandler = std::set_new_handler(newhandler);
      }
   }
   catch (std::exception& e)
   {
      std::cerr << "init() failed : " << e.what() << std::endl;
      throw;
   }
   catch (...)
   {
      std::cerr << "init() failed : unknown exception" << std::endl;
      throw;
   }
}

void extract_primes(char* str, Relation* r, prime_map_type& prime_map)
{
   // str is a set of strings separated by spaces
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
      long int primeNum = atol(c);
      if (primeNum > FiltMin)
      {
         int index = 0;
         Prime prime(c);
         if (prime_map.find(prime) != prime_map.end())
         {
            index = prime_map[prime];
         }
         else
         {
            index = frequencyTable->next_prime();
            prime_map[prime] = index;
            frequencyTable->add_new_prime(index);
         }
         if (r->prime_count_ > 0 && relationTable->get_prime(r->primes_index_, r->prime_count_ - 1) == index)
         {
            r->prime_count_--;
            frequencyTable->decrement_prime(index);
         }
         else
         {
            relationTable->set_prime(r->primes_index_, r->prime_count_, index);
            r->prime_count_++;
            frequencyTable->increment_prime(index);
         }
      }

      c = d;
   }
}

void parse(const std::string& str, char*& alg_str, char*& rat_str)
{
   long long int a;
   long long int b;
   Convert::parse_FBGNFS(str, a, b, alg_str, rat_str);
}

void remove_singletons()
{
   bool done = false;
   while (!done)
   {
      int relations_removed = 0;
      // for each relation in relation table
      int k = 0;
      Relation* iter = relationTable->begin();
      for (size_t j = 0; j < relationTable->size(); ++j, ++iter)
      {
         if (iter->is_clear()) continue;
         // if relations has a prime with freq == 1 in frequency table
         bool singleton = false;
         for (int i = 0; !singleton && i < iter->prime_count_; ++i)
         {
            if (frequencyTable->frequency(relationTable->get_prime(iter->primes_index_, i)) == 1)
            {
               singleton = true;
            }
         }
         if (singleton)
         {
            // for each prime in relation
            for (int i = 0; i < iter->prime_count_; ++i)
            {
               // decrement count in frequency table
               frequencyTable->decrement_prime(relationTable->get_prime(iter->primes_index_, i));
            }
            // remove relation from relation table
            iter->clear();
            ++relations_removed;
         }
         else
         {
            (*relationTable)[k] = (*relationTable)[j];
            k++;
         }
      }
      relationTable->set_size(k);
      if (relations_removed == 0) done = true;
      else
      {
         std::cerr << relations_removed << " singleton relations removed" << std::endl;
      }
   }
   std::cerr << relationTable->size() << " relations remaining" << std::endl;
   std::cerr << frequencyTable->prime_count() << " primes remaining" << std::endl;
   Excess = relationTable->size() - frequencyTable->prime_count();
   std::cerr << Excess << " excess relations over primes" << std::endl;
}

struct Hasher
{
   static unsigned long int hash(int i)
   {
      return i;
   };
   std::size_t operator()(int i) const
   {
       return i;
   }
};
void do_clique_processing()
{
   Graph<Relation*> cliqueGraph;
   int cliqueCount;

   {
      std::unordered_map<int, Relation*, Hasher> twoMergeTable;
      for (auto& rel: *relationTable)
      {
         if (rel.is_clear()) continue;
         for (int i = 0; i < rel.prime_count_; ++i)
         {
            // this relation can take part in a 2-merge
            if (frequencyTable->frequency(relationTable->get_prime(rel.primes_index_, i)) == 2)
            {
               int index = relationTable->get_prime(rel.primes_index_, i);
               if (twoMergeTable.find(index) == twoMergeTable.end())
               {
                  twoMergeTable[index] = &rel;
               }
               else
               {
                  Relation* rel1 = twoMergeTable[index];
                  Relation* rel2 = &rel;
                  cliqueGraph.connect(rel1, rel2);
               }
            }
         }
      }
   }

   cliqueGraph.connected_components();
   cliqueCount = cliqueGraph.get_component_count();

   std::cerr << "cliqueGraph has " << cliqueCount << " connected components" << std::endl;

   // now weigh each clique, to see which should be discarded.

   const int max_power_of_half = 16;
   static double powers_of_half[max_power_of_half]
   =
   {
      1.0,
      0.5,
      0.5*0.5,
      0.5*0.5*0.5,
      0.5*0.5*0.5*0.5,
      0.5*0.5*0.5*0.5*0.5,
      0.5*0.5*0.5*0.5*0.5*0.5,
      0.5*0.5*0.5*0.5*0.5*0.5*0.5,
      0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5,
      0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5,
      0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5,
      0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5,
      0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5,
      0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5,
      0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5,
      0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5
   };

   clique_queue_type cliqueQueue;
   for (int clique = 1; clique <= cliqueCount; ++clique)
   {
      // add 1 for each relation, to account for the small primes
      double weight = cliqueGraph.component_length(clique);
      for (auto ci = cliqueGraph.begin(clique);
            ci != cliqueGraph.end(clique);
            ++ci)
      {
         Relation* r = *ci;

         // ??? something here for free relations ???
         // ??? NOT YET IMPLEMENTED ???

         // add the contribution of each prime, according to frequency
         for (int i = 0; i < r->prime_count_; ++i)
         {
            int f = frequencyTable->frequency(relationTable->get_prime(r->primes_index_, i));
            if (f > 1 && f < max_power_of_half + 2) weight += powers_of_half[f - 2];
         }
      }
      cliqueQueue.add(clique, weight);
   }

   cliqueQueue.sort();
   // now remove the "heaviest" cliques

   size_t relations_removed = 0;
   bool done = false;
   for (auto cqi = cliqueQueue.begin();
         !done && cqi != cliqueQueue.end();
         ++cqi)
   {
      int clique = *cqi;
      for (auto ci = cliqueGraph.begin(clique);
            ci != cliqueGraph.end(clique);
            ++ci)
      {
         Relation* rpp = *ci;
         // To remove relation, remove entry from RelationTable
         // and decrement frequencies in frequencyTable
         // for each prime in relation
         for (int i = 0; i < rpp->prime_count_; ++i)
         {
            // decrement count in frequency table
            frequencyTable->decrement_prime(relationTable->get_prime(rpp->primes_index_, i));
         }
         // remove relation from relation table
         rpp->clear();
         ++relations_removed;
      }

      if (relationTable->size() - relations_removed - frequencyTable->prime_count() <= ExcessMin) done = true;
   }
   std::cerr << relationTable->size() - relations_removed << " relations remaining" << std::endl;
   std::cerr << frequencyTable->prime_count() << " primes remaining" << std::endl;
   std::cerr << relationTable->size() - relations_removed - frequencyTable->prime_count() << " excess relations over primes" << std::endl;

   // Tidy up RelationTable
   size_t k = 0;
   for (Relation* iter = relationTable->begin();
         iter != relationTable->end();
         ++iter)
   {
      if (!iter->is_clear())
      {
         (*relationTable)[k] = *iter;
         k++;
      }
   }
   if (relationTable->size() - relations_removed != k)
   {
      std::cerr << "Problem with RelationTable : relationTable->size() was " << relationTable->size() << ", relations_removed = " << relations_removed << ", k = " << k << std::endl;
      throw "RelationTable corruption";
   }
   relationTable->set_size(k);
   Excess = relationTable->size() - frequencyTable->prime_count();
}

void remove_duplicates()
{
   std::unordered_map<std::pair<long long int, long int>, char, RelationHasher> hashTable;
   std::string str;
   long int relations = 0;
   long int unique_relations = 0;
   std::cerr << "Reading relations and removing duplicates ..." << std::endl;
   MemoryMappedFile relmmfile(RelationsFile.c_str());
   while (getline(relmmfile, str))
   {
      if (relations % 10000L == 0L)
      {
         std::cerr << "row <" << relations << ">" << std::endl;
      }
      std::string str1(str);
      // First split line by colons
      long long int a = 0;
      long long int b = 0;
      Convert::parse_FBGNFS_relation(str, a, b);

      std::pair<long long int, long int> h(a,b);
      hashTable[h]++;
      if (hashTable[h] == 1 && Reloutfile)
      {
         *Reloutfile << str1 << std::endl;
         unique_relations++;
      }
      relations++;
   }
   std::cerr << relations << " relations read from " << RelationsFile.c_str() << std::endl;
   std::cerr << unique_relations << " unique relations written to " << Reloutfilename << std::endl;
}

void calculate_relation_table_size(MemoryMappedFile& relmmfile,
                                   unsigned long long int relations_file_size,
                                   long int& max_relations,
                                   long int& max_primes,
                                   long int& max_unique_primes)
{
   /*
     Figures from relation file:
     
     #Relations	#Primes		File size
     48116	1011764		6980764
     50839	1069066		7375622
     55248	1161808		8015635
     58470	1229642		8482977
     64002	1345846		9286033
     69039	1451802		10016126
     87746	1844443		12727218
     108071	2271625		15677206
     115016	2417469		16684197
     7162081	7027158		1462925822

     Roughly
     File size / 145 = #Relations
     File size / 7 = #Primes
    */
   //const double bytes_per_relation = 204.0; //145.0;
   //const double bytes_per_relation = 145.0;
   //const double bytes_per_prime = 20; //6.9;
   //const double bytes_per_prime = 6.9;
   //const double bytes_per_unique_prime = 70.0;
   //const double bytes_per_unique_prime = 65.0;
   const double bytes_per_unique_prime = 50.0;
   std::cerr << "Size of relations file is " << relations_file_size << " bytes" << std::endl;
   // Estimate number of bytes per relation, and hence the number of relations to expect,
   // by reading the beginning of the relation file
   std::string str;
   size_t relation_count = 0;
   size_t bytes = 0;
   size_t prime_count = 0;
   while (getline(relmmfile, str) && relation_count < 10000)
   {
      bytes += str.length();
      const char* p = str.c_str();
      while (p && *p)
      {
         if (*p == ' ')
         {
            ++prime_count;
         }
         ++p;
      }
      prime_count -= 4; // allow for a, b and : chars
      ++relation_count;
   }
   relmmfile.reset();
   double bytes_per_relation = (double)bytes / (double)relation_count;
   double bytes_per_prime = (double)bytes / (double)prime_count;

   max_relations = (double)relations_file_size / bytes_per_relation;
   max_relations += max_relations / 10;
   if (max_relations > 28000000L)
   {
      max_relations = 28000000L;
   }
   max_primes = (double)relations_file_size / bytes_per_prime;
   if (max_primes > 400000000L)
   {
      max_primes = 400000000L;
   }
   max_unique_primes = (double)relations_file_size / bytes_per_unique_prime;

   // Round up to nearest multiple of 100000
   const long int n = 100000L;
   max_relations /= n;
   ++max_relations;
   max_relations *= n;
   max_primes /= n;
   ++max_primes;
   max_primes *= n;
   max_unique_primes /= n;
   ++max_unique_primes;
   max_unique_primes *= n;
   std::cerr << "Relation table size : max_relations = " << max_relations << ", max_primes = " << max_primes << ", max_unique_primes = " << max_unique_primes << std::endl;
   std::cerr << "Bytes per relation = " << bytes_per_relation << ", bytes per prime = " << bytes_per_prime << std::endl;
}

void read_relations()
{
   // This function assumes no duplicates in the relations
   std::cerr << "Reading relations from " << RelationsFile << " ..." << std::endl;
   MemoryMappedFile relmmfile(RelationsFile.c_str());

   unsigned long long int relations_file_size = relmmfile.size();
   long int max_relations;
   long int max_primes;
   long int max_unique_primes;
   calculate_relation_table_size(relmmfile, relations_file_size, max_relations, max_primes, max_unique_primes);
   relationTable = new RelationTable(max_relations, max_primes);
   frequencyTable = new PrimeFrequencyTable(max_unique_primes);

   prime_map_type prime_map;

   long int relations = 0;
   std::string str;
   while (getline(relmmfile, str))
   {
      if (relations % 10000L == 0L)
      {
         std::cerr << "row <" << relations << ">" << std::endl;
      }
      std::string str1(str);
      // First split line by colons
      char* alg_str;
      char* rat_str;
      parse(str, alg_str, rat_str);

      // add relation to relation table
      // and for each (large) prime in relation
      // increment count in frequency table
      // Split the algebraic primes
      Relation r(relations);
      r.primes_index_ = relationTable->next_prime_index();
      extract_primes(alg_str, &r, prime_map);
      // Split the rational primes
      extract_primes(rat_str, &r, prime_map);
      (*relationTable)[relations] = r;
      relationTable->increment_next_prime(r.prime_count_);
      relations++;
   }
   std::cerr << relations << " relations read" << std::endl;
   std::cerr << frequencyTable->next_prime() << " unique primes read" << std::endl;
   relationTable->set_size(relations);
   relationTable->shrink_to_fit();
}

void write_relations()
{
   if (!Reloutfile)
   {
      std::cerr << "No relation output file defined - not writing" << std::endl;
      return;
   }
   std::string str;
   Relation* iter = relationTable->begin();
   long int row = 0L;
   MemoryMappedFile relmmfile(RelationsFile.c_str());
   while (iter != relationTable->end() && getline(relmmfile, str))
   {
      if (row % 10000L == 0L)
      {
         std::cerr << "row <" << row << ">" << std::endl;
      }
      if (iter->id_ == row)
      {
         *Reloutfile << str << std::endl;
         ++iter;
      }
      row++;
   }
}

//const int m_max = 7;
int m = 1;
bool inc = false;
void increment_m()
{
   if (inc)
   {
      ++m;
      inc = false;
   }
   else inc = true;
}

bool merge_ok(RelationManager& RelationSets, const std::vector<long int>& relation_sets, long int prime)
{
   // find lightest relation
   // "Suppose the lightest candidate relation-set has j primes (was relations), ...
   //  Let c be the number of relation-sets with exactly this minimum number j of primes (was relations)."
   int k = relation_sets.size();
   int j = 1000;
   int c = 0;
   for (auto& rs: relation_sets)
   {
      int w = RelationSets.prime_weight(rs);
      if (w < j)
      {
         j = w;
         c = 0;
      }
      if (w == j) ++c;
   }
   if (!inc) c = 0;
   if ((k - 2)*j <= 70 * (m - (c - 1) / 2)) return true;
   return false;
}

void merge_relation_sets(RelationManager& RelationSets, long int prime, const std::vector<long int>& relation_sets)
{
   LOG_DEBUG("merge_relation_sets : prime = " << prime << ", relation_sets.size() = " << relation_sets.size());
   if (relation_sets.size() == 2)
   {
      RelationSets.merge(relation_sets[0], relation_sets[1], prime);
      ++RelationsMergedInThisPass;
   }
   else
   {
      // Do we want to allow this merge
      if (!merge_ok(RelationSets, relation_sets, prime)) 
      {
          return;
      }
      LOG_DEBUG("merge_relation_sets : prime = " << prime << ", merge_ok");

      // construct a fully connected graph from the relation sets
      // and use minimum spanning tree algorithm
      Graph<long int, RelationSetWeight> g;
      for (auto it1 = relation_sets.begin();
            it1 != relation_sets.end();
            ++it1)
      {
         auto it2 = it1;
         ++it2;
         for (;it2 != relation_sets.end(); ++it2)
         {
            g.connect(*it1, *it2);
         }
      }
      Graph<long int> tree;
      g.minimum_spanning_tree(tree);

      // connections from tree give the pairs of relation sets to merge

      for (const auto& node: tree)
      {
         long int rs1 = node.first;
         long int rs2 = node.second;
         if (rs1 < rs2)
         {
            RelationSets.merge(rs1, rs2, prime);
            ++RelationsMergedInThisPass;
         }
      }
   }

   // remove the original relation sets
   for (auto& rs: relation_sets)
   {
      RelationSets.remove(rs);
   }
   if (frequencyTable->frequency(prime) != 0)
   {
      LOG_DEBUG("Problem : miscalculation of prime frequency for prime #" << prime << ", freq = " << frequencyTable->frequency(prime));
      throw "frequencyTable corruption";
   }
}

void do_merge_processing()
{
   LOG_DEBUG("Entering do_merge_processing() ...");
   LOG_DEBUG("Creating RelationSets ... ");
   RelationManager* RelationSets = 0;
   if (!frequencyTable)
   {
      frequencyTable = new PrimeFrequencyTable(1);
   }
   if (Relsetinmmfile)
   {
      if (Relsetprimeinmmfile)
      {
         // In this case frequencyTable is reset to size of relation set prime map (and hence cleared)
         RelationSets = new RelationManager(RelationSetCount, *Relsetinmmfile, *Relsetprimeinmmfile, *frequencyTable, ExcessMin);
      }
      else
      {
         // In this case frequencyTable is cleared and we have read relations, so relationTable is defined
         RelationSets = new RelationManager(*Relsetinmmfile, *relationTable, *frequencyTable, ExcessMin);
      }
   }
   else
   {
      if (!Relsetprimeinmmfile)
      {
         // In this case frequencyTable is simply checked for consistency, but otherwise unchanged
         RelationSets = new RelationManager(*relationTable, *frequencyTable, ExcessMin);
      }
      else
      {
         throw "Invalid combination of program arguments";
      }
   }
   delete Relsetprimeinmmfile;

   LOG_DEBUG("Created RelationSets");
   RelationSets->stats();
   for (long int pass = 0; pass < MaxPass; ++pass)
   {
      RelationsMergedInThisPass = 0;
      RelationSets->remove_heavy_relation_sets();
      long int merge_level = MergeLevel - MaxPass + pass;
      if (merge_level < static_cast<long int>(MinMergeLevel)) merge_level = MinMergeLevel;
      LOG_DEBUG("pass = " << pass << ", merge level = " << merge_level);

      RelationSets->build_prime_relation_set_map(merge_level);
      RelationSets->stats();
      for (size_t prime = 0; RelationsMergedInThisPass < CutoffForRelationsMergedInThisPass && !NewHandlerCalled && prime < frequencyTable->capacity(); ++prime)
      {
         if (frequencyTable->frequency(prime) > 0)
         {
            std::vector<long int> relation_sets;
            size_t count = RelationSets->sets_including_prime(prime, relation_sets);
            // LOG_DEBUG("# sets including prime " << prime << " is " << count);
            if (count > 1 && count <= MergeLevel)
            {
               if (count != frequencyTable->frequency(prime))
               {
                  LOG_ERROR("Problem: for prime #" << prime << " RelationSets says count = " << count << " but frequencyTable says " << frequencyTable->frequency(prime));
                  throw "frequencyTable out of synch with RelationSets";
               }
	           LOG_DEBUG("# sets including prime " << prime << " is " << count);

               merge_relation_sets(*RelationSets, prime, relation_sets);
            }
         }
         else
         {
             // LOG_DEBUG("frequencyTable->frequency(" << prime << ") <= 0");
         }
      }
      LOG_DEBUG("NewHandlerCalled = " << NewHandlerCalled);
      LOG_DEBUG("RelationsMergedInThisPass = " << RelationsMergedInThisPass);
      increment_m();
      RelationSets->stats();
   }
   // print out relation sets

   if (Relsetprimeoutfile) RelationSets->write_primes_to_file(*Relsetprimeoutfile);
   delete Relsetprimeoutfile;

   if (Relsetoutfile) RelationSets->write_to_file(*Relsetoutfile);
   delete Relsetoutfile;
   delete Relsetinmmfile;
   delete RelationSets;
}

}

int main(int argc, char* argv[])
{
   try
   {
      LogManager::instance().start_logging("filter.log");
      LogManager::instance().set_level(Logger::debug);
      init(argc, argv);

      if (MergeLevel == 0)
      {
         remove_duplicates();
      }
      else
      {
         if (!Relsetprimeinmmfile)
         {
            read_relations();
            if (!MergeOnly)
            {
               remove_singletons();
               long int prev_excess = Excess;
               bool done = false;

               while (MergeLevel > 1 && !done)
               {
                  do_clique_processing();
                  remove_singletons();
                  if (Excess <= static_cast<long int>(ExcessMin) || Excess >= prev_excess - 10)
                  {
                      done = true;
                  }
                  else
                  {
                      prev_excess = Excess;
                  } 
               }
               write_relations();
            }
            //write_relations();
         }

         if (MergeOnly)
         {
         }

         if (MergeLevel > 1 && MaxPass > 0)
         {
            do_merge_processing();
         }
      }
   }
   catch (const char* error)
   {
      std::cerr << "Fatal error : " << error << std::endl;
   }
   catch (std::string& error)
   {
      std::cerr << "Fatal error : " << error << std::endl;
   }
   catch (std::exception& e)
   {
      std::cerr << "Fatal error : " << e.what() << std::endl;
   }
   catch (...)
   {
      std::cerr << "Fatal error : unknown exception" << std::endl;
   }
   return 0;
}
