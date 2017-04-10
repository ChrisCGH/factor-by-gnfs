#include "SparseMatrix.h"
#include "RelationManager.h"
#include "MemoryMappedFile.h"
#include "VeryLong.h"
#include "Polynomial.h"
#include "SieveConfig.h"
#include "QuadraticCharacters.h"
#include "convert.h"
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <sstream>

namespace
{
std::string RelationsFile;
MemoryMappedFile* Relmmfile = 0;
MemoryMappedFile* Relsetmmfile = 0;
MemoryMappedFile* Relsetprimesmmfile = 0;
std::ostream* Matrixfile = 0;
std::string Matrixfilename;
std::string ConfigFile;
int FiltMin = 0;

void usage(const std::string& msg = "")
{
   if (msg != "") std::cerr << msg << std::endl;
   std::cerr << "usage: buildMatrix [-r[elations] relfile]" << std::endl;
   std::cerr << "                   [-rs|-relationsets] relsetfile" << std::endl;
   std::cerr << "                   [-m[atrix] matrix]" << std::endl;
   std::cerr << "                   [-c[onfig] config]" << std::endl;
   std::cerr << "                   [-f[iltmin] filtmin]" << std::endl;
   std::cerr << "                   [-h[elp]]" << std::endl;
   std::exit(-1);
}

void init(int argc, char* argv[])
{
   int arg = 1;
   while (arg < argc)
   {
      if (strcmp(argv[arg], "-relations") == 0 ||
            strcmp(argv[arg], "-r") == 0)
      {
         ++arg;
         RelationsFile = argv[arg];
         Relmmfile = new MemoryMappedFile(argv[arg]);
      }
      else if (strcmp(argv[arg], "-relationsets") == 0 ||
               strcmp(argv[arg], "-rs") == 0)
      {
         ++arg;
         Relsetmmfile = new MemoryMappedFile(argv[arg]);
      }
      else if (strcmp(argv[arg], "-relationsetprimes") == 0 ||
               strcmp(argv[arg], "-rsp") == 0)
      {
         ++arg;
         Relsetprimesmmfile = new MemoryMappedFile(argv[arg]);
      }
      else if (strcmp(argv[arg], "-matrix") == 0 ||
               strcmp(argv[arg], "-m") == 0)
      {
         ++arg;
         Matrixfilename = argv[arg];
         Matrixfile = new std::fstream(Matrixfilename.c_str(), std::ios::out);
      }
      else if (strcmp(argv[arg], "-config") == 0 ||
               strcmp(argv[arg], "-c") == 0)
      {
         ++arg;
         ConfigFile = argv[arg];
      }
      else if (strcmp(argv[arg], "-filtmin") == 0 ||
               strcmp(argv[arg], "-f") == 0)
      {
         ++arg;
         FiltMin = std::atoi(argv[arg]);
      }
      else if (strcmp(argv[arg], "-help") == 0 ||
               strcmp(argv[arg], "-h") == 0)
      {
         usage();
      }
      ++arg;
   }
   if (!Relmmfile) usage("Please specify relation file with -r option");
   if (!Matrixfile) Matrixfile = &std::cout;
   if (!Relsetmmfile) usage("Please specify relation sets file with -rs option");
   if (!Relsetprimesmmfile) usage("Please specify relation set primes file with -rsp option");
   if (ConfigFile == "") usage("Please specify config file with -c option");
}

void reopen()
{
   delete Relmmfile;
   Relmmfile = new MemoryMappedFile(RelationsFile.c_str());
}

void cleanup()
{
   delete Relmmfile;
   if (Matrixfile != &std::cout) delete Matrixfile;
   delete Relsetmmfile;
   delete Relsetprimesmmfile;
}

typedef std::pair<long int, std::string> relPrimePair;

char* buf = 0;
std::string::size_type buflen = 0;
void copy_str(const std::string& str)
{
   if (str.size() > buflen)
   {
      delete [] buf;
      buflen = str.size();
      buf = new char [ buflen + 1 ];
   }
   strcpy(buf, str.c_str());
}

void split1(const std::string& str, std::set<std::string>& base)
{
   // str is a set of strings separated by spaces
   if (str.empty()) return;
   copy_str(str);
   char* s = strtok(buf, " ");
   if (std::atoi(s) <= FiltMin)
   {
      std::string prime(s);
      base.insert(prime);
      while ((s = strtok(0, " ")))
      {
         if (std::atoi(s) <= FiltMin)
         {
            prime = s;
            base.insert(prime);
         }
      }
   }
}

void split2(const std::string& str, std::vector<std::string>& primes)
{
   primes.clear();
   if (str.empty()) return;
   copy_str(str);
   char* s = strtok(buf, " ");
   if (std::atoi(s) <= FiltMin)
   {
      std::string prime(s);
      primes.push_back(prime);
      while ((s = strtok(0, " ")))
      {
         if (std::atoi(s) <= FiltMin)
         {
            prime = s;
            primes.push_back(prime);
         }
      }
   }
   std::sort(primes.begin(), primes.end());
}

typedef SparseMatrix SPARSEMATRIX;
#define SMALL_PRIMES_AS_TRANSPOSE 1

void add_primes(std::vector<std::string>& primes,
                const std::map<std::string, long int>& base_index,
                long int row_index,
                long int colstart,
                SPARSEMATRIX& smat)
{
   std::string prev_prime = "";
   int keep = 0;
   for (size_t i = 0; i < primes.size(); i++)
   {
      std::string prime = primes[i];
      if (prime != prev_prime)
      {
         if (keep)
         {
            std::map<std::string, long int>::const_iterator found = base_index.find(prev_prime);
            if (found != base_index.end())
            {
               long int index = found->second;
               long int col = colstart + index;
#ifdef SMALL_PRIMES_AS_TRANSPOSE
               smat.xor(row_index, col);
#else
               smat.xor(col, row_index);
#endif
            }
            else
            {
               throw "Problem: cannot find prime in base_index";
            }
         }
         prev_prime = prime;
         keep = 1;
      }
      else
      {
         keep = 1 - keep;
      }
   }
   if (keep)
   {
      std::map<std::string, long int>::const_iterator found = base_index.find(prev_prime);
      if (found != base_index.end())
      {
         long int index = found->second;
         long int col = colstart + index;
#ifdef SMALL_PRIMES_AS_TRANSPOSE
         smat.xor(row_index, col);
#else
         smat.xor(col, row_index);
#endif
      }
      else
      {
         throw "Problem: cannot find prime in base_index";
      }
   }
}

void add_qcs(std::vector<char>& qcs,
             long int relation_set_index,
             int num_qcs,
             //SPARSEMATRIX& smat)
             BitMatrix64& qc_bm)
{
   for (int i = 0; i < num_qcs; ++i)
   {
      if (qcs[i] == 1)
      {
         BitOperations64::toggleBit(i, qc_bm.row_[relation_set_index]);
      }
   }
}

bool parse(const std::string& str, std::string& alg_str, std::string& rat_str)
{
    char* as;
    char* rs;
    long long int a;
    long long int b;
    Convert::parse_FBGNFS(str, a, b, as, rs);
    alg_str = as;
    rat_str = rs;
    return true;
}

bool parse(const std::string& str, VeryLong& a, VeryLong& b, std::string& alg_str, std::string& rat_str)
{
    long long int aa;
    long long int bb;
    char* as;
    char* rs;
    Convert::parse_FBGNFS(str, aa, bb, as, rs);
    a = aa;
    b = bb;
    alg_str = as;
    rat_str = rs;
    return true;
}

void firstPass(std::map<std::string, long int>& alg_base_index,
               std::map<std::string, long int>& rat_base_index,
               int num_qcs,
               long int& row_length,
               long int& rat_col_start)
{
   /*
    * Read relations to get the small primes (< Filtmin) which are not included in the
   * relation set prime map, since we missed out these primes in the filter stage
    */
   std::cerr << "First pass of relations ..." << std::endl;
   std::set<std::string> alg_base;
   std::set<std::string> rat_base;
   std::string str;
   long int relations = 0;
   while (::getline(*Relmmfile, str))
   {
      if (relations % 10000L == 0L)
      {
         std::cerr << "row <" << relations << ">" << std::endl;
      }
      // First split line by colons
      std::string alg_str;
      std::string rat_str;
      parse(str, alg_str, rat_str);

      // Split the algebraic primes
      split1(alg_str, alg_base);
      // Split the rational primes
      split1(rat_str, rat_base);
      ++relations;
   }

   size_t i = 0;
   for (auto& i1: alg_base)
   {
      alg_base_index[i1] = i;
      ++i;
   }
   int alg_base_size = i;

   i = 0;
   for (auto& i1: rat_base)
   {
      rat_base_index[i1] = i;
      ++i;
   }
   int rat_base_size = i;

   row_length = alg_base_size + rat_base_size + num_qcs;
   rat_col_start = alg_base_size;
   std::cerr << "alg_base.size() = " << alg_base_size << std::endl;
   std::cerr << "rat_base.size() = " << rat_base_size << std::endl;
   std::cerr << "num_qcs = " << num_qcs << std::endl;

   std::cerr << "row length = " << row_length << std::endl;
   std::cerr << "relations = " << relations << std::endl;
}

size_t writeMainMatrix(long int row_length)
{
   // Read relation sets primes file, transpose it and write out
   SparseMatrix relation_sets_primes(0,0);
   std::cerr << "reading relation set prime map ... " << std::flush;
   relation_sets_primes.read(*Relsetprimesmmfile);
   std::cerr << "done" << std::endl;
   size_t relation_set_count = relation_sets_primes.rows();

   SparseMatrix primes_relation_sets;
   std::cerr << "transposing ... " << std::flush;
   transpose(relation_sets_primes, primes_relation_sets, 0, true);
   std::cerr << "done" << std::endl;
   relation_sets_primes.clear();
   std::cerr << "removing empty rows ... " << std::flush;
   primes_relation_sets.removeEmptyRows();
   std::cerr << "done" << std::endl;
   primes_relation_sets.set_write_row_count(false);
   std::cerr << "writing primes -> relations sets ... " << std::flush;
   *Matrixfile << primes_relation_sets.rows() + row_length << std::endl;
   *Matrixfile << primes_relation_sets;
   std::cerr << "done" << std::endl;
   primes_relation_sets.clear();
   delete Relsetprimesmmfile;
   Relsetprimesmmfile = 0;
   return relation_set_count;
}

void generateMatrix()
{
   /*
    * A typical relations line is
   a b : algebraic primes : rational primes :
   e.g.
   -982 1 : 3/2 3/3 23/7 9421/8439 16127/15145 65651/64669 : 2 1741 :
   */
   std::map<std::string, long int> alg_base_index;
   std::map<std::string, long int> rat_base_index;
   int num_qcs = 50;
   long int row_length;
   long int rat_col_start;
   firstPass(alg_base_index, rat_base_index, num_qcs, row_length, rat_col_start);

   reopen();
   /*
    * At this point we have defined:
    *   alg_base_index : map from small algebraic primes to index
    *   rat_base_index : map from small rational primes to index
    *   num_qcs        : number of quadratic characters
    *   row_length     : total number of small primes and quadratic characters, which
    *                    is the number of rows which will be added to the final matrix
    *   rat_col_start  : position within row of first rational prime column
    */

   /*
    * write out main part of the matrix which is the transposed relation set / prime map produced from
    * the filtering
    */
   size_t relation_set_count = writeMainMatrix(row_length);

   // Read the relation sets
   RelationSetManager relation_sets(*Relsetmmfile);
   relation_sets.build_relation_relation_set_map(true);

   /*
    * Read the relations file again and build the rest of the matrix, consisting of one row for each small prime
    * and one row for each quadratic character.
    */
   //reopen();
   /*
    * We use a Sparse matrix to hold the rows corresponding to the small primes and a BitMatrix to hold the
    * quadratic characters.
    * The small prime rows are quite short in general, but the quadratic character rows are long
    * because roughly half of the relation sets (the columns) will be set, so it's better to hold them (transposed) as a
    * BitMatrix.
    */
#ifdef SMALL_PRIMES_AS_TRANSPOSE
   SPARSEMATRIX smat(relation_set_count);
#else
   SPARSEMATRIX smat(row_length - num_qcs);
#endif
   //size_t i = 0;
   /*
      for (i = 0; i < row_length - num_qcs; ++i)
      {
         smat.get_row(i).reserve(4000L);
         smat.get_row(i).set_compress_cutoff(4000L);
      }
      */

   BitMatrix64 qc_bm(relation_set_count, num_qcs);

   SieveConfig config(ConfigFile);
   Polynomial<VeryLong> f1 = config.f1();
   long int L1 = config.L1();
   // As long as f2 is linear, we only need to calculate quadratic characters for f1
   QuadraticCharacters qc1(f1, L1);
   std::vector<char> qcs1;

   std::string str;
   long int relation_index = 0L;
   while (::getline(*Relmmfile, str))
   {
      if (relation_index % 10000L == 0L)
      {
         std::cerr << "relation_index <" << relation_index << ">" << std::endl;
      }
      // First split line by colons
      VeryLong a;
      VeryLong b;
      std::string alg_str;
      std::string rat_str;
      parse(str, a, b, alg_str, rat_str);

      std::vector<std::string> alg_primes;
      split2(alg_str, alg_primes);

      std::vector<std::string> rat_primes;
      split2(rat_str, rat_primes);

      qc1.generate(a, b, qcs1);

      // find relation set containing relation
      std::vector<long int> rel_sets;
      relation_sets.relation_sets_for_relation(relation_index, rel_sets);
      for (auto& rel_set: rel_sets)
      {
         long int rs = rel_set;

         add_primes(alg_primes, alg_base_index, rs, 0L, smat);
         add_primes(rat_primes, rat_base_index, rs, rat_col_start, smat);
         add_qcs(qcs1, rs, num_qcs, qc_bm);
      }
      ++relation_index;
   }

   smat.compress();

   std::cerr << "smat has " << smat.rows() << " rows" << std::endl;
   std::cerr << "smat has " << smat.cols() << " cols" << std::endl;
   rat_base_index.clear();
   alg_base_index.clear();
   relation_sets.clear();
#ifdef SMALL_PRIMES_AS_TRANSPOSE
   SparseMatrix smat1;
   transpose(smat, smat1, 0, false);
   smat.clear();
   std::cerr << "Writing smat1 which has " << smat1.rows() << " rows and " << smat1.cols() << " cols" << std::endl;
   smat1.set_write_row_count(false);
   *Matrixfile << smat1;
   smat1.clear();
#else
   smat.set_write_row_count(false);
   *Matrixfile << smat;
   smat.clear();
#endif
   std::cerr << "Writing Quadratic Characters in BitMatrix with " << qc_bm.rows() << " rows and " << qc_bm.cols() << " cols" << std::endl;
   qc_bm.printTransposeAsSparseMatrix(*Matrixfile);

}

}

int main(int argc, char* argv[])
{
   try 
   {
      init(argc, argv);
      generateMatrix();
      cleanup();
   }
   catch (const std::string& s)
   {
      std::cerr << "Exception: " << s << std::endl;
   }
   return 0;
}
