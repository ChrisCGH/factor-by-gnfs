#ifndef SIEVE_H
#define SIEVE_H
#include "VeryLong.h"
#include "Polynomial.h"
#include "FactorBase.h"
#include "SieveConfig.h"
#include <utility>
#include <string>

class Siever
{
   public:
      Siever(const std::string& config_file = "sieve.cfg");
      ~Siever();
      void sieve(long int min_A, long int max_A, long int b);
      void checkRelations(const char* relfile);

   private:
#ifdef SHORT_SIEVE
      typedef unsigned short SIEVE_TYPE;
      typedef int BIGGER_THAN_SIEVE_TYPE;
      enum
      {
         SIEVE_MAX_VALUE = 65535
      };
#elsif UNSIGNED_SIEVE
      typedef unsigned char SIEVE_TYPE;
      typedef short BIGGER_THAN_SIEVE_TYPE;
      enum
      {
         SIEVE_MAX_VALUE = 255
      };
#else
      typedef signed char SIEVE_TYPE;
      typedef short BIGGER_THAN_SIEVE_TYPE;
      enum
      {
         SIEVE_MAX_VALUE = 127
      };
#endif
      Polynomial<VeryLong> f1_;
      Polynomial<double> f1d_;
      Polynomial<VeryLong> f2_;
      Polynomial<double> f2d_;
      VeryLong m_;
      long int B1_;
      long int L1_;
      long int LP1_;
      VeryLong L_LP_1_;
      long int B2_;
      long int L2_;
      long int LP2_;
      VeryLong L_LP_2_;
      long int SIEVE_BOUND_ADJUSTMENT1_;
      long int SIEVE_BOUND_ADJUSTMENT2_;
      long int SMALL_PRIME_BOUND_;
      bool RESIEVE_;
      std::string relation_file_;
      std::fstream* relfile_;
      FactorBase* alg_factor_base_;
      FactorBase* rat_factor_base_;

      long int chunk_size_;
      long int cache_chunks_;
      long int sieve_array_size_;
      bool dont_skip_[1024];
      long int first_a_[1024];
      SIEVE_TYPE* sieve_array_;

      void sieve(const Polynomial<double>& f, const VeryLong& c_d,
                 long int min_A, long int max_A, long int b, bool is_first_b, long int L, long int LP,
                 FactorBase* factor_base, bool first_pass);

      int check_interval(const Polynomial<double>& f, long int L, long int LP,
                         long int a1, long int a2, long int b,
                         long int min_A, SIEVE_TYPE projective_bump, bool first_pass, FactorBase* factor_base);

      int find_cutoff(const Polynomial<VeryLong>& f, long int L, long int LP,
                      long int a1, long int a2, long int b);

      int check_interval(int cutoff, long int a1, long int a2, long int b,
                         long int min_A, bool first_pass, FactorBase* factor_base);

      static int checkForRelation(const VeryLong& a, const VeryLong& b,
                                  const Polynomial<VeryLong>& f,
                                  long int B, long int L, long int LP, const VeryLong& L_LP,
                                  std::vector<VeryLong>& factors);

      static int checkForRelation(const VeryLong& a, const VeryLong& b,
                                  const Polynomial<VeryLong>& f,
                                  long int B, long int L, long int LP, const VeryLong& L_LP,
                                  std::vector<VeryLong>& factors, const std::vector<long int>& primes);

      static int find_cutoff(long int min_A, long int a1, long int a2, double log_mid_value, SIEVE_TYPE* sa);

};
#endif
