#ifndef SIEVECONFIG_H
#define SIEVECONFIG_H

#include <fstream>
#include <string>
// class for configuration information for
// sieving procedure

class SieveConfig
{
   public:
      SieveConfig(const std::string& filename) : LP1_(1), LP2_(1), SKEWEDNESS_(1.0)
      {
         SIEVE_ID_ = "ctc4";
         RELATION_FILE_ = "relations.out";
         SMALL_PRIME_BOUND1_ = 0L;
         SMALL_PRIME_BOUND2_ = 0L;
         RESIEVE_ = true;
         INITIAL_CUTOFF_ = 400L;
         MIN_A_ = -600000000L;
         MAX_A_ = 600000000L;
         MIN_B_ = 1L;
         MAX_B_ = 90000L;
         DEBUG_ = false;
         FIXED_SIEVE_REGION_ = false;

         std::ifstream config_file(filename.c_str(), std::ios::in);
         if (config_file)
         {
            std::string str;
            while (getline(config_file, str))
            {
               // nothing
               if (str.empty()) continue;
               // comments
               if (str[0] == '#') continue;
               // we expect to find =
               std::string::size_type eqpos = str.find('=');
               if (std::string::npos == eqpos) continue;

               // XYZ = 123
               // 0123456789
               // size => 9, eqpos + 2 => 6
               // XYZ = 1
               // 01234567
               // size => 7, eqpos + 2 => 6
               if (eqpos + 2 >= str.size()) continue;

               std::string s = str.substr(eqpos + 2);
               if (str.find("SIEVE_ID = ") == 0)
               {
                  SIEVE_ID_ = s;
               }
               else if (str.find("N = ") == 0)
               {
                  N_ = s;
               }
               else if (str.find("f1 = ") == 0)
               {
                  if (s == "Polynomial")
                  {
                     f1_ = Polynomial<VeryLong>::read_polynomial(config_file);
                  }
                  else
                  {
                     f1_ = Polynomial<VeryLong>::read_polynomial(s.c_str());
                  }
               }
               else if (str.find("f2 = ") == 0)
               {
                  if (s == "Polynomial")
                  {
                     f2_ = Polynomial<VeryLong>::read_polynomial(config_file);
                  }
                  else
                  {
                     f2_ = Polynomial<VeryLong>::read_polynomial(s.c_str());
                  }
               }
               else if (str.find("m = ") == 0)
               {
                  m_ = s;
               }
               else if (str.find("B1 = ") == 0)
               {
                  B1_ = std::atol(s.c_str());
               }
               else if (str.find("L1 = ") == 0)
               {
                  L1_ = std::atol(s.c_str());
               }
               else if (str.find("LP1 = ") == 0)
               {
                  LP1_ = std::atol(s.c_str());
               }
               else if (str.find("B2 = ") == 0)
               {
                  B2_ = std::atol(s.c_str());
               }
               else if (str.find("L2 = ") == 0)
               {
                  L2_ = std::atol(s.c_str());
               }
               else if (str.find("LP2 = ") == 0)
               {
                  LP2_ = std::atol(s.c_str());
               }
               else if (str.find("SIEVE_BOUND_ADJUSTMENT1 = ") == 0)
               {
                  SIEVE_BOUND_ADJUSTMENT1_ = std::atol(s.c_str());
               }
               else if (str.find("SIEVE_BOUND_ADJUSTMENT2 = ") == 0)
               {
                  SIEVE_BOUND_ADJUSTMENT2_ = std::atol(s.c_str());
               }
               else if (str.find("SKEWEDNESS = ") == 0)
               {
                  SKEWEDNESS_ = std::atof(s.c_str());
               }
               else if (str.find("SMALL_PRIME_BOUND1 = ") == 0)
               {
                  SMALL_PRIME_BOUND1_ = std::atol(s.c_str());
               }
               else if (str.find("SMALL_PRIME_BOUND2 = ") == 0)
               {
                  SMALL_PRIME_BOUND2_ = std::atol(s.c_str());
               }
               else if (str.find("INITIAL_CUTOFF = ") == 0)
               {
                  INITIAL_CUTOFF_ = std::atol(s.c_str());
               }
               else if (str.find("MIN_A = ") == 0)
               {
                  MIN_A_ = std::atof(s.c_str());
               }
               else if (str.find("MAX_A = ") == 0)
               {
                  MAX_A_ = std::atof(s.c_str());
               }
               else if (str.find("MIN_B = ") == 0)
               {
                  MIN_B_ = std::atol(s.c_str());
               }
               else if (str.find("MAX_B = ") == 0)
               {
                  MAX_B_ = std::atol(s.c_str());
               }
               else if (str.find("RESIEVE = ") == 0)
               {
                  if (s == "true" || s == "TRUE" || s == "1")
                  {
                     RESIEVE_ = true;
                  }
                  else
                  {
                     RESIEVE_ = false;
                  }
               }
               else if (str.find("DEBUG = ") == 0)
               {
                  if (s == "true" || s == "TRUE" || s == "1")
                  {
                     DEBUG_ = true;
                  }
                  else
                  {
                     DEBUG_ = false;
                  }
               }
               else if (str.find("FIXED_SIEVE_REGION = ") == 0)
               {
                  if (s == "true" || s == "TRUE" || s == "1")
                  {
                     FIXED_SIEVE_REGION_ = true;
                  }
                  else
                  {
                     FIXED_SIEVE_REGION_ = false;
                  }
               }
               else if (str.find("RELATION_FILE = ") == 0)
               {
                  RELATION_FILE_ = s;
               }
            }
         }
         else
         {
            throw std::string("Failed to open config file");
         }

         if (LP1_ < 0) LP1_ = 0;
         if (LP1_ > 5) LP1_ = 5;
         if (LP2_ < 0) LP2_ = 0;
         if (LP2_ > 5) LP2_ = 5;
         // check
         VeryLong x = f1_.evaluate(m_);
         VeryLong y = f2_.evaluate(m_);
         if (x % N_ != 0L || y % N_ != 0L)
         {
            std::cerr << "Problem: f1(m) = " << x << ", f2(m) = " << y << std::endl;
         }
      }

      std::string SIEVE_ID() const
      {
         return SIEVE_ID_;
      }
      VeryLong N() const
      {
         return N_;
      }
      Polynomial<VeryLong> f1() const
      {
         return f1_;
      }
      Polynomial<VeryLong> f2() const
      {
         return f2_;
      }
      VeryLong m() const
      {
         return m_;
      }
      long int B1() const
      {
         return B1_;
      }
      long int L1() const
      {
         return L1_;
      }
      long int LP1() const
      {
         return LP1_;
      }
      long int B2() const
      {
         return B2_;
      }
      long int L2() const
      {
         return L2_;
      }
      long int LP2() const
      {
         return LP2_;
      }
      double MIN_A() const
      {
         return MIN_A_;
      }
      double MAX_A() const
      {
         return MAX_A_;
      }
      long int MIN_B() const
      {
         return MIN_B_;
      }
      long int MAX_B() const
      {
         return MAX_B_;
      }
      long int SIEVE_BOUND_ADJUSTMENT1() const
      {
         return SIEVE_BOUND_ADJUSTMENT1_;
      }
      long int SIEVE_BOUND_ADJUSTMENT2() const
      {
         return SIEVE_BOUND_ADJUSTMENT2_;
      }
      double SKEWEDNESS() const
      {
         return SKEWEDNESS_;
      }
      long int SMALL_PRIME_BOUND1() const
      {
         return SMALL_PRIME_BOUND1_;
      }
      long int SMALL_PRIME_BOUND2() const
      {
         return SMALL_PRIME_BOUND2_;
      }
      long int INITIAL_CUTOFF() const
      {
         return INITIAL_CUTOFF_;
      }
      bool RESIEVE() const
      {
         return RESIEVE_;
      }
      std::string RELATION_FILE() const
      {
         return RELATION_FILE_;
      }

      bool DEBUG() const
      {
         return DEBUG_;
      }

      bool FIXED_SIEVE_REGION() const
      {
         return FIXED_SIEVE_REGION_;
      }

      void display() const
      {
         std::cerr << "# Configuration options: " << std::endl;
         std::cerr << "SIEVE_ID = " << SIEVE_ID() << std::endl;
         std::cerr << "N = " << N() << std::endl;
         std::cerr << "f1 = " << f1() << std::endl;
         std::cerr << "f2 = " << f2() << std::endl;
         std::cerr << "m = " << m() << std::endl;
         std::cerr << "SKEWEDNESS = " << SKEWEDNESS() << std::endl;
         std::cerr << "B1 = " << B1() << std::endl;
         std::cerr << "L1 = " << L1() << std::endl;
         std::cerr << "LP1 = " << LP1() << std::endl;
         std::cerr << "B2 = " << B2() << std::endl;
         std::cerr << "L2 = " << L2() << std::endl;
         std::cerr << "LP2 = " << LP2() << std::endl;
         std::cerr << "MIN_A = " << MIN_A() << std::endl;
         std::cerr << "MAX_A = " << MAX_A() << std::endl;
         std::cerr << "MIN_B = " << MIN_B() << std::endl;
         std::cerr << "MAX_B = " << MAX_B() << std::endl;
         std::cerr << "SIEVE_BOUND_ADJUSTMENT1 = " << SIEVE_BOUND_ADJUSTMENT1() << std::endl;
         std::cerr << "SIEVE_BOUND_ADJUSTMENT2 = " << SIEVE_BOUND_ADJUSTMENT2() << std::endl;
         std::cerr << "SMALL_PRIME_BOUND1 = " << SMALL_PRIME_BOUND1() << std::endl;
         std::cerr << "SMALL_PRIME_BOUND2 = " << SMALL_PRIME_BOUND2() << std::endl;
         std::cerr << "INITIAL_CUTOFF = " << INITIAL_CUTOFF() << std::endl;
         if (RESIEVE_) std::cerr << "RESIEVE = true" << std::endl;
         else std::cerr << "RESIEVE = false" << std::endl;
         std::cerr << "RELATION_FILE = " << RELATION_FILE() << std::endl;
         if (DEBUG_) std::cerr << "DEBUG = true" << std::endl;
         else std::cerr << "DEBUG = false" << std::endl;
         if (FIXED_SIEVE_REGION_) std::cerr << "FIXED_SIEVE_REGION = true" << std::endl;
         else std::cerr << "FIXED_SIEVE_REGION = false" << std::endl;
      }

   private:
      std::string SIEVE_ID_;
      VeryLong N_;
      Polynomial<VeryLong> f1_;
      Polynomial<VeryLong> f2_;
      VeryLong m_;
      long int B1_;
      long int L1_;
      long int LP1_;
      long int B2_;
      long int L2_;
      long int LP2_;
      double MIN_A_;
      double MAX_A_;
      long int MIN_B_;
      long int MAX_B_;
      long int SIEVE_BOUND_ADJUSTMENT1_;
      long int SIEVE_BOUND_ADJUSTMENT2_;
      long int SMALL_PRIME_BOUND1_;
      long int SMALL_PRIME_BOUND2_;
      long int INITIAL_CUTOFF_;
      double SKEWEDNESS_;
      bool RESIEVE_;
      std::string RELATION_FILE_;
      bool DEBUG_;
      bool FIXED_SIEVE_REGION_;
};

#endif
