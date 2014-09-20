#ifndef CONFIG_H
#define CONFIG_H

#include "VeryLong.h"
#include <fstream>
// class for configuration information for
// skewed polynomial selection procedure

class Skewed_selection_config
{
   public:
      Skewed_selection_config(const char* filename)
      {
         N_ = "188198812920607963838697239461650439807163563379417382700763356422988859715234665485319060606504743045317388011303396716199692321205734031879550656996221305168759307650257059";
         DEGREE_ = 5;
         MIN_AD_ = "10000000000";
         MAX_AD_ = "100000000000000";
         MAX_FRACTION_ = 5e-6;
         MAX_ALS_ = 54.5;
         MAX_J0_ = 10000;
         MAX_J1_ = 100;
         MAX_SMALL_PRIME_ = 100;
         ALPHA_CUTOFF_ = -4.5;
         C_START_ = 4500L;
         C_RESTART_ = 0L;
         C_FACTOR_ = 1000L;
         REPEAT_CUTOFF_ = 48.5;
         GOOD_M_CUTOFF_ = 0.6;
         OUTPUT_FILE_ = "skewed.out";
         PRINTING_BOUND_ = 47.0;
         NON_MONIC_ = true;
         DEBUG_ = false;

         std::ifstream config_file(filename, std::ios::in);
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

               std::string s = str.substr(eqpos + 2);
               if (str.find("N = ") == 0)
               {
                  N_ = s;
               }
               else if (str.find("DEGREE = ") == 0)
               {
                  DEGREE_ = std::atoi(s.c_str());
               }
               else if (str.find("MIN_AD = ") == 0)
               {
                  MIN_AD_ = s;
               }
               else if (str.find("MAX_AD = ") == 0)
               {
                  MAX_AD_ = s;
               }
               else if (str.find("MAX_FRACTION = ") == 0)
               {
                  MAX_FRACTION_ = std::atof(s.c_str());
               }
               else if (str.find("MAX_ALS = ") == 0)
               {
                  MAX_ALS_ = std::atof(s.c_str());
               }
               else if (str.find("MAX_J0 = ") == 0)
               {
                  MAX_J0_ = std::atol(s.c_str());
               }
               else if (str.find("MAX_J1 = ") == 0)
               {
                  MAX_J1_ = std::atol(s.c_str());
               }
               else if (str.find("MAX_SMALL_PRIME = ") == 0)
               {
                  MAX_SMALL_PRIME_ = std::atol(s.c_str());
               }
               else if (str.find("ALPHA_CUTOFF = ") == 0)
               {
                  ALPHA_CUTOFF_ = std::atof(s.c_str());
               }
               else if (str.find("C_START = ") == 0)
               {
                  C_START_ = s;
               }
               else if (str.find("C_RESTART = ") == 0)
               {
                  C_RESTART_ = s;
               }
               else if (str.find("C_FACTOR = ") == 0)
               {
                  C_FACTOR_ = std::atol(s.c_str());
               }
               else if (str.find("REPEAT_CUTOFF = ") == 0)
               {
                  REPEAT_CUTOFF_ = std::atof(s.c_str());
               }
               else if (str.find("PRINTING_BOUND = ") == 0)
               {
                  PRINTING_BOUND_ = std::atof(s.c_str());
               }
               else if (str.find("GOOD_M_CUTOFF = ") == 0)
               {
                  GOOD_M_CUTOFF_ = std::atof(s.c_str());
               }
               else if (str.find("OUTPUT_FILE = ") == 0)
               {
                  OUTPUT_FILE_ = s;
               }
               else if (str.find("NON_MONIC = ") == 0)
               {
                  if (s == "true") NON_MONIC_ = true;
                  else NON_MONIC_ = false;
               }
               else if (str.find("DEBUG = ") == 0)
               {
                  if (s == "true") DEBUG_ = true;
                  else DEBUG_ = false;
               }
            }
         }

      }
      VeryLong N() const
      {
         return N_;
      }
      int DEGREE() const
      {
         return DEGREE_;
      }
      VeryLong MIN_AD() const
      {
         return MIN_AD_;
      }
      VeryLong MAX_AD() const
      {
         return MAX_AD_;
      }
      double MAX_FRACTION() const
      {
         return MAX_FRACTION_;
      }
      double MAX_ALS() const
      {
         return MAX_ALS_;
      }
      long int MAX_J0() const
      {
         return MAX_J0_;
      }
      long int MAX_J1() const
      {
         return MAX_J1_;
      }
      long int MAX_SMALL_PRIME() const
      {
         return MAX_SMALL_PRIME_;
      }
      std::string OUTPUT_FILE() const
      {
         return OUTPUT_FILE_;
      }
      double REPEAT_CUTOFF() const
      {
         return REPEAT_CUTOFF_;
      }
      double GOOD_M_CUTOFF() const
      {
         return GOOD_M_CUTOFF_;
      }
      double ALPHA_CUTOFF() const
      {
         return ALPHA_CUTOFF_;
      }
      VeryLong C_START() const
      {
         return C_START_;
      }
      VeryLong C_RESTART() const
      {
         return C_RESTART_;
      }
      long int C_FACTOR() const
      {
         return C_FACTOR_;
      }
      double PRINTING_BOUND() const
      {
         return PRINTING_BOUND_;
      }
      bool NON_MONIC() const
      {
         return NON_MONIC_;
      }
      bool DEBUG() const
      {
         return DEBUG_;
      }

      void display() const
      {
         std::cout << "# Configuration options: " << std::endl;
         std::cout << "N = " << N() << std::endl;
         std::cout << "DEGREE = " << DEGREE() << std::endl;
         std::cout << "MIN_AD = " << MIN_AD() << std::endl;
         std::cout << "MAX_AD = " << MAX_AD() << std::endl;
         std::cout << "MAX_FRACTION = " << MAX_FRACTION() << std::endl;
         std::cout << "GOOD_M_CUTOFF = " << GOOD_M_CUTOFF() << std::endl;
         std::cout << "MAX_ALS = " << MAX_ALS() << std::endl;
         std::cout << "MAX_J0 = " << MAX_J0() << std::endl;
         std::cout << "MAX_J1 = " << MAX_J1() << std::endl;
         std::cout << "MAX_SMALL_PRIME  = " << MAX_SMALL_PRIME() << std::endl;
         std::cout << "ALPHA_CUTOFF = " << ALPHA_CUTOFF() << std::endl;
         std::cout << "PRINTING_BOUND = " << PRINTING_BOUND() << std::endl;
         std::cout << "REPEAT_CUTOFF = " << REPEAT_CUTOFF() << std::endl;
         std::cout << "C_START = " << C_START() << std::endl;
         std::cout << "C_RESTART = " << C_RESTART() << std::endl;
         std::cout << "C_FACTOR = " << C_FACTOR() << std::endl;
         std::cout << "NON_MONIC = " << NON_MONIC() << std::endl;
         std::cout << "DEBUG = " << DEBUG() << std::endl;
         std::cout << "OUTPUT_FILE = " << OUTPUT_FILE() << std::endl;
      }

   private:
      VeryLong N_;
      int DEGREE_;
      VeryLong MIN_AD_;
      VeryLong MAX_AD_;
      double MAX_FRACTION_;
      double MAX_ALS_;
      long int MAX_J0_;
      long int MAX_J1_;
      long int MAX_SMALL_PRIME_;
      std::string OUTPUT_FILE_;
      double ALPHA_CUTOFF_;
      double PRINTING_BOUND_;
      double REPEAT_CUTOFF_;
      VeryLong C_START_;
      VeryLong C_RESTART_;
      long int C_FACTOR_;
      double GOOD_M_CUTOFF_;
      bool NON_MONIC_;
      bool DEBUG_;
};

#endif
