#include "mt19937int.h"
#include "pselect.h"
#include <string.h>
#include <iostream>
#include <fstream>

int main(int argc, char** argv)
{
   if (argc <= 1)
   {
      std::vector<VeryLong> c(6);
      c[5] = "2494474680";
      c[4] = "-60025243754288938";
      c[3] = "-3373062207827200423874";
      c[2] = "35057246057064976169265284946";
      c[1] = "4645652748792665336283343904631161";
      c[0] = "-2514680780127600312183202859318919091241";
      double s = 896878.58;
      Polynomial<VeryLong> p(c);
      double als = average_log_size(p, s);
      double alpha = alpha_F(p, 2000, 200);
      std::cout << "p = " << p << std::endl;
      std::cout << "s = " << s << std::endl;
      std::cout << "als = " << als << std::endl;
      std::cout << "alpha = " << alpha << std::endl;
      std::cout << "E(F) = " << als + alpha << std::endl;

      double new_s = minimize_I_over_s(Polynomial<VeryLong>::convert_to_double(p), 0.0, 0.0, 0.0, 0.0, 0.0, 1000.0);
      std::cout << "new_s = " << new_s << std::endl;
      double new_als = average_log_size(p, new_s);
      std::cout << "new_als = " << new_als << std::endl;
      std::cout << "alpha = " << alpha << std::endl;
      std::cout << "new E(F) = " << new_als + alpha << std::endl;
   }
   else
   {
      std::fstream file(argv[1], std::ios::in);
      std::cout << "Reading polynomial from " << argv[1] << std::endl;
      Polynomial<VeryLong> p = Polynomial<VeryLong>::read_polynomial(file);
      std::cout << "p = " << p << std::endl;
      double alpha = alpha_F(p, 2000, 200);
      std::cout << "alpha = " << alpha << std::endl;
      double s = minimize_I_over_s(Polynomial<VeryLong>::convert_to_double(p), 0.0, 0.0, 0.0, 0.0, 0.0, 1000.0);
      std::cout << "s = " << s << std::endl;
      double als = average_log_size(p, s);
      std::cout << "als = " << als << std::endl;
      std::cout << "E(F) = " << als + alpha << std::endl;
   }
   return 0;
}
