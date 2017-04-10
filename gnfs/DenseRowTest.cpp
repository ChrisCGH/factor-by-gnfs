#include "SparseMatrix.h"
#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char* argv[])
{
#ifdef DENSEROW
   std::fstream fin(argv[1], std::ios::in);
   std::string str;
   while (std::getline(fin, str))
   {
      DenseRow dr(str);
      std::cout << dr.size();
      for (auto& r: dr)
      {
         std::cout << " " << r;
      }
      std::cout << std::endl;
   }
#endif
   return 0;
}
