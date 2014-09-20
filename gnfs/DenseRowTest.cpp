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
      for (DenseRow::const_iterator it = dr.begin();
            it != dr.end();
            ++it)
      {
         std::cout << " " << *it;
      }
      std::cout << std::endl;
   }
#endif
   return 0;
}
