#include <iostream>
#include <fstream>
#include <string>
#include "MemoryMappedFile.h"


int main(int argc, char* argv[])
{
   if (argc == 2)
   {
      MemoryMappedFile* Relmmfile = 0;
      Relmmfile = new MemoryMappedFile(argv[1]);

      std::string str;
      int count = 0;
      while (getline(*Relmmfile, str))
      {
         if (count%10000 == 0) std::cout << count << std::endl;
         ++count;
      }
      std::cout << count << std::endl;
      Relmmfile->reset();
      count = 0;
      while (getline(*Relmmfile, str))
      {
         if (count%10000 == 0) std::cout << count << std::endl;
         ++count;
      }
      std::cout << count << std::endl;
      delete Relmmfile;
   }
   else
   {
      std::fstream Relfile(argv[1]);

      std::string str;
      int count = 0;
      while (std::getline(Relfile, str))
      {
         if (count%10000 == 0) std::cout << count << std::endl;
         ++count;
      }
      std::cout << count << std::endl;
   }
   return 0;
}
