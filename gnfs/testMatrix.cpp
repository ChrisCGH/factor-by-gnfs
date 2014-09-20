#include "SparseMatrix.h"
#include "MemoryMappedFile.h"
#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char* argv[])
{
   if (argc > 1)
   {
      SparseMatrix3 A(argv[1]);
      std::cout << "A has " << A.rows() << " rows and " << A.cols() << " columns" << std::endl;

      std::fstream f("testMatrix.out", std::ios::out);
      f << A << std::flush;

      BitMatrix X2;
      {
         std::fstream fX2("save/X2.dat", std::ios::in);
         X2.read(fX2);
      }

      BitMatrix AX2;
      {
         std::fstream fAX2("save/AX2.dat", std::ios::in);
         AX2.read(fAX2);
      }

      BitMatrix actual_AX2;
      multiply(A, X2, actual_AX2);

      if (AX2 != actual_AX2)
      {
         std::cerr << "Problem" << std::endl;
         std::fstream result("save/actual_AX2.dat", std::ios::out);
         actual_AX2.write(result);
      }
      else
      {
         std::cerr << "multiply ok" << std::endl;
      }

      BitMatrix X;
      {
         std::fstream fX("save/X.dat", std::ios::in);
         X.read(fX);
      }

      BitMatrix AtX;
      {
         std::fstream fAtX("save/AtX.dat", std::ios::in);
         AtX.read(fAtX);
      }

      BitMatrix actual_AtX;


      multiplyt(A, X, actual_AtX);

      if (AtX != actual_AtX)
      {
         std::cerr << "Problem" << std::endl;
         std::fstream result("save/actual_AtX.dat", std::ios::out);
         actual_AtX.write(result);
      }
      else
      {
         std::cerr << "multiplyt ok" << std::endl;
      }

   }
   return 0;
}
