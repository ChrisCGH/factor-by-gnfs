#include "blockLanczos.h"
#include <iostream>
#include <fstream>
#include <string>

namespace
{
void usage()
{
   std::cerr << "Usage: bl -m matrix_file [ -c checkpoint_file ] [ -i checkpoint_interval ] [ -v validation_interval ]" << std::endl;
   std::exit(1);
}
};
int main(int argc, char* argv[])
{
   std::string matrix_file = "";
   std::string checkpoint_file = "";
   int checkpoint_interval = 0;
   int validation_interval = 0;
   // If set to true, split off the densest rows of the matrix, to process at the end
   bool split = false;
   int arg = 1;
   while (arg < argc)
   {
      if (strcmp(argv[arg], "-m") == 0)
      {
         ++arg;
         matrix_file = argv[arg];
      }
      else if (strcmp(argv[arg], "-c") == 0)
      {
         ++arg;
         checkpoint_file = argv[arg];
      }
      else if (strcmp(argv[arg], "-i") == 0)
      {
         ++arg;
         checkpoint_interval = std::atoi(argv[arg]);
      }
      else if (strcmp(argv[arg], "-v") == 0)
      {
         ++arg;
         validation_interval = std::atoi(argv[arg]);
      }
      else if (strcmp(argv[arg], "-split") == 0)
      {
         split = true;
      }
      else
      {
         usage();
      }
      arg++;
   }
   if (matrix_file == "") usage();

   BlockLanczos bl(matrix_file, checkpoint_file, checkpoint_interval, validation_interval, split);

   try
   {
      BITMATRIX kerL;
      BITMATRIX kerR;
      bl.kernel(kerL, kerR);

      kerL.printTranspose(std::cout);
      kerR.printTranspose(std::cout);
   }
   catch (const char* s)
   {
      std::cerr << "Caught exception: [" << s << "]" << std::endl;
   }
   catch (const std::string& s)
   {
      std::cerr << "Caught exception: [" << s << "]" << std::endl;
   }
   return 0;
}
