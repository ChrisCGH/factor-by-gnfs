#include <string>
#include <fstream>
#include "VeryLong.h"

int main(int argc, char* argv[])
{
    std::string filename("nums.txt");
    if (argc > 1)
    {
        filename = argv[1];
    }
    std::fstream infile(filename.c_str(), std::ios::in);

    std::string numstr;
    std::vector<long int> factors;
    while (std::getline(infile, numstr))
    {
        FastVeryLong num(const_cast<char*>(numstr.c_str()));
        num.factorise_no_trial(&factors);
    }
    return 0;
}
