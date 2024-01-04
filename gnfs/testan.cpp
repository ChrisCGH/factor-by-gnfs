#include "NumberField.h"
#include "AlgebraicNumber.h"
#include "AlgebraicNumber_in_O_pO.h"
#include "Ideal.h"
#include "pow.h"
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include "MPFloat.h"
#include <fstream>
#include "root.h"
#include "RootConfig.h"
#include <string>
#include <ostream>

template<> Matrix<Quotient<VeryLong> > AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular>::W_mult_;
template<> Matrix<Quotient<VeryLong> > AlgebraicNumber_in_O_pO_<long, VeryLong, LongModular>::W_mult_;
template<> Matrix<VeryLongModular> AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular>::M_;
template<> Matrix<LongModular> AlgebraicNumber_in_O_pO_<long, VeryLong, LongModular>::M_;
template<> VeryLong AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular>::p_;
template<> long AlgebraicNumber_in_O_pO_<long, VeryLong, LongModular>::p_;
template<> VeryLongModular AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular>::w01_;
template<> LongModular AlgebraicNumber_in_O_pO_<long, VeryLong, LongModular>::w01_;
template<> VeryLongModular AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular>::w11_;
template<> LongModular AlgebraicNumber_in_O_pO_<long, VeryLong, LongModular>::w11_;
template<> bool AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular>::optimisation_ok_;
template<> bool AlgebraicNumber_in_O_pO_<long, VeryLong, LongModular>::optimisation_ok_;

namespace
{
long long strtoll(const char* str)
{
    double x = std::atof(str);
    return (long long)x;
}
}

int main(int argc, char** argv)
{
    RootConfig config("root.cfg");
    config.display();
    std::string relfile_str = config.RELATION_FILE();
    std::string ratRelFile(relfile_str);
    const char* relfile = relfile_str.c_str();
    if (argc > 1) relfile = argv[1];

    Polynomial<VeryLong> f1 = config.f1();
    std::cout << "f1 = " << f1 << std::endl;

    for (int i = 0; i < config.EXTRA_PRIMES(); i++)
    {
        VeryLong::addPrime(config.EXTRA_PRIME(i));
    }
    char fbFile[132];
    strcpy(fbFile, config.ROOT_ID().c_str());
    strcat(fbFile, ".fb.dat");
    NumberField nf(f1, fbFile);

    AlgebraicNumber::setNumberField(nf);

    VeryLong p(23L);
    AlgebraicNumber_in_O_pO::set_basis(p);

    std::vector<std::pair<PrimeIdeal, int> > primeIdeal;
    PrimeIdeal::primeDecomposition(p, primeIdeal);

    for (size_t i = 0; i < primeIdeal.size(); ++i)
    {
        std::cout << primeIdeal[i].second << std::endl << primeIdeal[i].first << std::endl;
        std::cout << "********************************************" << std::endl;
    }

    long long int aa = -2251919LL;
    long int b = 48342L;
    AlgebraicNumber an(aa, b);
    AlgebraicNumber_in_O_pO an_(an);
    AlgebraicNumber_in_O_pO an1_(aa, b);

    if (an_ != an1_)
    {
        std::cout << "(a,b) = (" << aa << "," << b << ")" << std::endl;
        std::cout << "an = " << an << std::endl;
        std::cout << "an_ = " << an_ << std::endl;
        std::cout << "an1_ = " << an1_ << std::endl;
    }

    std::fstream infile(relfile, std::ios::in);
    // file format should be
    // a b
    std::string str;
    int done = 0;
    int line = 0;
    bool numeratorRelation = true;
    const std::string NumeratorStr("Numerator");
    const std::string DenominatorStr("Denominator");
    static char* buf = 0;
    static std::string::size_type buflen = 0;
    while (!done)
    {
        getline(infile, str);
        if (str.empty()) continue;
        if (str.size() > buflen)
        {
            delete [] buf;
            buflen = str.size();
            buf = new char [ buflen + 1 ];
        }
        strcpy(buf, str.c_str());
        if (buf[0] == '!') done = 1;
        else
        {
            line++;
            if (NumeratorStr == buf) continue;
            if (DenominatorStr == buf)
            {
                if (numeratorRelation)
                {
                    numeratorRelation = false;
                }
                continue;
            }
            // find space
            char* c = buf;
            while (c && *c && *c != ' ') c++;
            if (!c || !*c || !*(c+1))
            {
                std::cerr << "Problem: bad format in line " << line << ":" << std::endl;
                std::cerr << "[" << buf << "]" << std::endl;
                return 0;
            }
            *c = '\0';
            c++;
            long long int aa = strtoll(buf);
            long int b = std::atol(c);
            AlgebraicNumber an(aa, b);
            AlgebraicNumber_in_O_pO an_(an);
            AlgebraicNumber_in_O_pO an1_(aa, b);

            if (an_ != an1_)
            {
                std::cout << "(a,b) = (" << aa << "," << b << ")" << std::endl;
                std::cout << "an = " << an << std::endl;
                std::cout << "an_ = " << an_ << std::endl;
                std::cout << "an1_ = " << an1_ << std::endl;
            }
        }
    }
    return 0;
}

