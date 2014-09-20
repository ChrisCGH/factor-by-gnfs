#include "VeryLong.h"
#include "VeryLongModular.h"
#include "Polynomial.h"
#include "PolynomialOptimizer.h"
#include "MPFloat.h"
#include "pselect.h"
#include <iomanip>
#include <string>
#include <cstdlib>

int main(int argc, char* argv[])
{
    long int skew = 0;
    Polynomial<VeryLong> f1;
    Polynomial<VeryLong> f2;
    VeryLong N;

    int arg = 1;
    while (arg < argc)
    {
        if (argv[arg] == std::string("-s"))
        {
            ++arg;
            if (arg < argc)
            {
                skew = std::atol(argv[arg]);
            }
        } 
        else if (argv[arg] == std::string("-N"))
        {
            ++arg;
            if (arg < argc)
            {
                N = VeryLong(argv[arg]);
            }
        }
        else if (argv[arg] == std::string("-p1"))
        {
            ++arg;
            if (arg < argc)
            {
                f1 = Polynomial<VeryLong>::read_polynomial(argv[arg]);
            }
        }
        else if (argv[arg] == std::string("-p2"))
        {
            ++arg;
            if (arg < argc)
            {
                f2 = Polynomial<VeryLong>::read_polynomial(argv[arg]);
            }
        }
        ++arg;
    }

    std::cout << "f1 = " << f1 << std::endl;
    std::cout << "f2 = " << f2 << std::endl;
    MPFloat::set_precision(100);
    VeryLongModular::set_default_modulus(N);
    VeryLongModular tmp1 = VeryLongModular(f2.coefficient(0)) / VeryLongModular(-f2.coefficient(1));
    VeryLong m = tmp1.get_very_long();
    std::cout << "m = " << m << std::endl;
    VeryLong best_s(skew);
    double I_F_S = PolynomialOptimizer::average_log_size(f1, best_s);
    double alpha = PolynomialOptimizer::alpha_F(f1, 2000, 200);
    std::cout << "E(F) = " << I_F_S + alpha << std::endl;
    return 0;
}

