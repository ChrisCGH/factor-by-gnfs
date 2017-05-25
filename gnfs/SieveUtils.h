#ifndef SIEVEUTILS_H
#define SIEVEUTILS_H
#include <vector>
#include "VeryLong.h"
#include "Polynomial.h"

class SieveConfig;
namespace SieveUtils
{
bool checkRelations(const char* relfile, const SieveConfig& config);
bool checkForRelation(const VeryLong& a, const VeryLong& b,
                      const Polynomial<VeryLong>& f,
                      long int B, long int L, long int LP, const VeryLong& L_LP,
                      std::vector<VeryLong>& factors, const std::vector<long int>& primes);
};
#endif
