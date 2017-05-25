#ifndef _QUADRATICCHARACTERS_H
#define _QUADRATICCHARACTERS_H
#include <vector>
#include "VeryLong.h"
#include "Polynomial.h"

class QuadraticCharacters
{
public:
    QuadraticCharacters(const Polynomial<VeryLong>& f, long int B);
    ~QuadraticCharacters()
    {}
    void generate(const VeryLong& a, const VeryLong& b, std::vector<char>& qcs);

private:
    Polynomial<VeryLong> f_;
    long int B_;
    std::vector<std::pair<VeryLong, VeryLong> > qs_;
    void initialise();

};
#endif
