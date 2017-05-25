#ifndef EXCEPTIONAL_PRIMES_H
#define EXCEPTIONAL_PRIMES_H
#include "NumberField.h"
#include "Ideal.h"
#include <map>
#include <unordered_map>

namespace std
{
template <> struct hash<std::pair<PrimeIdeal*, int> >
{
    size_t operator()(const pair<PrimeIdeal*, int>& pip) const
    {
        size_t seed = hash<PrimeIdeal*>()(pip.first);
        seed ^= hash<int>()(pip.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};
};
class ExceptionalPrimes
{
public:
    ExceptionalPrimes(const NumberField* nf) : nf_(nf)
    {
        initialise();
    }
    ~ExceptionalPrimes()
    {}
    bool isExceptional(long int p)
    {
        if (index_ % p == 0L) return true;
        return false;
    }
    std::vector<PrimeIdealRep*>::const_iterator begin(long int p)
    {
        return exceptional_primes_[p].begin();
    }
    std::vector<PrimeIdealRep*>::const_iterator end(long int p)
    {
        return exceptional_primes_[p].end();
    }
    std::vector<long int> getListOfExceptionalPrimes()
    {
        std::vector<long int> theList;
        for (auto& p: exceptional_primes_)
        {
            theList.push_back(p.first);
        }
        return theList;
    }

    bool relationIsInExceptionalPrimePower(const VeryLong& a, const VeryLong& b, PrimeIdeal* pi, int e)
    {
        if (exceptional_primes_powers_.find(std::make_pair(pi, e)) ==
                exceptional_primes_powers_.end()) return false;
        Matrix<long int> H = exceptional_primes_powers_[std::make_pair(pi, e)];
        if (b % H(1,1) != 0L) return false;
        VeryLong v = -b / H(1,1);
        VeryLong u = a - v* H(0,1);
        if (u % H(0,0) != 0L) return false;
        return true;
    }

    int padicValuation(const VeryLong& a, const VeryLong& b, PrimeIdeal* pi)
    {
        int e = 0;
        const int max_power = 20;
        while (e < max_power && relationIsInExceptionalPrimePower(a, b, pi, e + 1)) ++e;
        if (e == max_power) e = -1;
        return e;
    }

private:
    typedef std::unordered_map<long int, std::vector<PrimeIdealRep*> > exceptional_primes_type;
    exceptional_primes_type exceptional_primes_;
    typedef std::unordered_map<std::pair<PrimeIdeal*, int>, Matrix<long int> > exceptional_primes_powers_type;
    exceptional_primes_powers_type exceptional_primes_powers_;

    const NumberField* nf_;
    VeryLong index_;
    AlgebraicNumber omega1_;
    long int x_;
    long int y_;

    void initialise()
    {
        const std::vector<AlgebraicNumber>& ib = AlgebraicNumber::integralBasis();
        omega1_ = ib[1];
        x_ = omega1_.coefficient(0).numerator().get_long();
        y_ = omega1_.coefficient(1).numerator().get_long();
        index_ = nf_->index();
        std::vector<VeryLong> factors;
        index_.factorise(&factors);
        VeryLong p;
        VeryLong prev_p(0L);
        const VeryLong LONG_MAX_VL(LONG_MAX);
        for (size_t i = 0; i < factors.size(); i++)
        {
            p = factors[i];
            if (p != prev_p && p < LONG_MAX_VL)
            {
                long int pp = p.get_long();
                std::vector<std::pair<PrimeIdeal, int> > primeIdeal;
                PrimeIdeal::primeDecomposition(p, primeIdeal);
                for (size_t j = 0; j < primeIdeal.size(); j++)
                {
                    PrimeIdeal* pi = new PrimeIdeal(primeIdeal[j].first);
                    //std::cout << "Exceptional prime pi = " << *pi << std::endl;
                    PrimeIdealRep* pir = new PrimeIdealRep(pi);
                    exceptional_primes_[pp].push_back(pir);
                    // generate powers of each prime ideal and store
                    const int max_power = 20;
                    Ideal prime_power(*pi);
                    for (int e = 1; e <= max_power; ++e)
                    {
                        Matrix<long int> H(2,2);
                        Matrix<VeryLong> hnf = prime_power.hnf_basis();
                        VeryLong den = prime_power.denominator();
                        VeryLong tmp = hnf(0,0) / den;
                        H(0,0) = tmp.get_long();
                        tmp = hnf(0,1) / den;
                        H(0,1) = tmp.get_long();
                        tmp = hnf(1,1) / den;
                        H(1,1) = tmp.get_long();

                        exceptional_primes_powers_[std::make_pair(pi, e)] = H;
                        if (e < max_power) prime_power *= *pi;
                    }
                }
            }
            prev_p = p;
        }
    }
};

#endif
