#include <iostream>
#include <fstream>
#include <vector>
#include "convert.h"
#include "SieveUtils.h"
#include "SieveConfig.h"

namespace SieveUtils
{
bool checkForRelation(const VeryLong& a, const VeryLong& b,
                      const Polynomial<VeryLong>& f,
                      long int B, long int L, long int LP, const VeryLong& L_LP,
                      std::vector<VeryLong>& factors, const std::vector<long int>& primes)
{
    VeryLong value = f.evaluate_homogeneous(a, b);
    //std::cerr << "f(a,b) = " << value << std::endl;
    factors.clear();
    VeryLong tmp = value;
    if (tmp < 0L) tmp = -tmp;
    int large_primes = 0;
    const VeryLong one(1L);
    // Zeroth of all remove primes from f(a,b)
    for (size_t i = 0; i < primes.size(); ++i)
    {
        if (tmp % primes[i] == 0L)
        {
            tmp /= primes[i];
            factors.push_back(VeryLong(primes[i]));
        }
    }

    // First of all use trial division
    VeryLong tmp1;
    std::vector<long int> ifactors;
    VeryLong::set_max_prime(B);
    if (B > VeryLong::get_max_prime())
    {
        tmp.factorise_trial_division(&ifactors, &tmp1, B);
    }
    else
    {
        tmp.factorise_trial_division(&ifactors, &tmp1);
    }
    if (tmp1 == one || (tmp1 <= B && tmp1.is_probable_prime()))
    {
        for (size_t i = 0; i < ifactors.size(); i++)
        {
            factors.push_back(VeryLong(ifactors[i]));
        }
        if (tmp1 != one && tmp1 <= B) factors.push_back(tmp1);
        return true;
    }
    if (tmp1 > B && tmp1 <= L && LP > 0 && tmp1.is_probable_prime())
    {
        for (size_t i = 0; i < ifactors.size(); i++)
        {
            factors.push_back(VeryLong(ifactors[i]));
        }
        if (tmp1 != one && tmp1 <= L) factors.push_back(tmp1);
        return true;
    }
    // here if tmp1 != 1 and (tmp1 > L or tmp1 is composite)
    if (LP == 0)
    {
        std::cerr << "LP == 0, tmp1 = " << tmp1 << ", B = " << B << ", L = " << L << std::endl;
        return false;
    }
    if (tmp1 <= L)
    {
        // tmp1 is a composite <= L, so split it
        // By definition, if we are here tmp1 has factors larger than MAX_PRIME
        // in VeryLong.cpp, but less than L
        for (size_t i = 0; i < ifactors.size(); i++)
        {
            factors.push_back(VeryLong(ifactors[i]));
        }
        std::cerr << "1. Calling factorise on " << tmp1 << std::endl;
        tmp1.factorise_no_trial(&factors);
        std::cerr << "1. done" << std::endl;
        return true;
    }
    // here if tmp1 > L
    // we are allowing LP large primes <= L
    if (tmp1.is_probable_prime())
    {
        std::cerr << tmp1 << " is probable prime > " << L << std::endl;
        return false;
    }
    // Note that we could have MAX_PRIME < B, in which case tmp1 may still contain some
    // factors > MAX_PRIME and < B
    if (tmp1 / L_LP > B)
    {
        std::cerr << "not smooth : tmp1 = " << tmp1 << ", L_LP = " << L_LP << ", tmp1 / L_LP = " << tmp1 / L_LP << ", B = " << B << std::endl;
        return false;
    }
    //std::cout << "tmp1 = " << tmp1 << ", L_LP = " << L_LP << ", tmp1 / L_LP = " << tmp1 / L_LP << ", B = " << B << std::endl;
    std::vector<VeryLong> factors1;

    std::cerr << "2. Calling factorise on " << tmp1 << std::endl;
    tmp1.factorise_no_trial(&factors1);
    std::cerr << "2. done" << std::endl;

    for (size_t i = 0; i < factors1.size(); i++)
    {
        if (factors1[i] > L)
        {
            std::cerr << std::endl;
            std::cerr << "not smooth : " << L << " < " << factors1[i] << std::endl;
            //cout << "not smooth : " << L << " < " << factors1[i] << std::endl;
            return false;
        }
        if (factors1[i] > B && factors1[i] <= L) large_primes++;
    }
    //if (large_primes <= LP)
    {
        std::cerr << "smooth : " << std::endl;
        //cout << "smooth : " << std::endl;
        for (size_t i = 0; i < ifactors.size(); i++)
        {
            factors.push_back(VeryLong(ifactors[i]));
        }
        for (size_t i = 0; i < factors1.size(); i++)
        {
            factors.push_back(factors1[i]);
        }
        return true;
    }
    std::cerr << "falling through, large_primes = " << large_primes << ", LP = " << LP << std::endl;
    return false;
}

bool checkRelations(const char* relfile, const SieveConfig& config)
{
    std::fstream rels(relfile, std::ios::in);
    std::vector<VeryLong> factors1;
    std::vector<VeryLong> factors2;
    std::string str;
    Polynomial<VeryLong> f1_ = config.f1();
    Polynomial<VeryLong> f2_ = config.f2();
    long int B1_ = config.B1();
    long int L1_ = config.L1();
    long int LP1_ = config.LP1();
    long int B2_ = config.B2();
    long int L2_ = config.L2();
    long int LP2_ = config.LP2();
    VeryLong L_LP_1_ = L1_;
    for (int i = 0; i < LP1_; i++) L_LP_1_ *= L1_;
    VeryLong L_LP_2_ = L2_;
    for (int i = 0; i < LP2_; i++) L_LP_2_ *= L2_;
    while (getline(rels, str))
    {
        long long int a;
        long long int b;
        char* alg_str;
        char* rat_str;
        Convert::parse_FBGNFS(str, a, b, alg_str, rat_str);
        std::vector<long int> alg_primes;
        std::vector<long int> rat_primes;
        Convert::extract_FBGNFS_primes(alg_str, alg_primes);
        Convert::extract_FBGNFS_primes(rat_str, rat_primes);

        if (!checkForRelation(VeryLong(a), VeryLong(b), f2_, B2_, L2_, LP2_, L_LP_2_, factors2, rat_primes) ||
                !checkForRelation(VeryLong(a), VeryLong(b), f1_, B1_, L1_, LP1_, L_LP_1_, factors1, alg_primes))
        {
            return false;
        }
    }
    return true;
}

};
