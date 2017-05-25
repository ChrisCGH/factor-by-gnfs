#ifndef ROOTCONFIG_H
#define ROOTCONFIG_H

#include "Polynomial.h"
#include <fstream>
#include <string>
// class for configuration information for
// square root procedure

class RootConfig
{
public:
    RootConfig(const char* filename)
    {
        ROOT_ID_ = "ctc4";
        N_ = "188198812920607963838697239461650439807163563379417382700763356422988859715234665485319060606504743045317388011303396716199692321205734031879550656996221305168759307650257059";
        RELATION_FILE_ = "ctc4.rel1";
        DEBUG_ = false;
        DUMP_FILE_ = "";

        std::ifstream config_file(filename, std::ios::in);
        if (config_file)
        {
            std::string str;
            while (getline(config_file, str))
            {
                // nothing
                if (str.empty()) continue;
                // comments
                if (str[0] == '#') continue;
                // we expect to find =
                std::string::size_type eqpos = str.find('=');
                if (std::string::npos == eqpos) continue;

                if (eqpos + 2 >= str.size()) continue;
                std::string s = str.substr(eqpos + 2);

                if (str.find("ROOT_ID = ") == 0)
                {
                    ROOT_ID_ = s;
                }
                else if (str.find("N = ") == 0)
                {
                    N_ = s;
                }
                else if (str.find("f1 = ") == 0)
                {
                    f1_ = Polynomial<VeryLong>::read_polynomial(s.c_str());
                }
                else if (str.find("f2 = ") == 0)
                {
                    f2_ = Polynomial<VeryLong>::read_polynomial(s.c_str());
                }
                else if (str.find("m = ") == 0)
                {
                    m_ = s;
                }
                else if (str.find("EXTRA_PRIME = ") == 0)
                {
                    EXTRA_PRIMES_.push_back(VeryLong(s));
                }
                else if (str.find("RELATION_FILE = ") == 0)
                {
                    RELATION_FILE_ = s;
                }
                else if (str.find("DEBUG = ") == 0)
                {
                    if (s == "true") DEBUG_ = true;
                }
                else if (str.find("DUMP_FILE = ") == 0)
                {
                    DUMP_FILE_ = s;
                }
            }
        }
        // check
        VeryLong x = f1_.evaluate(m_);
        VeryLong y = f2_.evaluate(m_);
        if (x % N_ != 0L || y % N_ != 0L)
        {
            std::cerr << "Problem: f1(m) = " << x << ", f2(m) = " << y << std::endl;
        }
    }

    std::string ROOT_ID()
    {
        return ROOT_ID_;
    }
    VeryLong N()
    {
        return N_;
    }
    Polynomial<VeryLong> f1()
    {
        return f1_;
    }
    Polynomial<VeryLong> f2()
    {
        return f2_;
    }
    VeryLong m()
    {
        return m_;
    }
    std::string RELATION_FILE()
    {
        return RELATION_FILE_;
    }
    int EXTRA_PRIMES()
    {
        return EXTRA_PRIMES_.size();
    }
    VeryLong EXTRA_PRIME(int i)
    {
        return EXTRA_PRIMES_[i];
    }
    bool DEBUG()
    {
        return DEBUG_;
    }
    std::string DUMP_FILE()
    {
        return DUMP_FILE_;
    }

    void display()
    {
        std::cout << "# Configuration options: " << std::endl;
        std::cout << "ROOT_ID = " << ROOT_ID() << std::endl;
        std::cout << "N = " << N() << std::endl;
        std::cout << "f1 = " << f1() << std::endl;
        std::cout << "f2 = " << f2() << std::endl;
        std::cout << "m = " << m() << std::endl;
        std::cout << "RELATION_FILE = " << RELATION_FILE() << std::endl;
        for (size_t i = 0; i < EXTRA_PRIMES_.size(); i++)
        {
            std::cout << "EXTRA_PRIME = " << EXTRA_PRIMES_[i] << std::endl;
        }
        if (DEBUG_)
        {
            std::cout << "DEBUG = true" << std::endl;
        }
        else
        {
            std::cout << "DEBUG = false" << std::endl;
        }
        std::cout << "DUMP_FILE = " << DUMP_FILE() << std::endl;
    }

private:
    std::string ROOT_ID_;
    VeryLong N_;
    Polynomial<VeryLong> f1_;
    Polynomial<VeryLong> f2_;
    VeryLong m_;
    std::string RELATION_FILE_;
    std::vector<VeryLong> EXTRA_PRIMES_;
    bool DEBUG_;
    std::string DUMP_FILE_;
};

#endif
