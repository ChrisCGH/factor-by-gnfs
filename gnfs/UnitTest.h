#ifndef UNITTEST_H
#define UNITTEST_H

#include <limits>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

class UnitTest
{
public:
    UnitTest() : passes_(0), failures_(0) {}

    void check(bool expression, const std::string& msg)
    {
        if (expression)
        {
            ++passes_;
        }
        else
        {
            std::cout << "Test failed: " << msg << std::endl;
            ++failures_;
        }
    }

    void test_summary()
    {
        std::cout << "Tests passed : " << passes_ << ", tests failed : " << failures_ << std::endl;
    }

    static bool compare_double(double d1, double d2)
    {
        if (std::fabs(d1 - d2) < std::numeric_limits<float>::epsilon())
        {
            return true;
        }
        else
        {
            std::cerr << std::setprecision(15) << "d1 = " << d1 << ", d2 = " << d2 << std::endl;
            return false;
        }
    }

private:
    size_t passes_;
    size_t failures_;
};

#endif
