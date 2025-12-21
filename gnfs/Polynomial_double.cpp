#include "Polynomial.h"
#include <cmath>

// Optimized Horner's method for double polynomial evaluation
template <>
double Polynomial<double>::evaluate(const double& value) const
{
    if (_coefficients.empty()) return 0.0;
    
    double result = _coefficients[_coefficients.size() - 1];
    for (int i = static_cast<int>(_coefficients.size()) - 2; i >= 0; --i)
    {
        result = result * value + _coefficients[i];
    }
    return result;
}

// Specialized homogeneous evaluation for common GNFS polynomial degrees
template <>
double Polynomial<double>::evaluate_homogeneous(const double& a, const double& b) const
{
    const int d = deg();
    
    // For homogeneous polynomial F(a,b), we want F(a/b) * b^deg
    // But for lattice sieving, b is often 1, so optimize for that
    if (b == 1.0)
    {
        return evaluate(a);
    }
    
    // Unrolled evaluation for common degrees
    switch (d)
    {
        case 1:
            return _coefficients[0] + _coefficients[1] * a;
            
        case 2:
        {
            double a2 = a * a;
            return _coefficients[0] + _coefficients[1] * a + _coefficients[2] * a2;
        }
        
        case 3:
        {
            double a2 = a * a;
            double a3 = a2 * a;
            return _coefficients[0] + _coefficients[1] * a + 
                   _coefficients[2] * a2 + _coefficients[3] * a3;
        }
        
        case 4:
        {
            double a2 = a * a;
            double a3 = a2 * a;
            double a4 = a2 * a2;
            return _coefficients[0] + _coefficients[1] * a + 
                   _coefficients[2] * a2 + _coefficients[3] * a3 + 
                   _coefficients[4] * a4;
        }
        
        case 5:
        {
            // Use FMA if available - Horner's method
            #ifdef __FMA__
            double result = _coefficients[5];
            result = std::fma(result, a, _coefficients[4]);
            result = std::fma(result, a, _coefficients[3]);
            result = std::fma(result, a, _coefficients[2]);
            result = std::fma(result, a, _coefficients[1]);
            result = std::fma(result, a, _coefficients[0]);
            return result;
            #else
            double result = _coefficients[5];
            result = result * a + _coefficients[4];
            result = result * a + _coefficients[3];
            result = result * a + _coefficients[2];
            result = result * a + _coefficients[1];
            result = result * a + _coefficients[0];
            return result;
            #endif
        }
        
        case 6:
        {
            // Use FMA if available - Horner's method
            #ifdef __FMA__
            double result = _coefficients[6];
            result = std::fma(result, a, _coefficients[5]);
            result = std::fma(result, a, _coefficients[4]);
            result = std::fma(result, a, _coefficients[3]);
            result = std::fma(result, a, _coefficients[2]);
            result = std::fma(result, a, _coefficients[1]);
            result = std::fma(result, a, _coefficients[0]);
            return result;
            #else
            double result = _coefficients[6];
            result = result * a + _coefficients[5];
            result = result * a + _coefficients[4];
            result = result * a + _coefficients[3];
            result = result * a + _coefficients[2];
            result = result * a + _coefficients[1];
            result = result * a + _coefficients[0];
            return result;
            #endif
        }
        
        default:
            // Fall back to generic implementation for unusual degrees
            return evaluate_homogeneous_generic(a, b);
    }
}
