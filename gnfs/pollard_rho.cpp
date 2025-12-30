// Pollard's rho factorization algorithm
// Fast algorithm for finding small factors of composite numbers
#include <cstdlib>
#include <climits>
#include <iostream>
#include "pollard_rho.h"
#include "gcd.h"

namespace {
    // Compute (a * b) % mod avoiding overflow using 128-bit arithmetic if available
    inline unsigned long long int mulmod(unsigned long long int a, unsigned long long int b, unsigned long long int mod) {
#if defined(__SIZEOF_INT128__)
        // Use 128-bit arithmetic if available
        return (unsigned long long int)(((unsigned __int128)a * b) % mod);
#else
        // Fallback: use slow but safe method
        unsigned long long int result = 0;
        a %= mod;
        while (b > 0) {
            if (b & 1) {
                result = (result + a) % mod;
            }
            a = (a + a) % mod;
            b >>= 1;
        }
        return result;
#endif
    }

    // Polynomial function: f(x) = (x^2 + c) % n
    // Using c=1 is common and works well in practice
    inline unsigned long long int f(unsigned long long int x, unsigned long long int c, unsigned long long int n) {
        return (mulmod(x, x, n) + c) % n;
    }
}

bool pollard_rho(unsigned long long int N, long int& factor, bool debug)
{
    // Handle trivial cases
    if (N <= 1) return false;
    if (N == 2) {
        factor = 2;
        return true;
    }
    if (N % 2 == 0) {
        factor = 2;
        return true;
    }
    
    // For small N, just check a few small primes
    if (N < 100) {
        for (long int p = 3; p * p <= N; p += 2) {
            if (N % p == 0) {
                factor = p;
                return true;
            }
        }
        return false; // N is prime
    }
    
    // Try different starting points and constants if first attempt fails
    const int max_attempts = 5;
    for (int attempt = 0; attempt < max_attempts; ++attempt) {
        unsigned long long int x = 2;
        unsigned long long int y = 2;
        unsigned long long int c = 1 + attempt; // Try different c values
        unsigned long long int d = 1;
        
        // Floyd's cycle detection
        int iterations = 0;
        const int max_iterations = 100000; // Prevent infinite loops
        
        while (d == 1 && iterations < max_iterations) {
            x = f(x, c, N);
            y = f(f(y, c, N), c, N);
            
            // Compute GCD
            unsigned long long int abs_diff = (x > y) ? (x - y) : (y - x);
            d = gcd<unsigned long long int>(abs_diff, N);
            
            iterations++;
            
            // Batch GCD computation every 100 iterations for efficiency
            if (iterations % 100 == 0 && debug) {
                std::cerr << "Pollard rho: iteration " << iterations 
                          << ", attempt " << (attempt + 1) << std::endl;
            }
        }
        
        if (d != N && d != 1) {
            // Found a non-trivial factor
            if (d <= (unsigned long long int)LONG_MAX) {
                factor = (long int)d;
                if (debug) {
                    std::cerr << "Pollard rho found factor " << factor 
                              << " of " << N << " in " << iterations 
                              << " iterations (attempt " << (attempt + 1) << ")" << std::endl;
                }
                return true;
            } else {
                // Factor too large for long int, but we found it
                if (debug) {
                    std::cerr << "Pollard rho found factor too large for long int: " 
                              << d << std::endl;
                }
                return false;
            }
        }
        
        if (debug && d == N) {
            std::cerr << "Pollard rho: attempt " << (attempt + 1) 
                      << " failed (found trivial factor)" << std::endl;
        }
    }
    
    // All attempts failed
    if (debug) {
        std::cerr << "Pollard rho failed to factor " << N 
                  << " after " << max_attempts << " attempts" << std::endl;
    }
    return false;
}
