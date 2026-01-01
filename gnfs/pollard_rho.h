#ifndef POLLARD_RHO_H
#define POLLARD_RHO_H

// Pollard's rho algorithm for factoring small numbers (up to unsigned long long)
// Returns true if a non-trivial factor is found, false otherwise
// factor will contain a factor of N (not necessarily prime)
// This is faster than QS for small numbers (typically < 10^12)
bool pollard_rho(unsigned long long int N, long int& factor, bool debug = false);

#endif // POLLARD_RHO_H
