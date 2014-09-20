#ifndef _CONVERT_H
#define _CONVERT_H
namespace Convert
{
char* parse_FBGNFS_relation(const std::string& str, long long int& a, long long int&b);
void parse_FBGNFS(const std::string& str, long long int& a, long long int& b, char*& alg_str, char*& rat_str);
void extract_FBGNFS_primes(char* str, std::vector<long int>& primes);
};
#endif
