// Quadratic Sieve Attempt
#include <vector>
#include <iostream>
#include <algorithm>
#include "gmp.h"
#include "gcd.h"
#include "SparseMatrix.h"
#include "VeryLong.h"
#include "VeryLongModular.h"
#include "timings.h"
#include <limits.h>

extern "C"
{
#include "lip.h"
}

namespace
{
Timing timer("qs.tim");
bool Debug = false;

static long int primes[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999 };

const size_t primes_size = sizeof(primes) / sizeof(long int);

static long int primes_n[primes_size] = { 0, 2, 2, 3, 2, 2, 3, 2, 5, 2, 3, 2, 3, 2, 5, 2, 2, 2, 2, 7, 5, 3, 2, 3, 5, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 2, 2, 5, 2, 2, 2, 7, 5, 2, 3, 2, 3, 2, 2, 3, 7, 7, 2, 3, 5, 2, 3, 2, 3, 2, 2, 2, 11, 5, 2, 2, 5, 2, 2, 3, 7, 3, 2, 2, 5, 2, 2, 3, 7, 2, 2, 7, 5, 3, 2, 3, 5, 2, 3, 2, 13, 3, 2, 2, 5, 2, 3, 2, 2, 2, 2, 2, 3, 2, 5, 2, 3, 7, 7, 3, 2, 3, 2, 3, 3, 2, 5, 2, 2, 2, 5, 2, 2, 2, 2, 2, 11, 3, 2, 2, 5, 3, 2, 3, 7, 2, 2, 2, 3, 2, 2, 3, 2, 2, 11, 2, 3, 2, 5, 2, 3, 2, 5, 2, 7, 3, 3, 5, 2, 2, 3, 3, 2, 3, 5, 3, 2, 11, 2, 2, 2, 7, 5, 3, 3, 2, 2, 3, 2, 3, 2, 2, 3, 5, 2, 2, 2, 11, 13, 5, 2, 2, 2, 2, 3, 11, 2, 3, 5, 2, 3, 2, 7, 2, 2, 3, 2, 3, 2, 5, 2, 3, 2, 13, 7, 3, 3, 5, 2, 2, 3, 3, 3, 2, 2, 3, 7, 3, 2, 2, 2, 3, 3, 2, 5, 7, 2, 2, 11, 2, 2, 3, 2, 3, 17, 3, 2, 2, 5, 2, 3, 5, 7, 2, 2, 2, 2, 2, 5, 3, 2, 2, 2, 3, 2, 2, 3, 2, 2, 2, 2, 5, 3, 5, 3, 2, 2, 11, 2, 5, 3, 5, 2, 2, 7, 5, 2, 3, 3, 2, 2, 3, 2, 2, 2, 3, 2, 2, 2, 5, 2, 3 };

static long int primes_q[primes_size] = { 0, 1, 1, 3, 5, 3, 1, 9, 11, 7, 15, 9, 5, 21, 23, 13, 29, 15, 33, 35, 9, 39, 41, 11, 3, 25, 51, 53, 27, 7, 63, 65, 17, 69, 37, 75, 39, 81, 83, 43, 89, 45, 95, 3, 49, 99, 105, 111, 113, 57, 29, 119, 15, 125, 1, 131, 67, 135, 69, 35, 141, 73, 153, 155, 39, 79, 165, 21, 173, 87, 11, 179, 183, 93, 189, 191, 97, 99, 25, 51, 209, 105, 215, 27, 219, 221, 7, 57, 115, 231, 233, 239, 243, 245, 249, 251, 127, 65, 261, 135, 273, 139, 281, 71, 285, 9, 293, 37, 299, 75, 303, 153, 77, 309, 315, 5, 321, 323, 163, 329, 165, 21, 169, 341, 345, 175, 177, 359, 363, 183, 369, 371, 375, 189, 95, 3, 193, 393, 199, 101, 405, 205, 411, 413, 207, 419, 213, 107, 429, 431, 219, 55, 441, 443, 453, 455, 459, 29, 117, 235, 473, 119, 483, 485, 61, 491, 495, 249, 63, 253, 509, 255, 515, 129, 519, 131, 525, 265, 531, 267, 543, 545, 273, 137, 551, 277, 279, 561, 141, 575, 9, 581, 585, 295, 593, 149, 75, 303, 19, 611, 307, 615, 309, 39, 629, 319, 639, 641, 161, 645, 81, 325, 651, 653, 659, 165, 663, 85, 683, 343, 345, 699, 11, 711, 713, 357, 179, 719, 723, 725, 363, 729, 735, 185, 741, 743, 93, 373, 749, 755, 761, 765, 771, 387, 97, 779, 783, 785, 789, 791, 399, 25, 803, 201, 403, 809, 405, 813, 409, 207, 831, 833, 417, 423, 53, 849, 427, 215, 861, 433, 435, 873, 219, 879, 111, 891, 893, 447, 225, 905, 911, 915, 923, 465, 933, 935, 117, 469, 939, 59, 475, 953, 239, 965, 483, 487, 975, 493, 989, 993, 249, 499, 999 };

static long int primes_e[primes_size] = { 0, 1, 2, 1, 1, 2, 4, 1, 1, 2, 1, 2, 3, 1, 1, 2, 1, 2, 1, 1, 3, 1, 1, 3, 5, 2, 1, 1, 2, 4, 1, 1, 3, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 6, 2, 1, 1, 1, 1, 2, 3, 1, 4, 1, 8, 1, 2, 1, 2, 3, 1, 2, 1, 1, 3, 2, 1, 4, 1, 2, 5, 1, 1, 2, 1, 1, 2, 2, 4, 3, 1, 2, 1, 4, 1, 1, 6, 3, 2, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 1, 2, 1, 3, 1, 6, 1, 4, 1, 3, 1, 2, 3, 1, 1, 7, 1, 1, 2, 1, 2, 5, 2, 1, 1, 2, 2, 1, 1, 2, 1, 1, 1, 2, 3, 8, 2, 1, 2, 3, 1, 2, 1, 1, 2, 1, 2, 3, 1, 1, 2, 4, 1, 1, 1, 1, 1, 5, 3, 2, 1, 3, 1, 1, 4, 1, 1, 2, 4, 2, 1, 2, 1, 3, 1, 3, 1, 2, 1, 2, 1, 1, 2, 3, 1, 2, 2, 1, 3, 1, 7, 1, 1, 2, 1, 3, 4, 2, 6, 1, 2, 1, 2, 5, 1, 2, 1, 1, 3, 1, 4, 2, 1, 1, 1, 3, 1, 4, 1, 2, 2, 1, 7, 1, 1, 2, 3, 1, 1, 1, 2, 1, 1, 3, 1, 1, 4, 2, 1, 1, 1, 1, 1, 2, 4, 1, 1, 1, 1, 1, 2, 6, 1, 3, 2, 1, 2, 1, 2, 3, 1, 1, 2, 2, 5, 1, 2, 3, 1, 2, 2, 1, 3, 1, 4, 1, 1, 2, 3, 1, 1, 1, 1, 2, 1, 1, 4, 2, 1, 5, 2, 1, 3, 1, 2, 2, 1, 2, 1, 1, 3, 2, 1 };

static long int primes_n_q_p[primes_size] =
{ 0, 2, 2, 6, 10, 8, 3, 18, 22, 12, 30, 31, 38, 42, 46, 30, 58, 11, 66, 70, 10, 78, 82, 37, 28, 10, 102, 106, 33, 40, 126, 130, 127, 138, 105, 150, 129, 162, 166, 80, 178, 162, 190, 125, 183, 198, 210, 222, 226, 122, 221, 238, 111, 250, 3, 262, 187, 270, 60, 60, 282, 138, 306, 310, 188, 203, 330, 191, 346, 213, 294, 358, 366, 104, 378, 382, 115, 63, 268, 343, 418, 29, 430, 238, 438, 442, 391, 207, 48, 462, 466, 478, 486, 490, 498, 502, 301, 315, 522, 52, 546, 118, 562, 277, 570, 557, 586, 384, 598, 59, 606, 578, 478, 618, 630, 243, 642, 646, 149, 658, 555, 118, 26, 682, 690, 566, 96, 718, 726, 380, 738, 742, 750, 87, 135, 343, 317, 786, 215, 239, 810, 295, 822, 826, 583, 838, 333, 669, 858, 862, 151, 767, 882, 886, 906, 910, 918, 701, 67, 844, 946, 797, 966, 970, 620, 982, 990, 161, 179, 45, 1018, 374, 1030, 802, 1038, 461, 1050, 958, 1062, 820, 1086, 1090, 530, 486, 1102, 354, 214, 1122, 692, 1150, 1096, 1162, 1170, 243, 1186, 524, 473, 495, 910, 1222, 632, 1230, 691, 672, 1258, 113, 1278, 1282, 792, 1290, 157, 1250, 1302, 1306, 1318, 950, 1326, 574, 1366, 668, 366, 1398, 1022, 1422, 1426, 809, 1091, 1438, 1446, 1450, 497, 1458, 1470, 826, 1482, 1486, 143, 432, 1498, 1510, 1522, 1530, 1542, 88, 251, 1558, 1566, 1570, 1578, 1582, 610, 828, 1606, 979, 1486, 1618, 166, 1626, 316, 239, 1662, 1666, 1449, 92, 69, 1698, 1319, 232, 1722, 1323, 59, 1746, 190, 1758, 912, 1782, 1786, 724, 524, 1810, 1822, 1830, 1846, 61, 1866, 1870, 1496, 137, 1878, 1456, 1683, 1906, 922, 1930, 598, 589, 1950, 259, 1978, 1986, 960, 1585, 1998 };

static long int primes_log[primes_size] =
{ 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };
//{ 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 };

struct Factor_base_item
{
    long int p_;
    long int p2_;
    long int n_;
    long int e_;
    long int q_;
    long int n_q_p_;
    long int log_; // nearest integer to log(p_)
    long int tmem_;
    long int soln1_;
    long int soln2_;
};

size_t factor_base_size = 0;
size_t next_q_index = 0;
bool last_polynomial = false;

Factor_base_item factor_base[primes_size];

mpz_t N_;
//
//                                   2
// coefficients of g   (x) = (ax + b) - N
//                  a,b
//
long int a;
long int b;
long int two_b;
//
//       2
// c = (b - N) / a
//
long long int c;
struct Prime
{
    size_t i_;        // Index into factor_base
    unsigned char e_; // multiplicity
};

struct Relation
{
    void initialize(long int a, long int b, long int x, size_t q_index, bool negative)
    {
        x_ = x;
        a_ = a;
        b_ = b;
        q_index_ = q_index;
        negative_ = negative;
    }

    long long int u()
    {
        long long int u_ = x_;
        u_ *= a_;
        u_ += b_;
        if (u_ < 0) u_ = -u_;
        return u_;
    }

    void display()
    {
        for (size_t i = 0; i < primes_count_; ++i)
        {
            std::cerr << "(" << primes_[i].i_ << ", " << (int)primes_[i].e_ << ")" << std::endl;
        }
    }
    long int a_;
    long int b_;
    long int x_;
    size_t q_index_;
    bool negative_;
    size_t primes_count_;
    Prime primes_[40];
};

size_t smooth_count = 0;
const long int MAX_RELATIONS = 300L;
Relation relations[MAX_RELATIONS];
#if 0
void print_relation(size_t r, Relation* rel)
{
    std::cerr << r << " : ";
    for (size_t j = 0; j < rel->primes_count_; ++j)
    {
        size_t i_ = rel->primes_[j].i_;
        long int e = rel->primes_[j].e_;
        long int p = factor_base[i_].p_;
        std::cerr << "(" << i_ << ", " << p << ", " << e << "), ";
    }
    std::cerr << std::endl;
}
void print_relations()
{
    std::cerr << smooth_count << " relations :" << std::endl;
    for (size_t i = 0; i < smooth_count; ++i)
    {
        print_relation(i, relations + i);
    }
}
#endif

//const long int F = 2000L;
//const long int M = 5000L;
const long int F = 1500L;
const long int M = 9300L;
typedef signed char SIEVE_TYPE;
//typedef int SIEVE_TYPE;

SIEVE_TYPE sieve_array[2*M + 1];
const SIEVE_TYPE* sieve_array_end = sieve_array + sizeof(sieve_array) / sizeof(SIEVE_TYPE);

unsigned long int myzmulmods(unsigned long int a, unsigned long int b, unsigned long int n)
{
    unsigned long int i =  (unsigned long int)( ((unsigned long long)a * (unsigned long long)b) % (unsigned long long)n);
    //if (i < 0) i += n;
    return i;
}

long int find_modular_square_root(long int a, long int p, long int n, long int e, long int q, long int n_q_p)
{
    if (p == 2)
        return a;

    // 1. [Find generator]
    long int z = n_q_p;

    // 2. [Initialize]
    long int y = z;
    long int r = e;
    long int x = zexpmods(a, (q - 1L) / 2L, p);
    long int b = myzmulmods(a, x, p);
    b = myzmulmods(b, x, p);
    x = myzmulmods(a, x, p);

    while (true)
    {
        // 3. [Find exponent]
        if (b == 1L)
        {
#if 0
            // check
            long int check = (x * x - a) % p;
            if (check != 0)
            {
                std::cerr << "Problem : find_modular_square_root : x = " << x << " is not a modular square root of a = " << a << std::endl;
            }
#endif
            return x;
        }

        // Find the smallest m >= 1 such that b^2^m = 1 (mod p)
        long int m = 1;
        long int b2m = myzmulmods(b, b, p);
        while (b2m != 1L)
        {
            b2m = myzmulmods(b2m, b2m, p);
            ++m;
        }

        if (m == r)
        {
            return 0L;
        }

        // 4. [Reduce exponent]
        long int rm1 = r - m - 1;
        long int power = 1 << rm1;
        long int t = zexpmods(y, power, p);
        y = myzmulmods(t, t, p);
        r = m;
        x = myzmulmods(x, t, p);
        b = myzmulmods(b, y, p);
    }
}
#if 0
void print_factor_base()
{
    Factor_base_item* fbi = factor_base;
    for (size_t i = 0; i < factor_base_size; ++i)
    {
        std::cerr << "i : " << fbi->p_ << ", " << fbi->log_ << ", " << fbi->n_ << ", " << fbi->e_ << ", " << fbi->q_ << ", " << fbi->n_q_p_ << ", " << fbi->tmem_ << std::endl;
        ++fbi;
    }
}
#endif
void compute_startup_data(unsigned long long int N)
{
    if (Debug) std::cerr << "Computing start up data ... ";
    unsigned long int q = (unsigned long int)(N >> 32);
    unsigned long int r = (unsigned long int)(N & 0x00000000ffffffff);
    mpz_init_set_ui(N_, q);
    mpz_mul_2exp(N_, N_, 32);
    mpz_add_ui(N_, N_, r);

    factor_base_size = 0;
    for (size_t i = 0; i < primes_size; ++i)
    {
        long int p = primes[i];
        if (p >= F)
            break;
        if (p == 2 || mpz_kronecker_si(N_, p) == 1)
        {
            factor_base[factor_base_size].p_ = p;
            factor_base[factor_base_size].p2_ = p * p;
            factor_base[factor_base_size].log_ = primes_log[i];
            factor_base[factor_base_size].n_ = primes_n[i];
            factor_base[factor_base_size].e_ = primes_e[i];
            factor_base[factor_base_size].q_ = primes_q[i];
            factor_base[factor_base_size].n_q_p_ = primes_n_q_p[i];
            factor_base[factor_base_size].tmem_ =
                find_modular_square_root(N % p, p, primes_n[i], primes_e[i], primes_q[i], primes_n_q_p[i]);
            ++factor_base_size;
        }
    }
    if (Debug) std::cerr << "done" << std::endl;
    //print_factor_base();
    smooth_count = 0;
    next_q_index = factor_base_size - 1;
}

void initialization_stage(unsigned long long int N)
{
    if (Debug) std::cerr << "Initializing ... ";
    // Find a prime q ~ sqrt(sqrt(2N)/M) such that N is a quadratic residue mod q

    if (last_polynomial)
    {
        throw "Run out of polynomials";
    }

    if (next_q_index == 0)
    {
        last_polynomial = true;
    }
    size_t i = next_q_index--;
    Factor_base_item* fbi = factor_base + i;
    long int q = factor_base[i].p_;
    a = q * q;
    // find modular square root of N mod q and then lift via Hensel's lemma
    // to get a modular square root mod q^2
    b = find_modular_square_root(N % q, q, fbi->n_, fbi->e_, fbi->q_, fbi->n_q_p_);
#if 0
    unsigned long long int b_ll = b;
    long int r_ = ((b_ll * b_ll - N) / q) % q;
#else
    VeryLong b_vl(b);
    VeryLong N_vl(N);
    long int r_ = ((b_vl * b_vl - N_vl) / q) % q;
#endif
    long int s = (::inverse(2*b, q) * r_) % q;
    b -= s * q;
    b %= a;
    two_b = b + b;
    // Check
#if 0
    b_ll = b;
#else
    b_vl = b;
#endif
#if 0
    if ((b_ll * b_ll - N) % a != 0)
    {
        std::cerr << "Problem: b = " << b << " is not a modular square root of " << N << "mod " << a << std::endl;
    }
#else
    if ((b_vl * b_vl - N_vl) % a != 0)
    {
        std::cerr << "Problem: b = " << b << " is not a modular square root of " << N << "mod " << a << std::endl;
    }
#endif
#if 0
    c = (b_ll * b_ll - N) / a;
#else
    VeryLong c_vl = (b_vl * b_vl - N_vl) / a;
    c = c_vl.get_long_long();
#endif
    //std::cerr << "a = " << a << ", b = " << b << ", c = " << c << std::endl;

    //                                        2
    // compute solutions to g   (x) = (ax + b) - N = 0 (mod p)
    //                       a,b
    // for each p in the factor base
    //
    for (size_t i = 0; i < factor_base_size; ++i)
    {
        long int log = factor_base[i].log_;
        if (log == 0)
            continue;
        long int p = factor_base[i].p_;
        long int tmem = factor_base[i].tmem_;
        long int a_inv = ::inverse(a, p);
        factor_base[i].soln1_ = (a_inv * (tmem - b)) % p + M;
        factor_base[i].soln2_ = (a_inv * (-tmem - b)) % p + M;
    }

    if (Debug) std::cerr << "done" << std::endl;
}

//#define USE_CACHE 1
#ifdef USE_CACHE
struct SieveCache
{
    SIEVE_TYPE* ptr_;
    long int log_;
    bool operator<(const SieveCache& sc) const
    {
        return (ptr_ < sc.ptr_);
    }
};
#endif

void sieve_stage()
{
    if (Debug) std::cerr << "Sieving ... ";
    ::memset(sieve_array, 0, sizeof(sieve_array));
    size_t i = 0;
#ifdef USE_CACHE
    const size_t sieve_cache_size = 16L;
    size_t cache_index = 0;
    SieveCache cache[sieve_cache_size];
#endif

    for (; i < factor_base_size; ++i)
    {
        long int log = factor_base[i].log_;
        if (log == 0)
            continue;
        long int p = factor_base[i].p_;
        long int soln1 = factor_base[i].soln1_;
        long int soln2 = factor_base[i].soln2_;

        long int k_min = soln1 / p;
        k_min = -k_min;
        SIEVE_TYPE* ptr = sieve_array + soln1 + k_min * p;
        SIEVE_TYPE* end = sieve_array + M + M + 1;
        while (ptr < end)
        {
#ifndef USE_CACHE
            *ptr += log;
#else
            cache[cache_index].ptr_ = ptr;
            cache[cache_index].log_ = log;
            ++cache_index;
            if (cache_index == sieve_cache_size)
            {
                std::sort(cache, cache + cache_index);
                for (size_t i = 0; i < cache_index; ++i)
                {
                    SIEVE_TYPE* ptr1 = cache[i].ptr_;
                    long int log1 = cache[i].log_;
                    *ptr1 += log1;
                }
                cache_index = 0;
            }
#endif
            ptr += p;
        }

        k_min = soln2 / p;
        k_min = -k_min;
        ptr = sieve_array + soln2 + k_min * p;
        while (ptr < end)
        {
#ifndef USE_CACHE
            *ptr += log;
#else
            cache[cache_index].ptr_ = ptr;
            cache[cache_index].log_ = log;
            ++cache_index;
            if (cache_index == sieve_cache_size)
            {
                std::sort(cache, cache + cache_index);
                for (size_t i = 0; i < cache_index; ++i)
                {
                    SIEVE_TYPE* ptr1 = cache[i].ptr_;
                    long int log1 = cache[i].log_;
                    *ptr1 += log1;
                }
                cache_index = 0;
            }
#endif
            ptr += p;
        }
    }
#ifdef USE_CACHE
    if (cache_index)
    {
        std::sort(cache, cache + cache_index);
        for (size_t i = 0; i < cache_index; ++i)
        {
            SIEVE_TYPE* ptr1 = cache[i].ptr_;
            long int log1 = cache[i].log_;
            *ptr1 += log1;
        }
    }
#endif
    if (Debug) std::cerr << "done" << std::endl;
}

long long int g(long int x)
{
    //
    //                   2                      2
    // g   (x) = (ax + b) - N = 0 (mod p) = a(ax + 2bx + c)
    //  a,b
    //
    //             2
    // where c = (b - N) / a
    //
    // (ax + 2b)x + c
    long long int result = x;
    result *= a;
    result += two_b;
    result *= x;
    result += c;
    return result;
}

/*inline*/ void check_for_smoothness(long int x)
{
    // trial divide g(x)
    long long int value = g(x);
    //long int value_l = 0;
    unsigned long int value_l = 0;

    int sign = 1;
    if (value < 0)
    {
        value = -value;
        sign = -1;
    }

    Relation* relation = relations + smooth_count;
    relation->primes_count_ = 0;
    long int max_p = factor_base[factor_base_size - 1].p_;

    size_t i = 0;
    long int p;
    for (; (unsigned long long int)value > ULONG_MAX && i < factor_base_size; ++i)
    {
        unsigned char multiplicity = 0;
        p = factor_base[i].p_;
        while (value % p == 0)
        {
            value /= p;
            ++multiplicity;
        }
        if (multiplicity > 0)
        {
            relation->primes_[relation->primes_count_].i_ = i;
            relation->primes_[relation->primes_count_].e_ = multiplicity;
            relation->primes_count_++;
        }
#if 1
        if (value < factor_base[i].p2_)
        {
            //std::cerr << "value = " << value << ", value / p = " << value / p << ", p = " << p << std::endl;
            return;
        }
#endif
    }

    value_l = value;
    for (; value_l > 1L && i < factor_base_size; ++i)
    {
        p = factor_base[i].p_;
        unsigned char multiplicity = 0;
        while (value_l % p == 0)
        {
            value_l /= p;
            ++multiplicity;
        }
        if (multiplicity > 0)
        {
            relation->primes_[relation->primes_count_].i_ = i;
            relation->primes_[relation->primes_count_].e_ = multiplicity;
            relation->primes_count_++;
        }
#if 1
        if (value_l > (unsigned long int)max_p && value_l < (unsigned long int)factor_base[i].p2_)
        {
            //std::cerr << "value_l = " << value_l << ", value_l / p = " << value_l / p << ", p = " << p << std::endl;
            return;
        }
#endif
    }
    if (value_l == 1UL)
    {
        // we have a smooth value
        relation->initialize(a, b, x, next_q_index + 1, (sign == -1));
        ++smooth_count;
    }
}

bool trial_division_stage(unsigned long long int N)
{
    if (Debug) std::cerr << "Checking for smooth values ... ";
    // Scan sieve array for locations x that have accumulated a value of at least log(M sqrt(N)) minus
    // a small error term
    const long int fudge = 1;
    long int cutoff = static_cast<long int>(log10(static_cast<double>(N))/2 + log10(static_cast<double>(M))) - fudge;

    const size_t excess = 5;
    long int x = -M;

    for (SIEVE_TYPE* ptr = sieve_array; ptr != sieve_array_end; ++ptr, ++x)
    {
        if (*ptr > cutoff)
        {
            check_for_smoothness(x);
            if (smooth_count >= factor_base_size + excess)
                break;
        }
    }

    //std::cerr << smooth_count << " smooth relations found with a = " << a << ", b = " << b << ", c = " << c << std::endl;
    //std::cerr << "Factor base size = " << factor_base_size << std::endl;
    //std::cerr << "Excess = " << (int)smooth_count - (int)factor_base_size << std::endl;
    //std::cerr << "done" << std::endl;
    return (smooth_count >= factor_base_size + excess);
}

bool linear_algebra_stage(unsigned long long int N, long int& factor)
{
    bool factor_found = false;
    timer.start("linear_algebra_stage 1 - remove singletons");
    if (Debug) std::cerr << "Solving matrix ... ";
    // Build a set of BitMatrix objects to represent the matrix.
    // Each row corresponds to a prime, and each column to a relation
    // Remove singletons first, i.e. relations which contain a prime
    // which only appear in one relation
    Relation* rel;

    std::vector<int> relation_removed(smooth_count, 0);
    std::vector<int> relations_containing_prime(factor_base_size, 0);
    size_t column_count = smooth_count;

    bool singletons = true;
    while (singletons)
    {
        relations_containing_prime.resize(factor_base_size, 0);
        // find out how many unremoved relations each prime currently exists in
        rel = relations;
        for (size_t i = 0; i < smooth_count; ++i, ++rel)
        {
            if (!relation_removed[i])
            {
                for (size_t j = 0; j < rel->primes_count_; ++j)
                {
                    if (rel->primes_[j].e_ % 2)
                    {
                        relations_containing_prime[rel->primes_[j].i_]++;
                    }
                }
            }
        }

        // detect singleton relations, i.e. relations which contain a prime
        // that is in no other relation
        singletons = false;
        rel = relations;
        for (size_t i = 0; i < smooth_count; ++i, ++rel)
        {
            if (!relation_removed[i])
            {
                for (size_t j = 0; j < rel->primes_count_; ++j)
                {
                    if (rel->primes_[j].e_ % 2)
                    {
                        if (relations_containing_prime[rel->primes_[j].i_] == 1)
                        {
                            // this is a singleton relation, so mark for removal
                            relation_removed[i] = 1;
                            --column_count;
                            singletons = true;
                            break; // escape from loop round primes in relation
                        }
                    }
                }
            }
        }
    }

    std::vector<size_t> inverse_relation_index;
    for (size_t i = 0; i < smooth_count; ++i)
    {
        if (!relation_removed[i])
        {
            inverse_relation_index.push_back(i);
        }
    }

    // set up prime_index which gives
    // the mapping between primes in the factor base and rows of the
    // final matrix. If prime_index for a prime is -1, then that prime
    // does not correspond to a row in the final matrix
    std::vector<long int> prime_index(factor_base_size, -1);
    size_t next_prime_index = 0;
    size_t row_count = factor_base_size + 1;
    for (size_t i = 0; i < relations_containing_prime.size(); ++i)
    {
        if (relations_containing_prime[i] <= 1)
        {
            // only one relation contains this prime, so mark this prime
            // as not corresponding to a row in the final matrix ...
            prime_index[i] = -1;
            // ... and hence reduce the number of rows in the final matrix ...
            --row_count;
        }
        else
        {
            // this prime will appear in the final matrix as row next_prime_index
            prime_index[i] = next_prime_index;
            ++next_prime_index;
        }
    }

    if (Debug)
    {
        std::cerr << "row_count = " << row_count << std::endl;
        std::cerr << "column_count = " << column_count << std::endl;
    }
    timer.stop();

    timer.start("linear_algebra_stage 2 - build matrix");

    size_t N1 = BitOperations::BITS_IN_WORD;
    int bm_count = column_count / N1;
    if (column_count % N1 != 0)
    {
        ++bm_count;
    }
    std::vector<BitMatrix> M(bm_count);
    for (auto& bmrow: M)
    {
        bmrow.row_.resize(row_count);
        for (size_t i = 0; i < bmrow.row_.size(); ++i) bmrow.row_[i] = 0UL;
        bmrow.cols_ = N1;
    }
    if (column_count % N1 != 0)
    {
        M[bm_count - 1].cols_ = column_count % N1;
    }

    rel = relations;
    size_t next_relation_index = 0;
    int u = 0;
    for (size_t i = 0; i < smooth_count; ++i, ++rel)
    {
        if (!relation_removed[i])
        {
            BitMatrix& bm = M[u];
            for (size_t j = 0; j < rel->primes_count_; ++j)
            {
                if (rel->primes_[j].e_ % 2 && prime_index[rel->primes_[j].i_] >= 0)
                {
                    BitOperations::setBit(next_relation_index, bm.row_[prime_index[rel->primes_[j].i_]]);
                }
            }
            if (rel->negative_)
            {
                BitOperations::setBit(next_relation_index, bm.row_[row_count - 1]);
            }
            ++next_relation_index;
            if (next_relation_index >= N1)
            {
                ++u;
                next_relation_index = 0;
            }
        }
    }
    timer.stop();

    timer.start("linear_algebra_stage 3 - kernel");
    BitMatrix kerM;

    ::kernel(M, kerM);
    timer.stop();
    //print_relations();

    timer.start("linear_algebra_stage 4 - find dependency");
    if (Debug)
    {
        std::cerr << "kerM is " << kerM.rows() << " x " << kerM.cols() << std::endl;
        std::cerr << kerM << std::endl;
    }
    if (Debug) std::cerr << "done" << std::endl;

    VeryLongModular::set_default_modulus(N);

    // Now check the dependencies
    BitMatrix& bm = kerM;
    for (size_t c = 0; !factor_found && c < bm.cols(); ++c)
    {
        // each column corresponds to a dependency, i.e. picks out the relations we need to
        // multiply together to get a square
        VeryLongModular x_(1L);
        std::vector<int> prime_count(factor_base_size, 0);

        for (size_t r = 0; r < bm.rows(); ++r)
        {
            if (BitOperations::bitSet(c, bm.row_[r]))
            {
                //Relation* rel = relations + r;
                Relation* rel = relations + inverse_relation_index[r];
                //print_relation(r, rel);
                x_ *= rel->u();

                for (size_t j = 0; j < rel->primes_count_; ++j)
                {
                    prime_count[rel->primes_[j].i_] += rel->primes_[j].e_;
                }
                prime_count[rel->q_index_] += 2; // take account of the extra factor of a = q^2
            }
        }
        VeryLongModular y_(1L);
        for (size_t j = 0; j < factor_base_size; ++j)
        {
            if (prime_count[j] % 2)
            {
                std::stringstream ss;
                ss << "multiplicity of p = " << factor_base[j].p_ << " (" << prime_count[j] << ") is not even in dependency " << c << std::endl;
                std::cerr << ss.str() << std::endl;
                throw "multiplicity of p is not even in dependency";
            }
            else
            {
                int pc = prime_count[j] / 2;
                long int p = factor_base[j].p_;
                for (int ii = 0; ii < pc; ++ii)
                {
                    y_ *= p;
                }
            }
        }

        long long int x = x_.get_very_long().get_long_long();
        long long int y = y_.get_very_long().get_long_long();

        long long int w = x - y;
        if (w < 0) w = -w;
        unsigned long long int f = ::gcd((unsigned long long int)w, N);
        if (f > 1 && f < N)
        {
            long long int g = N / f;
            factor = f;
            if (g < (long long int)f)
                factor = g;
            if (Debug)
            {
                std::cerr << "N = " << N << std::endl;
                std::cerr << "f = " << f << std::endl;
                std::cerr << "g = N / f = " << g << std::endl;
                std::cerr << "N - fg = " << N - f * g << std::endl;
            }
            factor_found = true;
            timer.stop();
        }
    }
    timer.stop();
    return factor_found;
}

};

bool QS(unsigned long long int N, long int& factor, bool debug)
{
    try
    {
        Debug = debug;
        timer.start("compute_startup_data");
        compute_startup_data(N);
        timer.stop();

        bool sieving_done = false;
        while (!sieving_done)
        {
            timer.start("initialization_stage");
            initialization_stage(N);
            timer.stop();

            timer.start("sieve_stage");
            sieve_stage();
            timer.stop();

            timer.start("trial_division_stage");
            sieving_done = trial_division_stage(N);
            timer.stop();
        }

        bool done = false;
        //timer.start("linear_algebra_stage");
        //for (int i = 0; i < 10; ++i)
        {
            done = linear_algebra_stage(N, factor);
        }
        //timer.stop();

        //timer.summary();
        mpz_clear(N_);
        return done;
    }
    catch (const char* m)
    {
        std::cerr << "Problem: " << m << std::endl;
        mpz_clear(N_);
        return false;
    }
}

//#define MAINQS
#ifdef MAINQS
int main(int argc, char* argv[])
{
    if (argc > 1) Debug = true;
    long int factor;
    //long long int N = 485320736930657LL;
    //unsigned long long int N = 181442377297009ULL;
    //unsigned long long int N = 17562016984893972197ULL;
    //unsigned long long int N = 5957934429021177763ULL;
    unsigned long long int N = 140666453917430783ULL;
    //long long int N = 8862154124940029LL;
    //long long int N = 1141154556981208597LL;
    //long long int N = 34061096999LL;
    //long long int N = 8157068382101LL;
    QS(N, factor, Debug);
    return 0;
}
#endif
