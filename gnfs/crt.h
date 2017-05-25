#ifndef __CRT_H
#define __CRT_H
#include <vector>
#include <iostream>

// Algorithm 1.3.12 (Inductive Chinese)
template <class Z>
Z crt(const std::vector<Z>& xx, const std::vector<Z>& mm)
{
    int k = xx.size();
    if (k != (int)mm.size())
    {
        std::cout << "Problem: xx.size() = " << xx.size() << ", mm.size = " << mm.size() << std::endl;
    }
    const Z one(1L);

    // Step 1. [Initialize]
    int i = 0;
    Z m = mm[0];
    Z x = xx[0];

    // Step 2. [Finished?]
    while (i < k - 1)
    {
        i++;
        Z u;
        Z v;
        Z m_i = mm[i];
        Z g = extended_gcd(m, m_i, u, v);
        if (g != one)
        {
            std::cout << "Problem: gcd(" << m << "," << m_i << ") = " << g << ", g != 1" << std::endl;
        }
        // Step 3. [Compute next x]
        Z tmp = u;
        tmp *= m;
        tmp *= xx[i];
        x *= v;
        x *= m_i;
        x += tmp;
        m *= m_i;
        x = x % m;
    }
    const Z two(2L);
    if (m - x < m / two) x -= m;
    return x;
}

template <class Z, class Z1>
Z crt(const std::vector<Z>& xx, const std::vector<Z1>& mm)
{
    int k = xx.size();
    if (k != (int)mm.size())
    {
        std::cout << "Problem: xx.size() = " << xx.size() << ", mm.size = " << mm.size() << std::endl;
    }
    const Z one(1L);

    // Step 1. [Initialize]
    int i = 0;
    Z m = mm[0];
    Z x = xx[0];

    // Step 2. [Finished?]
    while (i < k - 1)
    {
        i++;
        Z u;
        Z v;
        Z m_i = mm[i];
        Z g = extended_gcd(m, m_i, u, v);
        //std::cerr << "(m, m_i, u, v) = (" << m << ", " << m_i << ", " << u << ", " << v << "), g = " << g << std::endl;
        if (g != one)
        {
            std::cout << "Problem: gcd(" << m << "," << m_i << ") = " << g << ", g != 1" << std::endl;
        }
        // Step 3. [Compute next x]
        Z new_m = m * m_i;
        Z tmp = u;
        tmp *= m;
        tmp %= new_m;
        tmp *= xx[i];
        tmp %= new_m;
        x *= v;
        x %= new_m;
        x *= m_i;
        x %= new_m;
        x += tmp;
        x %= new_m;
        //std::cerr << "x = " << x << std::endl;
        m = new_m;
        //std::cerr << "m = " << m << std::endl;
        if (x < 0L) x += m;
        //std::cerr << "x = " << x << std::endl;
    }
    const Z two(2L);
    if (m - x < m / two) x -= m;
    // check
    for (size_t i = 0; i < mm.size(); ++i)
    {
        Z p = mm[i];
        Z b = xx[i];
        if ((x - b) % p != 0L)
        {
            std::cerr << "Problem! : p = " << p << ", b = " << b << ", x = " << x << std::endl;
            for (size_t j = 0; j < mm.size(); ++j)
            {
                std::cerr << "mm[" << j << "] = " << mm[j] << ", xx[" << j << "] = " << xx[j] << std::endl;
            }
        }
    }
    return x;
}
#endif
