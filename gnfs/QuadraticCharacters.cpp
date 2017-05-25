#include "QuadraticCharacters.h"
#include "LongModular.h"
#include <vector>
#include "lip.h"

QuadraticCharacters::QuadraticCharacters(const Polynomial<VeryLong>& f, long int B)
    : f_(f), B_(B)
{
    initialise();
}

// This function generates a list of (q, s) pairs satisfying
// q > B, f(s) = 0 mod q, f'(s) != 0 mod q, q prime
void QuadraticCharacters::initialise()
{
    const int max_qc = 100;
    long int q = zpnextb(B_);
    VeryLong qq(q);
    std::vector<LongModular> roots;
    VeryLong c_d = f_.coefficient(f_.deg());
    Polynomial<VeryLong> f1 = f_.derivative();
    const VeryLong zero(0L);
    int count = 0;
    while (count < max_qc)
    {
        if (c_d % qq != zero)
        {
            roots.clear();
            find_roots_mod_p<VeryLong, long int, LongModular>(f_, q, roots);
            for (size_t i = 0; i < roots.size(); i++)
            {
                VeryLong s(roots[i].get_long());
                int repeated_root = (f1.evaluate(s) % qq == zero);
                if (!repeated_root)
                {
                    qs_.push_back(std::pair<VeryLong, VeryLong>(qq, s));
                    count++;
                }
            }
        }
        q = zpnext();
        qq = VeryLong(q);
    }
}

void QuadraticCharacters::generate(const VeryLong& a, const VeryLong& b, std::vector<char>& qcs)
{
    qcs.clear();
    for (size_t i = 0; i < qs_.size(); i++)
    {
        VeryLong q(qs_[i].first);
        VeryLong s(qs_[i].second);
        int x = jacobi_symbol(a - s * b, q);
        //std::cerr << "QCS|" << a << "|" << b << "|" << q << "|" << s << "|" << x << std::endl;
        if (x == 1) qcs.push_back(0);
        else qcs.push_back(1);
    }
}
