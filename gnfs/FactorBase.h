#ifndef __FACTORBASE_H
#define __FACTORBASE_H
#include "VeryLong.h"
#include "Polynomial.h"
#include "LongModular.h"
#include <vector>
#include <unordered_map>
#include <cstdlib>
#include <string.h>

typedef signed char log_type;
class FactorBase
{
public:
    struct R_p
    {
        friend class FactorBase;
        bool operator<(const R_p& r_p) const
        {
            return (p < r_p.p);
        }
        R_p(int32_t pp) : p(pp)
        {}
        R_p() : p(0)
        {}
        int32_t get_p() const
        {
            return p;
        }
        log_type get_logp() const
        {
            return logp;
        }
//private:
        size_t root_info_index;
        int32_t p;
        log_type logp;
        unsigned char count;
    };
    FactorBase(const Polynomial<VeryLong>& f, long int bound);
    FactorBase(const Polynomial<VeryLong>& f, long int bound, const char* filename = "factorbase.dat");
    FactorBase(const char* filename = "factorbase.dat");
    ~FactorBase();

    typedef const R_p* a_const_iterator;
    typedef R_p* a_iterator;
    typedef const int32_t* a_const_root_iterator;

    a_iterator begin()
    {
        return factor_base_;
    }

    a_iterator end()
    {
        return factor_base_ + factor_base_size_;
    }

    a_const_iterator begin() const
    {
        return factor_base_;
    }

    a_const_iterator end() const
    {
        return factor_base_ + factor_base_size_;
    }

    a_const_root_iterator begin(int32_t p) const;

    a_const_root_iterator end(int32_t p) const;

    a_const_root_iterator begin(a_const_iterator iter) const
    {
        return root_array_ + iter->root_info_index;
    }

    a_const_root_iterator end(a_const_iterator iter) const
    {
        return root_array_ + iter->root_info_index + iter->count;
    }

    bool exists(int32_t p) const
    {
        R_p r_p(p);
        //if (p < B_)
        if (p <= highest_prime())
        {
            return std::binary_search(begin(), end(), r_p);
        }
        else
        {
            return (factor_base_overflow_.find(p) != factor_base_overflow_.end());
        }
    }

    bool exists_extra(int32_t p) const
    {
        return (factor_base_overflow_.find(p) != factor_base_overflow_.end());
    }

    std::vector<int32_t>::const_iterator begin_inert()
    {
        return inert_primes_.begin();
    }

    std::vector<int32_t>::const_iterator end_inert()
    {
        return inert_primes_.end();
    }

    std::vector<int32_t>::const_reverse_iterator rbegin_inert() const
    {
        return inert_primes_.rbegin();
    }

    std::vector<int32_t>::const_reverse_iterator rend_inert() const
    {
        return inert_primes_.rend();
    }

    void add_extra(int32_t p);
    void queue_extra(int32_t p);
    void add_extra_queue();

    void write(const char* filename = "factorbase.dat");
    void read(const char* filename = "factorbase.dat");
    void dump(const char* filename = "factorbase.dmp") const;

    int32_t highest_prime() const
    {
        if (factor_base_size_ <= 0 || factor_base_size_ > factor_base_capacity_) return 0L;
        return factor_base_[factor_base_size_ - 1].p;
    }

    void clear();

private:
    void generate();
    void write_as_we_generate(const char* filename);

    void init()
    {
        const long int initial_allocation = 1800000L;
        factor_base_capacity_ = initial_allocation;
        factor_base_ = reinterpret_cast<R_p*>(std::malloc(factor_base_capacity_ * sizeof(R_p)));
        factor_base_size_ = 0;
        root_array_capacity_ = initial_allocation;
        root_array_ = reinterpret_cast<int32_t*>(std::malloc(root_array_capacity_ * sizeof(int32_t)));
        root_array_size_ = 0;
    }

    void extend()
    {
        const long int incr = 100000L;
        if (factor_base_size_ >= factor_base_capacity_)
        {
            if (factor_base_size_ == 0) init();
            else
            {
                factor_base_capacity_ += incr;
                R_p* new_fb = reinterpret_cast<R_p*>(std::realloc(factor_base_, factor_base_capacity_ * sizeof(R_p)));
                //if (new_fb != factor_base_) std::free(factor_base_);
                factor_base_ = new_fb;
            }
        }
        if (root_array_size_ + 10 >= root_array_capacity_)
        {
            root_array_capacity_ += incr;
            int32_t* new_ra = reinterpret_cast<int32_t*>(std::realloc(root_array_, root_array_capacity_ * sizeof(int32_t)));
            //if (new_ra != root_array_) std::free(root_array_);
            root_array_ = new_ra;
        }
    }

    void shrink()
    {
        if (factor_base_capacity_ > factor_base_size_ + 1)
        {
            factor_base_capacity_ = factor_base_size_ + 1;
            R_p* new_fb = reinterpret_cast<R_p*>(std::realloc(factor_base_, factor_base_capacity_ * sizeof(R_p)));
            //if (new_fb != factor_base_) std::free(factor_base_);
            factor_base_ = new_fb;
        }

        if (root_array_capacity_ > root_array_size_ + 1)
        {
            root_array_capacity_ = root_array_size_ + 1;
            int32_t* new_ra = reinterpret_cast<int32_t*>(std::realloc(root_array_, root_array_capacity_ * sizeof(int32_t)));
            //if (new_ra != root_array_) std::free(root_array_);
            root_array_ = new_ra;
        }
    }

    void add(const char* buffer)
    {
        static char buf[1024];
        extend();
        a_iterator fba_ptr = end();

        strcpy(buf, buffer);
        char* s = strtok(buf, " ");
        int32_t p = std::atol(s);
        fba_ptr->p = p;
        fba_ptr->logp = static_cast<log_type>(log10((double)p) + 0.5);
        if (fba_ptr->logp == 0) fba_ptr->logp = 1;

        unsigned char count = 0;
        size_t index = root_array_size_;
        fba_ptr->root_info_index = index;
        while ((s = strtok(0, " ")))
        {
            int32_t r = std::atol(s);
            root_array_[index] = r;
            ++count;
            ++index;
        }
        fba_ptr->count = count;
        root_array_size_ = index;
        ++factor_base_size_;
    }

public:
    void add_extra(int32_t p, std::vector<LongModular>& roots)
    {
        for (auto& r: roots)
        {
            int32_t root = r.get_long();
            size_t i = 0;
            while (i < factor_base_overflow_[p].size() && factor_base_overflow_[p][i] != root)
            {
                i++;
            }
            if (i >= factor_base_overflow_[p].size())
            {
                factor_base_overflow_[p].push_back(root);
            }
        }
        // Projective roots
        if (f_.coefficient(f_.deg()) % static_cast<long>(p) == 0)
        {
            factor_base_overflow_[p].push_back(p);
        }
    }
private:
    void add(int32_t p, std::vector<LongModular>& roots)
    {
        extend();
        a_iterator fba_ptr = end();
        fba_ptr->p = p;
        fba_ptr->logp = static_cast<log_type>(log10((double)p) + 0.5);
        if (fba_ptr->logp == 0) fba_ptr->logp = 1;

        size_t count = 0;
        size_t index = root_array_size_;
        fba_ptr->root_info_index = index;
        for (size_t j = 0; j < roots.size(); j++)
        {
            int32_t r = roots[j].get_long();
            size_t k = 0;
            while (k < count && root_array_[k + root_array_size_] != r) ++k;
            if (k >= count)
            {
                root_array_[index] = r;
                ++count;
                ++index;
            }
        }
        // Projective roots
        if (f_.coefficient(f_.deg()) % static_cast<long>(p) == 0)
        {
            root_array_[index] = p;
            ++count;
            ++index;
        }
        fba_ptr->count = static_cast<unsigned char>(count);
        root_array_size_ = index;
        ++factor_base_size_;
    }

    typedef std::unordered_map<int32_t, std::vector<int32_t> > fb_overflow;
    fb_overflow factor_base_overflow_;
    std::vector<int32_t> inert_primes_;
    Polynomial<VeryLong> f_;
    long int B_;
    size_t factor_base_size_;
    size_t factor_base_capacity_;
    R_p* factor_base_;
    size_t root_array_size_;
    size_t root_array_capacity_;
    int32_t* root_array_;
    std::vector<int32_t> extra_prime_queue_;
};

#endif
