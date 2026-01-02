#ifndef LATTICE_SIEVER_H
#define LATTICE_SIEVER_H
#include "VeryLong.h"
#include "Polynomial.h"
#include "FactorBase.h"
#include "SieveConfig.h"
#include <utility>
#include <string>
#include <sstream>
#include <exception>
#include <vector>
#include "PointerHashTable.h"
#include "Parallelogram.h"
#include "BitOperations.h"
#include "timings.h"
//#define RESIEVE1 1
//#define OVERFLOW_CHECK 1
class LatticeSiever
{
public:
    friend class LatticeSieverTest;
    LatticeSiever(const std::string& config_file = "sieve.cfg");
    ~LatticeSiever();
    bool sieve(long int q);

    typedef signed char SIEVE_TYPE;
private:
    typedef std::vector<long int> FactorList;

    struct PrimeFactor
    {
        uint32_t offset_;
        int32_t p_;
        PrimeFactor() : offset_(0), p_(0)
        {}
        bool operator< (const PrimeFactor& t) const
        {
            if (offset_ < t.offset_) return true;
            return false;
        }
    };
#if 0
    class IPrimeFactorList
    {
    public:
        virtual void add(LatticeSiever::SIEVE_TYPE* ptr, long int p) = 0;
    };

    template <int bucket_bits, int sieve_array_size>
    class BucketedPrimeFactorList : public IPrimeFactorList
    {
    public:
        enum
        {
            bucket_count = sieve_array_size >> bucket_bits
        };
        BucketedPrimeFactorList(LatticeSiever::SIEVE_TYPE* const sieve_array, long int list_size, bool debug = false)
            : sieve_array_(sieve_array), debug_(debug)
        {
            reset();
        }
        BucketedPrimeFactorList()
            : debug_(false)
        {}
        ~BucketedPrimeFactorList()
        {
        }
        void reset()
        {
            for (size_t i = 0; i < bucket_count; ++i)
            {
                buckets_[i].clear();
            }
        }
        void add(LatticeSiever::SIEVE_TYPE* ptr, long int p)
        {
            std::vector<PrimeFactor>& pfv = buckets_[(ptr - sieve_array_) >> bucket_bits];
            pfv.push_back(PrimeFactor(ptr, p));
        }
        void set_end()
        {
        }

        struct const_iterator
        {
            const_iterator(const BucketedPrimeFactorList<bucket_bits, sieve_array_size>& bpfl, size_t bucket, size_t index)
                : bpfl_(bpfl), bucket_(bucket), index_(index)
            {
                while (bucket_ < bucket_count &&
                        index_ >= bpfl_.buckets_[bucket_].size())
                {
                    ++bucket_;
                    index_ = 0;
                }
            }

            const PrimeFactor* operator*() const
            {
                if (bucket_ < bucket_count)
                {
                    return &bpfl_.buckets_[bucket_][index_];
                }
                return 0;
            }
            const_iterator& operator++()
            {
                ++index_;
                while (bucket_ < bucket_count &&
                        index_ >= bpfl_.buckets_[bucket_].size())
                {
                    ++bucket_;
                    index_ = 0;
                }
                return *this;
            }

            bool operator!=(const const_iterator& rhs) const
            {
                return (bucket_ != rhs.bucket_ || index_ != rhs.index_);
            }
            const BucketedPrimeFactorList<bucket_bits, sieve_array_size>& bpfl_;
            size_t bucket_;
            size_t index_;
        };

        const_iterator begin() const
        {
            return const_iterator(*this, 0, 0);
        }
        const_iterator end() const
        {
            return const_iterator(*this, bucket_count, 0);
        }

    private:
        LatticeSiever::SIEVE_TYPE* const sieve_array_;
        std::vector<PrimeFactor> buckets_[bucket_count];
        bool debug_;
    };
#endif

    class PrimeFactorList /*: public IPrimeFactorList*/
    {
    public:
        PrimeFactorList(long int list_size, LatticeSiever::SIEVE_TYPE* const sieve_array, bool debug = false)
            : prime_factor_(0), pf_ptr_(0), pf_end_(0), list_size_(list_size), sieve_array_(sieve_array), debug_(debug)
        {
            prime_factor_ = new PrimeFactor [ list_size ];
            reset();
        }
        ~PrimeFactorList()
        {
            delete [] prime_factor_;
        }
        void reset()
        {
            pf_ptr_ = prime_factor_;
            pf_end_ = prime_factor_ + list_size_;
        }
        void add(LatticeSiever::SIEVE_TYPE* ptr, int32_t p)
        {
            //pf_ptr_->set(ptr, p);
            //pf_ptr_->ptr_ = ptr;
            pf_ptr_->offset_ = ptr - sieve_array_;
            pf_ptr_->p_ = p;
            ++pf_ptr_;
#ifdef OVERFLOW_CHECK
            if (pf_ptr_ >= pf_end_)
            {
                std::stringstream ss;
                ss << "PrimeFactorList::add(" << ptr << ", " << p << "), pf_list_ overflowed";
                throw std::out_of_range(ss.str());
            }
#endif
        }
        void set_end()
        {
            pf_end_ = pf_ptr_;
            size_t pf_list_size = pf_end_ - prime_factor_;
            if (debug_)
            {
                std::cerr << "Size of pf_list_ = " << pf_list_size << std::endl;
            }
        }
        void sort()
        {
            std::sort(prime_factor_, pf_end_);
        }
        PrimeFactor* begin() const
        {
            return prime_factor_;
        }
        PrimeFactor* end() const
        {
            return pf_end_;
        }
        //private:
        PrimeFactor* prime_factor_;
        PrimeFactor* pf_ptr_;
        PrimeFactor* pf_end_;
        long int list_size_;
        LatticeSiever::SIEVE_TYPE* const sieve_array_;
        bool debug_;
    };

    typedef uint32_t bucket_size_t;
    struct SieveCacheItem
    {
        bucket_size_t offset_;
        int32_t p_;
        signed char logp_;

        SieveCacheItem()
            : offset_(0), p_(0), logp_(0) {}
        static void set_pf_list(LatticeSiever::PrimeFactorList* pf_list)
        {
            pf_list_ = pf_list;
        }
        static LatticeSiever::PrimeFactorList* pf_list_;
    };

    template <int cache_size>
    struct SieveCacheBucket
    {
        SieveCacheItem* next_cache_;
        LatticeSiever::SIEVE_TYPE* base_;
        SieveCacheItem cache_[cache_size];

        SieveCacheBucket() : next_cache_(cache_) {}

        void init(LatticeSiever::SIEVE_TYPE* base)
        {
            base_ = base;
            next_cache_ = cache_;
        }
    };
//#define BUCKET_BITS 1

#ifdef BUCKET_BITS
    template <int cache_size, int bucket_bits, int sieve_array_size>
#else
    template <int cache_size, uint32_t bucket_size, int sieve_array_size>
#endif
    class SieveCache
    {
//#define DEBUG_SIEVE_CACHE 1
    public:
#ifdef BUCKET_BITS
        static const size_t bucket_count = sieve_array_size >> bucket_bits;
        static const uint32_t bucket_size = 1 << bucket_bits;
#else
        /*
            sieve_array_size                    bucket_count            (sieve_array_size - 1)/ bucket_size
            1 - bucket_size                                1                                              0
            (bucket_size + 1) - 2*bucket_size              2                                              1
            (2*bucket_size + 1) - 3*bucket_size            3                                              2
            ...
            ((k-1)*bucket_size + 1) - k*bucket_size        k                                            k-1
            ...


         */
        static const size_t bucket_count = ((sieve_array_size - 1) / bucket_size) + 1;
        // Precompute reciprocal for fast division by bucket_size
        // Using: (offset * reciprocal) >> 32 ≈ offset / bucket_size
        static constexpr uint64_t bucket_size_reciprocal = 
            (uint64_t)((1ULL << 32) + bucket_size - 1) / bucket_size;
#endif
        //SieveCache(LatticeSiever::SIEVE_TYPE* const sieve_array, const BitArray64<sieve_array_size>& sieve_bit_array)
        SieveCache(SIEVE_TYPE* const sieve_array, const BitArray64<sieve_array_size>& sieve_bit_array)
            : sieve_array_(sieve_array), sieve_bit_array_(sieve_bit_array)
        {
            non_empty_buckets_.reserve(bucket_count / 4);  // Reserve 25% capacity
            LatticeSiever::SIEVE_TYPE* base = sieve_array_;
            for (size_t i = 0; i < bucket_count; ++i)
            {
                buckets_[i].init(base);
                base += bucket_size;
            }
#ifdef DEBUG_SIEVE_CACHE
            std::ostringstream debug_file_name;
            debug_file_name << "SieveCache_" << (long int)(this) << ".dbg";
            debug_file_.open(debug_file_name.str().c_str());
#endif
        }
        ~SieveCache()
        {
#ifdef DEBUG_SIEVE_CACHE
            debug_file_.close();
#endif
        }

        void add(uint32_t offset, int32_t count, int32_t inc, FactorBase::a_iterator iter)
        {
            SIEVE_TYPE* __restrict__ sieve_array = sieve_array_;
            while (count)
            {
                --count;
                offset += inc;
#ifdef BUCKET_BITS
                size_t bucket_idx = offset >> bucket_bits;
                SieveCacheBucket<cache_size>& scb = buckets_[bucket_idx];
#else
                size_t bucket_idx = (uint32_t)(((uint64_t)offset * bucket_size_reciprocal) >> 32);
                SieveCacheBucket<cache_size>& scb = buckets_[bucket_idx];
#endif
                SieveCacheItem* const & item = scb.next_cache_;
                
                // Track if this bucket was previously empty
                if (item == scb.cache_)
                {
                    non_empty_buckets_.push_back(bucket_idx);
                }
                
                item->offset_ = offset;
                item->p_ = iter->p;
                item->logp_ = iter->logp;
                if (item == scb.cache_ + cache_size - 1)
                {
                    const SieveCacheItem* it = scb.cache_;
                    for (long int i = 0; i < cache_size; ++i, ++it)
                    {
                        if (!sieve_bit_array_.isSet(it->offset_))
                        {
#ifdef DEBUG_SIEVE_CACHE
                            debug_file_ << std::hex << size_t(it->offset_) << std::dec << std::endl;
#endif
                            *(it->offset_ + sieve_array) += it->logp_;
                            
                            // Cache the pointer for faster access
                            PrimeFactor* pf = SieveCacheItem::pf_list_->pf_ptr_;
                            pf->offset_ = it->offset_;
                            pf->p_ = it->p_;
                            SieveCacheItem::pf_list_->pf_ptr_ = pf + 1;
                        }
                    }
                    scb.next_cache_ = scb.cache_;
                }
                else
                {
                    ++scb.next_cache_;
                }
            }
        }

#ifdef RESIEVE1
        void add1(uint32_t offset, int32_t count, int32_t inc, FactorBase::a_iterator iter)
        {
            while (count)
            {
                --count;
                offset += inc;
                SieveCacheBucket<cache_size>& scb = buckets_[offset >> bucket_bits];
                SieveCacheItem* const & item = scb.next_cache_;
                item->offset_ = offset;
                item->p_ = iter->p;
                item->logp_ = iter->logp;
                if (item == scb.cache_ + cache_size - 1)
                {
                    SieveCacheItem* it = scb.cache_;
                    for (long int i = 0; i < cache_size; ++i, ++it)
                    {
                        if (!sieve_bit_array_.isSet(it->offset_))
                        {
                            *(it->offset_ + sieve_array_) += it->logp_;
                        }
                    }
                    scb.next_cache_ = scb.cache_;
                }
                else
                {
                    ++scb.next_cache_;
                }
            }
        }

        void add1_again(uint32_t offset, int32_t count, int32_t inc, FactorBase::a_iterator iter)
        {
            while (count)
            {
                --count;
                offset += inc;
#if 1
                SieveCacheBucket<cache_size>& scb = buckets_[offset >> bucket_bits];
                SieveCacheItem* const & item = scb.next_cache_;
                item->offset_ = offset;
                item->p_ = iter->p;
                item->logp_ = 0;
                if (item == scb.cache_ + cache_size - 1)
                {
                    SieveCacheItem* it = scb.cache_;
                    for (long int i = 0; i < cache_size; ++i, ++it)
                    {
                        if (!sieve_bit_array_.isSet(it->offset_))
                        {
                            SieveCacheItem::pf_list_->pf_ptr_->offset_ = it->offset_;
                            SieveCacheItem::pf_list_->pf_ptr_->p_ = it->p_;
                            ++SieveCacheItem::pf_list_->pf_ptr_;
                        }
                    }
                    scb.next_cache_ = scb.cache_;
                }
                else
                {
                    ++scb.next_cache_;
                }
#else
                if (!sieve_bit_array_.isSet(offset))
                {
                    SieveCacheItem::pf_list_->pf_ptr_->offset_ = offset;
                    SieveCacheItem::pf_list_->pf_ptr_->p_ = iter->p;
                    ++SieveCacheItem::pf_list_->pf_ptr_;
                }
#endif
            }
        }
#endif
        void dump(bool add_to_pf_list = true)
        {
            for (size_t bucket_index : non_empty_buckets_)
            {
                SieveCacheBucket<cache_size>& scb = buckets_[bucket_index];
                SieveCacheItem* it = scb.cache_;
                while (it != scb.next_cache_)
                {
                    if (!sieve_bit_array_.isSet(it->offset_))
                    {
                        *(it->offset_ + sieve_array_) += it->logp_;
                        if (add_to_pf_list)
                        {
                            SieveCacheItem::pf_list_->add(it->offset_ + sieve_array_, it->p_);
                        }
                    }
                    ++it;
                }
                scb.next_cache_ = scb.cache_;
            }
            non_empty_buckets_.clear();
        }
        
        // Block-based dump: only process entries in the specified offset range
        // This enables cache blocking during sieving
        void dump_block(size_t block_start, size_t block_end, bool add_to_pf_list = true)
        {
            // Process only non-empty buckets, but filter by block range
            std::vector<size_t> buckets_to_keep;
            buckets_to_keep.reserve(non_empty_buckets_.size());
            
            for (size_t bucket_index : non_empty_buckets_)
            {
                SieveCacheBucket<cache_size>& scb = buckets_[bucket_index];
                SieveCacheItem* read_ptr = scb.cache_;
                SieveCacheItem* write_ptr = scb.cache_;
                bool bucket_still_has_items = false;
                
                while (read_ptr != scb.next_cache_)
                {
                    // Check if this item is in the current block
                    if (read_ptr->offset_ >= block_start && read_ptr->offset_ < block_end)
                    {
                        // Process this item
                        if (!sieve_bit_array_.isSet(read_ptr->offset_))
                        {
                            *(read_ptr->offset_ + sieve_array_) += read_ptr->logp_;
                            if (add_to_pf_list)
                            {
                                SieveCacheItem::pf_list_->add(read_ptr->offset_ + sieve_array_, read_ptr->p_);
                            }
                        }
                        // Don't keep this item
                    }
                    else
                    {
                        // Keep this item for later blocks
                        if (write_ptr != read_ptr)
                        {
                            *write_ptr = *read_ptr;
                        }
                        ++write_ptr;
                        bucket_still_has_items = true;
                    }
                    ++read_ptr;
                }
                
                // Update next_cache_ to reflect removed items
                scb.next_cache_ = write_ptr;
                
                // If bucket still has items, keep it in the list
                if (bucket_still_has_items)
                {
                    buckets_to_keep.push_back(bucket_index);
                }
            }
            
            // Update non_empty_buckets_ to only contain buckets that still have unprocessed items
            non_empty_buckets_ = std::move(buckets_to_keep);
        }

    private:
        SieveCacheBucket<cache_size> buckets_[bucket_count];
        LatticeSiever::SIEVE_TYPE* const sieve_array_;
        const BitArray64<sieve_array_size>& sieve_bit_array_;
        std::vector<size_t> non_empty_buckets_;
#ifdef DEBUG_SIEVE_CACHE
        std::ofstream debug_file_;
#endif
    };

    struct PartialFactorisation
    {
        FastVeryLong remaining_quotient_;
        FactorList factor_;

        PartialFactorisation(const VeryLong& v) : remaining_quotient_(v)
        {
            factor_.clear();
        }
        PartialFactorisation() : remaining_quotient_(0L)
        {
            factor_.clear();
        }
    };

    struct PotentiallySmoothPoint
    {
        PotentiallySmoothPoint* next_;
        int32_t c_;
        int32_t d_;
        SIEVE_TYPE* ptr_;
        PartialFactorisation partial1_;
        PartialFactorisation partial2_;
        PotentiallySmoothPoint() : next_(0), c_(0), d_(0), ptr_(0)
        {}
        PotentiallySmoothPoint(long int c, long int d, SIEVE_TYPE* ptr, const VeryLong& v)
            : next_(0), c_(c), d_(d), ptr_(ptr), partial2_(v)
        {}
        int operator<(const PotentiallySmoothPoint& psp) const
        {
            return (ptr_ < psp.ptr_);
        }
        int operator<(const SIEVE_TYPE* ptr) const
        {
            return (ptr_ < ptr);
        }
        static int compare(PotentiallySmoothPoint* psp, SIEVE_TYPE* ptr)
        {
            return (psp->ptr_ < ptr);
        }
        int operator==(const PotentiallySmoothPoint& psp) const
        {
            return (ptr_ == psp.ptr_);
        }
        int operator==(const SIEVE_TYPE* ptr) const
        {
            return (ptr_ == ptr);
        }
        static SIEVE_TYPE* get_sieve_ptr(PotentiallySmoothPoint* psp)
        {
            return psp->ptr_;
        }
        void add_factor1(int32_t p)
        {
            partial1_.remaining_quotient_ /= p;
            partial1_.factor_.push_back(p);
        }
        void add_factor2(int32_t p)
        {
            partial2_.remaining_quotient_ /= p;
            partial2_.factor_.push_back(p);
        }
    };

    void generate_lattice(long int q, long int s);
    void generate_ef_lattice(int32_t p, int32_t r1, std::pair<int32_t, int32_t>& e1, std::pair<int32_t, int32_t>& e2);
    void remove_sieved_factors1();
    void remove_sieved_factors2();
    void eliminate1(long int q);
    void eliminate2();
    void divide_by_small_primes1(PotentiallySmoothPoint* smooth_iter);
    void divide_by_small_primes2(PotentiallySmoothPoint* smooth_iter);

    static bool is_likely_to_be_smooth(const FastVeryLong& remaining_quotient,
                                       long int B, long int L, long int LP, const VeryLong& L_LP);

    static bool is_actually_smooth(FastVeryLong& remaining_quotient, FactorList& factor,
                                   long int B, long int L, long int LP, bool debug);

    void print_relation(long int c, long int d, FactorList& factors1, FactorList& factors2);

    int check_for_remaining_relations();

    void sieve_by_vectors(long int q, long int s);
    void sieve_by_vectors1();
    void sieve_by_vectors1_again();
    void sieve_by_vectors2();
    void sieve1(FactorBase::a_iterator iter, long int r1);
#ifdef RESIEVE1
    void sieve1_again(FactorBase::a_iterator iter, long int r1);
#endif
    void sieve2(FactorBase::a_iterator iter, long int r1);
    long int check_interval1(long int q);
    void check_interval2();
    std::pair<long int, long int> offset_to_c_d(size_t offset);
    size_t c_d_to_offset(const std::pair<long int, long int>& cd);
    bool allocate_c_d_region();
    bool in_range(long int c, long int d);
    void validate_sieve_array(const char* checkpoint_name);

    Polynomial<VeryLong> f1_;
    Polynomial<double> f1d_;
    Polynomial<VeryLong> f2_;
    Polynomial<double> f2d_;
    long int B1_;
    long int L1_;
    long int LP1_;
    VeryLong L_LP_1_;
    long int B2_;
    long int L2_;
    long int LP2_;
    double MIN_A_;
    double MAX_A_;
    long int MIN_B_;
    long int MAX_B_;
    bool debug_;
    bool fixed_sieve_region_;
    VeryLong L_LP_2_;
    long int SIEVE_BOUND_ADJUSTMENT1_;
    long int SIEVE_BOUND_ADJUSTMENT2_;
    long int SMALL_PRIME_BOUND1_;
    long int SMALL_PRIME_BOUND2_;
    long int INITIAL_CUTOFF_;
    double SKEWEDNESS_;
    
    // Precomputed values for check_interval functions
    double L1_pow_LP1_;  // L1^LP1
    double L2_pow_LP2_;  // L2^LP2
    double log_L2_pow_LP2_;  // log10(L2^LP2)

    std::string relation_file_;
    std::fstream* relfile_;
    FactorBase* alg_factor_base_;
    FactorBase* rat_factor_base_;
    
    // Cached small primes for trial division
    static std::vector<long int> small_primes_1_;
    static std::vector<long int> small_primes_2_;
    static bool small_primes_initialized_;
    static void initialize_small_primes(long int bound1, long int bound2);

    static const int LOGQ_BASE = 10;
    static const size_t rat_pf_list_size = 700000L;
#ifdef RESIEVE1
    static const size_t alg_pf_list_size = 50000L;
#else
    static const size_t alg_pf_list_size = 50000000L;
#endif
    static const int max_potentially_smooth = 200000;
#ifdef BUCKET_BITS
    static const size_t bucket_bits = 16;
#else
    static const size_t bucket_size = 57052;
#endif
    static const int c_span_bits = 11;
    //static const int min_c = -8192;
    static const int min_c = -(1 << (c_span_bits - 1));
    //static const int max_c = 8191;
    static const int max_c = (1 << (c_span_bits - 1)) - 1;
    static const int min_d = 0;
    static const int max_d = 4095;
    static const int c_span = max_c - min_c + 1; // 2^14
    static const int fixed_sieve_array_size = c_span * (max_d - min_d + 1); // 2^26
    //static const int sieve_cache_size = 7424;
    static const int sieve_cache_size = 7424;
    
    // Cache blocking parameters
    // Block size chosen to fit in L2 cache (256KB typical)
    // Each sieve entry is 1 byte, so 256K entries = 256KB
    static const size_t CACHE_BLOCK_SIZE = 262144;  // 256KB worth of sieve entries
    static const size_t BLOCKS_PER_SIEVE = (fixed_sieve_array_size + CACHE_BLOCK_SIZE - 1) / CACHE_BLOCK_SIZE;

    SIEVE_TYPE fixed_sieve_array_[fixed_sieve_array_size];
    BitArray64<fixed_sieve_array_size> sieve_bit_array_;

    typedef PointerHashTable<PotentiallySmoothPoint*, SIEVE_TYPE*, 1024> PSPHashTable;
    PotentiallySmoothPoint* potentially_smooth_point_;
    PotentiallySmoothPoint* head_psp_;
    int number_potentially_smooth_;
    std::pair<long int, long int> c1_;
    std::pair<long int, long int> c2_;
    const Parallelogram c_region_;

    PrimeFactorList rat_pf_list_;
    PrimeFactorList alg_pf_list_;
    //BucketedPrimeFactorList<bucket_bits, fixed_sieve_array_size> alg_pf_list_;
#ifdef BUCKET_BITS
    SieveCache<sieve_cache_size, bucket_bits, fixed_sieve_array_size> sieveCache_;
#else
    SieveCache<sieve_cache_size, bucket_size, fixed_sieve_array_size> sieveCache_;
#endif
    Timing timer_;
    static long int total_relations_;
    static double total_sieving_time_;
};
#endif
