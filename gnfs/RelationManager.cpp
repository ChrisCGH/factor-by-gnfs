#include "RelationManager.h"
#include "MemoryMappedFile.h"
#include <ctype.h>
#include <unordered_set>
#include <string>
#include "Logger.h"

namespace
{
char* buf = 0;
std::string::size_type buflen = 0;
void copy_str(const std::string& str)
{
    if (str.empty()) return;
    if (str.size() > buflen)
    {
        buflen = str.size();
        delete [] buf;
        buf = new char [ buflen + 1 ];
    }
    strcpy(buf, str.c_str());
}
}

RelationManager* RelationSetWeight::relation_sets_ = 0;

RelationSetManager::RelationSetManager() : next_relation_set_(0), relation_set_count_(0), merge_record_file(0), input_relation_set_file(0)
{}

RelationSetManager::RelationSetManager(long int relation_count)
    : next_relation_set_(relation_count), relation_set_count_(relation_count), merge_record_file(0), input_relation_set_file(0)
{}

RelationSetManager::RelationSetManager(long int relation_count, MemoryMappedFile& mmf)
    : next_relation_set_(relation_count), relation_set_count_(relation_count), merge_record_file(0), input_relation_set_file(&mmf)
{
}

RelationSetManager::RelationSetManager(MemoryMappedFile& mmf) : next_relation_set_(0), relation_set_count_(0), merge_record_file(0), input_relation_set_file(0)
{
    LOG_DEBUG("Reading relation sets ...");
    relation_set_relation_map_.read(mmf);

    next_relation_set_ = relation_set_relation_map_.rows();
    relation_set_count_ = next_relation_set_;
    LOG_DEBUG(next_relation_set_ << " relation sets read");
}

RelationSetManager::~RelationSetManager()
{
    relation_set_relation_map_.clear();
    relation_relation_set_map_.clear();
    delete merge_record_file;
}

size_t RelationSetManager::build_relation_relation_set_map(bool clear_relation_set_relation_map)
{
    transpose(relation_set_relation_map_, relation_relation_set_map_, 0, false);
    if (clear_relation_set_relation_map)
        relation_set_relation_map_.clear();
    return relation_relation_set_map_.rows();
}

void RelationSetManager::relation_sets_for_relation(long int relation, std::vector<long int>& relation_sets)
{
    relation_sets.clear();
    if (relation_relation_set_map_.row_size(relation))
    {
        for (auto it = relation_relation_set_map_.begin(relation);
                it != relation_relation_set_map_.end(relation);
                ++it)
        {
            relation_sets.push_back(*it);
        }
    }
}

int RelationSetManager::relation_weight(long int relation_set)
{
    // number of relations in relation set
    int w = relation_set_relation_map_.row_size(relation_set);
    return w;
}

void RelationSetManager::remove(long int relation_set, bool record_only)
{
    if (!record_only)
    {
        relation_set_relation_map_.clear_row(relation_set);
    }
    free_relation_sets_.push_back(relation_set);
    if (record_only)
    {
        record_remove(relation_set);
    }
}

void RelationSetManager::display(long int relation_set)
{
    std::cerr << "Relation set #" << relation_set << " :";
    for (auto it = relation_set_relation_map_.begin(relation_set);
            it != relation_set_relation_map_.end(relation_set);
            ++it)
    {
        std::cerr << " " << *it;
    }
    std::cerr << std::endl;
}

bool RelationSetManager::check(long int rs)
{
    std::unordered_set<long int> s;
    for (auto it = relation_set_relation_map_.begin(rs);
            it != relation_set_relation_map_.end(rs);
            ++it)
    {
        if (s.find(*it) != s.end())
        {
            return false;
        }
        s.insert(*it);
    }
    return true;
}

struct xorer
{
    xorer(SparseMatrix& sm): sm_(sm) {}
    ISparseRow::xor_status operator() (size_t row, size_t col) const
    {
        return sm_.xor(row, col);
    }
    SparseMatrix& sm_;
};

long int RelationSetManager::merge(long int rs1, long int rs2)
{
    std::unordered_set<long int> s;
    for (auto it = relation_set_relation_map_.begin(rs1);
            it != relation_set_relation_map_.end(rs1);
            ++it)
    {
        s.insert(*it);
    }
    for (auto it = relation_set_relation_map_.begin(rs2);
            it != relation_set_relation_map_.end(rs2);
            ++it)
    {
        if (s.find(*it) == s.end())
        {
            s.insert(*it);
        }
        else
        {
            s.erase(*it);
        }
    }

    long int new_rs = next_relation_set_;
    if (!free_relation_sets_.empty())
    {
        new_rs = free_relation_sets_.front();
        free_relation_sets_.pop_front();
    }
    else
    {
        ++next_relation_set_;
    }
    std::for_each(s.begin(), s.end(), [this, new_rs](size_t col) {
        return relation_set_relation_map_.xor(new_rs, col);
    });
    return new_rs;
}

long int RelationSetManager::merge_record_only(long int rs1, long int rs2)
{
    long int new_rs = next_relation_set_;
    if (!free_relation_sets_.empty())
    {
        new_rs = free_relation_sets_.front();
        free_relation_sets_.pop_front();
    }
    else
    {
        ++next_relation_set_;
    }

    record_merge(rs1, rs2, new_rs);
    return new_rs;
}

void RelationSetManager::write_to_file(std::ostream& os)
{
    LOG_DEBUG("Entering ...");
    build_relation_set_relation_map_from_merge_records();
    LOG_DEBUG("relation_set_relation_map_ has " << relation_set_relation_map_.rows() << " rows and " << relation_set_relation_map_.cols() << " columns");
    relation_set_relation_map_.removeEmptyRows();
    LOG_DEBUG("After removing empty rows, relation_set_relation_map_ has " << relation_set_relation_map_.rows() << " rows and " << relation_set_relation_map_.cols() << " columns");
    os << relation_set_relation_map_;
    LOG_DEBUG("Done");
}

void RelationSetManager::record_merge(long int rs1, long int rs2, long int new_rs)
{
    start_merge_recording();
    *merge_record_file << "M:" << rs1 << ":" << rs2 << ":" << new_rs << std::endl;
}

void RelationSetManager::record_remove(long int rs)
{
    start_merge_recording();
    *merge_record_file << "R:" << rs << std::endl;
}

void RelationSetManager::parse_merge_record(const std::string& merge_record)
{
    copy_str(merge_record);
    char* c = buf;
    if (!c) return;
    switch (*c)
    {
    case 'M':
    {
        // M:rs1:rs2:new_rs
        // merge relation sets rs1 and rs2, into new relation set new_rs
        ++c;
        if (!(*c) || *c != ':' || !(*(c+1)))
        {
            LOG_ERROR("Invalid merge record: [" << merge_record << "]");
        }
        ++c;
        long int rs1 = std::atol(c);
        while (*c && *c != ':') ++c;
        if (!(*c) || *c != ':' || !(*(c+1)))
        {
            LOG_ERROR("Invalid merge record: [" << merge_record << "]");
        }
        ++c;
        long int rs2 = std::atol(c);
        while (*c && *c != ':') ++c;
        if (!(*c) || *c != ':' || !(*(c+1)))
        {
            LOG_ERROR("Invalid merge record: [" << merge_record << "]");
        }
        ++c;
        long int new_rs = std::atol(c);
        long int merge_new_rs = RelationSetManager::merge(rs1, rs2);
        if (new_rs != merge_new_rs)
        {
            LOG_ERROR("Problem: new_rs = [" << new_rs << "], merge_new_rs = [" << merge_new_rs << "]");
        }
    }
    break;

    case 'R':
    {
        // R:rs
        // remove relation set rs
        ++c;
        if (!(*c) || *c != ':' || !(*(c+1)))
        {
            LOG_ERROR("Invalid merge record: [" << merge_record << "]");
        }
        ++c;
        long int rs = std::atol(c);
        remove(rs, false);
    }
    break;

    default:
    {
        LOG_ERROR("Invalid merge record: [" << merge_record << "]");
    }
    break;
    }
}


void RelationSetManager::start_merge_recording()
{
    if (!merge_record_file)
    {
        // hard code name for the moment
        merge_record_file = new std::fstream("merge_record.txt", std::ios::out);
    }
}

void RelationSetManager::build_relation_set_relation_map_from_merge_records()
{
    if (input_relation_set_file)
    {
        LOG_DEBUG("Reading input_relation_set_file into relation_set_relation_map_ ...");
        relation_set_relation_map_.read(*input_relation_set_file);
    }
    else
    {
        LOG_DEBUG("Creating relation_set_relation_map_ ...");
        relation_set_relation_map_.set_size(relation_set_count_, relation_set_count_);
        for (int relation = 0; relation < relation_set_count_; ++relation)
        {
            relation_set_relation_map_.xor(relation, relation);
        }
    }
    next_relation_set_ = relation_set_relation_map_.rows();
    relation_set_count_ = next_relation_set_;
    free_relation_sets_.clear();
    if (!merge_record_file) return;
    delete merge_record_file;
    merge_record_file = 0;

    MemoryMappedFile merge_record_mmfile("merge_record.txt");
    LOG_DEBUG("Reading merge_record.txt as a memory mapped file");
    std::string str;
    while (getline(merge_record_mmfile, str))
    {
        parse_merge_record(str);
    }
}

RelationManager::RelationManager(size_t relation_set_count, MemoryMappedFile& relation_sets, MemoryMappedFile& relation_set_primes,
                                 PrimeFrequencyTable& frequency_table, int excess_min)
    : RelationSetManager(relation_set_count, relation_sets),
      relation_set_prime_map_(0, 0),
      prime_relation_set_map_(),
      frequency_table_(frequency_table),
      excess_min_(excess_min)
{
    LOG_DEBUG("relation_set_relation_map_ has " << relation_set_relation_map_.rows() << " rows and " << relation_set_relation_map_.cols() << " columns");
    // read relation_set_prime_map_ from relation_set_primes
    relation_set_prime_map_.read(relation_set_primes);
    LOG_DEBUG("relation_set_prime_map_ has " << relation_set_prime_map_.rows() << " rows and " << relation_set_prime_map_.cols() << " columns");

    // populate frequency table, by iterating through relation_set_prime_map_
    frequency_table_.reset(relation_set_prime_map_.cols());

    for (size_t relation_set = 0; relation_set < relation_set_prime_map_.rows(); ++relation_set)
    {
        for (auto it = relation_set_prime_map_.begin(relation_set);
                it != relation_set_prime_map_.end(relation_set);
                ++it)
        {
            long int prime = *it;
            frequency_table_.increment_prime(prime);
        }
    }
    LOG_DEBUG("PrimeCount = " << frequency_table_.prime_count());
    frequency_table_.check();
    RelationSetWeight::relation_sets_ = this;
}

RelationManager::RelationManager(MemoryMappedFile& is, RelationTable& relation_table, PrimeFrequencyTable& frequency_table, int excess_min)
    : RelationSetManager(is),
      relation_set_prime_map_(relation_set_relation_map_.rows(), 0),
      prime_relation_set_map_(),
      frequency_table_(frequency_table),
      excess_min_(excess_min)
{
    frequency_table_.clear();
    // build relation_set_prime_map_ from relation_set_relation_map_
    // and relation_table
    LOG_DEBUG("relation_set_relation_map_ has " << relation_set_relation_map_.rows() << " rows and " << relation_set_relation_map_.cols() << " columns");
    for (size_t relation_set = 0; relation_set < relation_set_relation_map_.rows(); ++relation_set)
    {
        for (auto it = relation_set_relation_map_.begin(relation_set);
                it != relation_set_relation_map_.end(relation_set);
                ++it)
        {
            Relation& r = relation_table[*it];
            for (int i = 0; i < r.prime_count_; ++i)
            {
                ISparseRow::xor_status rc = relation_set_prime_map_.xor(relation_set, relation_table.get_prime(r.primes_index_, i));
                if (rc == ISparseRow::XOR_REMOVED)
                {
                    frequency_table_.decrement_prime(relation_table.get_prime(r.primes_index_, i));
                }
                else
                {
                    frequency_table_.increment_prime(relation_table.get_prime(r.primes_index_, i));
                }
            }
        }
    }
    relation_table.clear();
    relation_set_prime_map_.set_cols();
    LOG_DEBUG("relation_set_prime_map_ has " << relation_set_prime_map_.rows() << " rows and " << relation_set_prime_map_.cols() << " columns");
    LOG_DEBUG("PrimeCount = " << frequency_table_.prime_count());
    frequency_table_.check();
    RelationSetWeight::relation_sets_ = this;
}

RelationManager::RelationManager(RelationTable& relation_table, PrimeFrequencyTable& frequency_table, int excess_min)
    : RelationSetManager(relation_table.size()),
      relation_set_prime_map_(relation_table.size(), 0),
      prime_relation_set_map_(),
      frequency_table_(frequency_table),
      excess_min_(excess_min)
{
    frequency_table_.check();

    for (size_t relation = 0; relation < relation_table.size(); ++relation)
    {
        long int rs = relation;
        Relation& r = relation_table[relation];
        relation_set_prime_map_.set_row(rs, relation_table.get_prime_start(r.primes_index_), r.prime_count_);
    }
    relation_table.clear();
    LOG_DEBUG("relation_set_prime_map_ has " << relation_set_prime_map_.rows() << " rows and " << relation_set_prime_map_.cols() << " columns");
    RelationSetWeight::relation_sets_ = this;
}

RelationManager::~RelationManager()
{
    relation_set_prime_map_.clear();
    prime_relation_set_map_.clear();
}

void RelationManager::build_prime_relation_set_map(long int merge_level)
{
    transpose(relation_set_prime_map_, prime_relation_set_map_, merge_level, false);

    LOG_DEBUG("prime_relation_set_map_ has " << prime_relation_set_map_.rows() << " rows and " << prime_relation_set_map_.cols() << " columns");
}

void RelationManager::remove(long int relation_set)
{
    RelationSetManager::remove(relation_set);
    SparseRow sr(relation_set_prime_map_.row_size(relation_set));
    relation_set_prime_map_.copy_row(relation_set, sr);
    for (const auto& pr: sr)
    {
        if (prime_relation_set_map_.row_size(pr))
        {
            prime_relation_set_map_.xor(pr, relation_set);
        }
        frequency_table_.decrement_prime(pr);
    }
    relation_set_prime_map_.clear_row(relation_set);
}

long int RelationManager::merge(long int rs1, long int rs2, long int prime)
{
    long int new_rs = RelationSetManager::merge_record_only(rs1, rs2);

    SparseRow sr1(relation_set_prime_map_.row_size(rs1));
    relation_set_prime_map_.copy_row(rs1, sr1);
    SparseRow sr2(relation_set_prime_map_.row_size(rs2));
    relation_set_prime_map_.copy_row(rs2, sr2);

    auto it1 = sr1.begin();
    auto it1_end = sr1.end();
    auto it2 = sr2.begin();
    auto it2_end = sr2.end();

    while (it1 != it1_end && it2 != it2_end)
    {
        if (*it1 < *it2)
        {
            relation_set_prime_map_.xor(new_rs, *it1);
            if (prime_relation_set_map_.row_size(*it1))
            {
                prime_relation_set_map_.xor(*it1, new_rs);
            }
            frequency_table_.increment_prime(*it1);
            ++it1;
        }
        else if (*it2 < *it1)
        {
            relation_set_prime_map_.xor(new_rs, *it2);
            if (prime_relation_set_map_.row_size(*it2))
            {
                prime_relation_set_map_.xor(*it2, new_rs);
            }
            frequency_table_.increment_prime(*it2);
            ++it2;
        }
        else // it1->second == it2->second
        {
            ++it1;
            ++it2;
        }
    }
    while (it1 != it1_end)
    {
        relation_set_prime_map_.xor(new_rs, *it1);
        if (prime_relation_set_map_.row_size(*it1))
        {
            prime_relation_set_map_.xor(*it1, new_rs);
        }
        frequency_table_.increment_prime(*it1);
        ++it1;
    }
    while (it2 != it2_end)
    {
        relation_set_prime_map_.xor(new_rs, *it2);
        if (prime_relation_set_map_.row_size(*it2))
        {
            prime_relation_set_map_.xor(*it2, new_rs);
        }
        frequency_table_.increment_prime(*it2);
        ++it2;
    }

    return new_rs;
}

int RelationManager::combined_weight(int rs1, int rs2)
{
    return prime_weight(rs1) + prime_weight(rs2) - 1;
}

size_t RelationManager::prime_weight(long int relation_set)
{
    // number of primes in relation set
    int w = relation_set_prime_map_.row_size(relation_set);
    return w;
}

size_t RelationManager::sets_including_prime(int prime, std::vector<long int>& relation_sets)
{
    relation_sets.clear();
    if (!prime_relation_set_map_.row_size(prime)) return 0;
    for (auto it = prime_relation_set_map_.begin(prime);
            it != prime_relation_set_map_.end(prime);
            ++it)
    {
        relation_sets.push_back(*it);
    }
    return relation_sets.size();
}

void RelationManager::stats()
{
    long int rs_count = 0;
    long int min_w = 1000000L;
    long int max_w = 0L;
    long long int total_w = 0L;
    for (size_t rs = 0; rs < relation_set_prime_map_.rows(); ++rs)
    {
        if (relation_set_prime_map_.row_size(rs) > 0)
        {
            long int w = prime_weight(rs);
            if (w < min_w) min_w = w;
            if (w > max_w) max_w = w;
            total_w += w;
            ++rs_count;
        }
    }
    double av_w = (double)total_w / (double)rs_count;

    LOG_DEBUG("We have " << rs_count << " relation sets, minimum weight = " << min_w << ", maximum weight = " << max_w << ", average weight = " << av_w);
    LOG_DEBUG("total weight = " << total_w);
    LOG_DEBUG("PrimeCount = " << frequency_table_.prime_count());
    LOG_DEBUG("Memory usage for RelationManager : " << memory_usage());
    minimum_weight_ = min_w;
    maximum_weight_ = max_w;
    relation_set_count_ = rs_count;
    prime_count_ = frequency_table_.prime_count();
}

long int RelationSetManager::memory_usage()
{
    long int s = relation_set_relation_map_.memory_usage("relation_set_relation_map_");
    s += relation_relation_set_map_.memory_usage("relation_relation_set_map_");
    long int free_relation_sets_size = free_relation_sets_.size() * sizeof(long int) + sizeof(free_relation_sets_);
    LOG_DEBUG("RelationSetManager::memory_usage() [free_relation_sets_] : " << free_relation_sets_size << " bytes used");
    s += free_relation_sets_size;
    LOG_DEBUG("RelationSetManager::memory_usage() : " << s << " bytes used");
    return s;
}

void RelationSetManager::clear()
{
    relation_set_relation_map_.clear();
    relation_relation_set_map_.clear();
}

long int RelationManager::memory_usage()
{
    long int s = sizeof(*this);
    s += relation_set_prime_map_.memory_usage("relation_set_prime_map_");
    s += prime_relation_set_map_.memory_usage("prime_relation_set_map_");
    s += RelationSetManager::memory_usage();
    return s;
}

void RelationManager::compress()
{
    LOG_DEBUG("Compressing prime_relation_set_map_ ...");
    prime_relation_set_map_.compress();
    prime_relation_set_map_.memory_usage();
}

void RelationManager::remove_singletons()
{
    bool done = false;
    int relation_set_count;
    while (!done)
    {
        int relation_sets_removed = 0;
        relation_set_count = 0;
        // for each relation set
        for (size_t relation_set = 0; relation_set < relation_set_prime_map_.rows(); ++relation_set)
        {
            if (prime_weight(relation_set) == 0) continue;
            // if relation has a prime with freq == 1 in frequency table
            bool singleton = false;
            for (auto it = relation_set_prime_map_.begin(relation_set);
                    !singleton && it != relation_set_prime_map_.end(relation_set);
                    ++it)
            {
                if (frequency_table_.frequency(*it) == 1)
                {
                    singleton = true;
                }
            }
            if (singleton)
            {
                // remove relation set
                remove(relation_set);
                ++relation_sets_removed;
            }
            else
            {
                ++relation_set_count;
            }
        }
        if (relation_sets_removed == 0) done = true;
        else
        {
            LOG_DEBUG(relation_sets_removed << " singleton relation sets removed");
        }
    }
    LOG_DEBUG(relation_set_count << " relation sets remaining");
    LOG_DEBUG(frequency_table_.prime_count() << " primes remaining");
    long int excess = relation_set_count - frequency_table_.prime_count();
    LOG_DEBUG(excess << " excess relation sets over primes");
}

void RelationManager::remove_heavy_relation_sets()
{
    int excess = relation_set_count_ - prime_count_;
    if (excess < excess_min_ * 2) return;
    // only remove at most 50% of heaviest relations
    LOG_DEBUG("RelationManager::remove_heavy_relation_sets : relation_set_count_ = " << relation_set_count_ << ", prime_count_ = " << prime_count_);
    size_t max_to_remove = (excess) / 2;
    LOG_DEBUG("max_to_remove = " << max_to_remove);
    const int percent = 20;
    long int cutoff_weight = (maximum_weight_ * percent + minimum_weight_ * (100 - percent)) / 100;
    LOG_DEBUG("cutoff_weight = " << cutoff_weight);
    std::vector<long int> to_remove;
    for (size_t rs = 0; to_remove.size() < max_to_remove && rs < relation_set_relation_map_.rows(); ++rs)
    {
        if (relation_set_relation_map_.row_size(rs) == 0) continue;
        if (prime_weight(rs) > static_cast<size_t>(cutoff_weight))
        {
            to_remove.push_back(rs);
        }
    }
    LOG_DEBUG("About to remove " << to_remove.size() << " relation sets");
    for (auto& tr: to_remove)
    {
        remove(tr);
    }
    LOG_DEBUG("About to remove singletons");
    remove_singletons();
}

void RelationManager::write_primes_to_file(std::ostream& os)
{
    LOG_DEBUG("relation_set_prime_map_ has " << relation_set_prime_map_.rows() << " rows and " << relation_set_prime_map_.cols() << " columns");
    relation_set_prime_map_.removeEmptyRows();
    LOG_DEBUG("After removing empty rows, relation_set_prime_map_ has " << relation_set_prime_map_.rows() << " rows and " << relation_set_prime_map_.cols() << " columns");
    os << relation_set_prime_map_;
}

void RelationManager::clear()
{
    relation_set_prime_map_.clear();
    prime_relation_set_map_.clear();
}
