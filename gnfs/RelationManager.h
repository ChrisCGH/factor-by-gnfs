#ifndef RELATIONSETMANAGER_H
#define RELATIONSETMANAGER_H
#include "MemoryMappedFile.h"
#include <iostream>
#include <vector>
#include <queue>
#include <map>
#include <algorithm>
#include <stdlib.h>
#include "SparseMatrix.h"

// manages prime frequency table
const size_t MAX_PRIMES = 40000000;
class PrimeFrequencyTable
{
   public:
      PrimeFrequencyTable(size_t max_primes = MAX_PRIMES) : prime_index_(0), prime_count_(0), max_primes_(max_primes), max_frequency_(0)
      {
         frequency_table_ = new freq_type [max_primes];
         for (size_t i = 0; i < max_primes; ++i) frequency_table_[i] = 0;
      }

      ~PrimeFrequencyTable()
      {
         delete [] frequency_table_;
      }

      void clear()
      {
         prime_count_ = 0;
         max_frequency_ = 0;
         delete [] frequency_table_;
         frequency_table_ = new freq_type [prime_index_];
         max_primes_ = prime_index_;
         for (size_t i = 0; i < prime_index_; ++i) frequency_table_[i] = 0;
      }

      void reset(long int prime_count)
      {
         prime_index_ = prime_count;
         clear();
      }

      void add_new_prime(long int prime)
      {
         if (prime_index_ >= max_primes_)
         {
            std::cerr << "Problem: FrequencyTable full, prime_index = " << prime_index_ << ", max_primes_ = " << max_primes_ << std::endl;
            throw "FrequencyTable full";
         }
         frequency_table_[prime_index_] = 0;
         ++prime_index_;
      }

      void increment_prime(long int prime)
      {
         long freq = ++frequency_table_[prime];
         if (frequency_table_[prime] == 1) ++prime_count_;
         if (freq > max_frequency_)
         {
            max_frequency_ = freq;
            if (max_frequency_ == 0xFFFF)
            {
               std::cerr << "Look out! : prime = " << prime << ", max_frequency_ = " << max_frequency_ << std::endl;
               throw "frequency overflow";
            }
         }
      }

      void decrement_prime(long int prime)
      {
         if (frequency_table_[prime] == 0)
         {
            std::cerr << "Look out! : prime = " << prime << ", negative frequency" << std::endl;
            throw "negative frequency";
         }
         --frequency_table_[prime];
         if (frequency_table_[prime] == 0)
            --prime_count_;
      }

      size_t next_prime() const
      {
         return prime_index_;
      }

      size_t prime_count() const
      {
         return prime_count_;
      }

      size_t frequency(long int prime) const
      {
         return frequency_table_[prime];
      }

      size_t capacity() const
      {
         return prime_index_;
      }

      bool check() const
      {
         size_t pc = 0;
         for (size_t i = 0; i < prime_index_; ++i)
         {
            if (frequency_table_[i] > 0) ++pc;
         }
         if (pc != prime_count_)
         {
            std::cerr << "Problem: PrimeFrequencyTable::check : prime_count_ = " << prime_count_ << ", pc = " << pc << std::endl;
            return false;
         }
         return true;
      }

   private:
      typedef unsigned short int freq_type;
      freq_type* frequency_table_;
      size_t prime_index_;
      size_t prime_count_;
      size_t max_primes_;
      long int max_frequency_;
};

struct Relation
{
   int32_t id_;
   int32_t primes_index_;
   unsigned char prime_count_;
   Relation(int32_t id) : id_(id), primes_index_(-1), prime_count_(0)
   {}
   Relation() : id_(-1), primes_index_(-1), prime_count_(0)
   {}
   void clear()
   {
      id_ = -1;
      primes_index_ = -1;
      prime_count_ = 0;
   }
   bool is_clear()
   {
      return id_ == -1;
   }
};

class RelationTable
{
   public:
      RelationTable(size_t max_relations, size_t max_primes) : max_relations_(max_relations), max_primes_(max_primes)
      {
         //relation_table_ = new Relation [ max_relations ];
         relation_table_ = static_cast<Relation*>(::malloc(max_relations * sizeof(Relation)));
         //relation_primes_ = new int [ max_primes ];
         relation_primes_ = static_cast<int*>(::malloc(max_primes * sizeof(int)));
         relation_primes_ptr_ = relation_primes_;
      }

      Relation* begin()
      {
         return relation_table_;
      }

      Relation* end()
      {
         return relation_table_ + relation_count_;
      }

      Relation& operator[](size_t i)
      {
         if (i >= max_relations_)
         {
            std::cerr << "Problem: relationTable.relation_table_ is full" << std::endl;
            throw "relationTable.relation_table_ is full";
         }
         return relation_table_[i];
      }

      int32_t next_prime_index() const
      {
         return static_cast<int32_t>(relation_primes_ptr_ - relation_primes_);
      }

      int get_prime(int32_t start, unsigned char index) const
      {
         return (relation_primes_ + start)[index];
      }

      const int* get_prime_start(int32_t start) const
      {
         return relation_primes_ + start;
      }

      void set_prime(int32_t start, unsigned char index, int prime)
      {
         (relation_primes_ + start)[index] = prime;
      }

      void increment_next_prime(size_t count)
      {
         if (relation_primes_ptr_ + count >= relation_primes_ + max_primes_)
         {
            std::cerr << "Problem: relationTable.relation_primes_ is full : count = " << count << ", max_primes_ = " << max_primes_ << std::endl;
            throw "relationTable.relation_primes_ is full";
         }
         relation_primes_ptr_ += count;
      }

      void clear()
      {
         //delete [] relation_table_;
         free(relation_table_);
         //delete [] relation_primes_;
         free(relation_primes_);
      }

      size_t size() const
      {
         return relation_count_;
      }

      void set_size(size_t s)
      {
         relation_count_ = s;
      }

      void shrink_to_fit()
      {
         size_t prime_count = relation_primes_ptr_ - relation_primes_;
         if (prime_count < max_primes_)
         {
            std::cerr << "Shrinking relation_primes_, from " << max_primes_ << " to " << prime_count << std::endl;
            relation_primes_ = static_cast<int*>(::realloc(relation_primes_, prime_count * sizeof(int)));
         }
         if (relation_count_ < max_relations_)
         {
            std::cerr << "Shrinking relation_table_, from " << max_relations_ << " to " << relation_count_ << std::endl;
            relation_table_ = static_cast<Relation*>(::realloc(relation_table_, relation_count_ * sizeof(Relation)));
         }
      }

   private:
      int* relation_primes_;
      int* relation_primes_ptr_;
      Relation* relation_table_;
      size_t relation_count_;
      size_t max_relations_;
      size_t max_primes_;
};

// manages relation sets
class RelationSetManager
{
   public:
      RelationSetManager();
      RelationSetManager(long int relation_count);
      RelationSetManager(MemoryMappedFile& mmf);
      RelationSetManager(long int relation_count, MemoryMappedFile& mmf);
      ~RelationSetManager();
      size_t build_relation_relation_set_map(bool clear_relation_set_relation_map = false);
      int relation_weight(long int relation_set);
      void remove(long int relation_set, bool record_only = true);
      long int merge(long int rs1, long int rs2);
      long int merge_record_only(long int rs1, long int rs2);
      void write_to_file(std::ostream& os);
      long int number_of_relation_sets() const
      {
         return next_relation_set_;
      }
      void relation_sets_for_relation(long int relation, std::vector<long int>& relation_sets);
      void display(long int relation_set);
      bool check(long int relation_set);
      long int memory_usage();
      void clear();
   private:
      long int next_relation_set_;
      long int relation_set_count_;
      std::deque<long int> free_relation_sets_;
      void record_merge(long int rs1, long int rs2, long int new_rs);
      void record_remove(long int rs);
      void start_merge_recording();
      void build_relation_set_relation_map_from_merge_records();
      void parse_merge_record(const std::string& merge_record);
      std::fstream* merge_record_file;
      MemoryMappedFile* input_relation_set_file;
   protected:
      SparseMatrix relation_set_relation_map_; // maps relation set -> relations
      SparseMatrix relation_relation_set_map_; // maps relations -> relation sets
};

class RelationManager : public RelationSetManager
{
   public:
      RelationManager(RelationTable& relation_table, PrimeFrequencyTable& frequency_table, int excess_min);
      RelationManager(MemoryMappedFile& relation_sets,
                      RelationTable& relation_table, PrimeFrequencyTable& frequency_table, int excess_min);
      RelationManager(size_t relation_set_count, MemoryMappedFile& relation_sets, MemoryMappedFile& relation_set_primes, PrimeFrequencyTable& frequency_table, int excess_min);
      ~RelationManager();
      void remove_singletons();
      void remove_heavy_relation_sets();
      void remove(long int relation_set);
      long int merge(long int rs1, long int rs2, long int prime);
      int combined_weight(int rs1, int rs2);
      size_t prime_weight(long int relation_set);
      size_t sets_including_prime(int prime, std::vector<long int>& relation_sets);
      void stats();
      long int memory_usage();
      void compress();
      void build_prime_relation_set_map(long int merge_level);
      void write_primes_to_file(std::ostream& os);
      void clear();
   private:
      SparseMatrix relation_set_prime_map_; // maps relations sets -> primes
      SparseMatrix prime_relation_set_map_; // maps primes -> relations sets

      PrimeFrequencyTable& frequency_table_;
      int excess_min_;
      long int minimum_weight_;
      long int maximum_weight_;
      long int relation_set_count_;
      long int prime_count_;
};

struct RelationSetWeight : public std::binary_function<long int, long int, double>
{
   double operator() (long int rs1, long int rs2)
   {
      return relation_sets_->combined_weight(rs1, rs2);
   }
   static RelationManager* relation_sets_;
};
#endif
