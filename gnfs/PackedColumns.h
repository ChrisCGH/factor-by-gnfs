#ifndef PACKEDCOLUMNS_H
#define PACKEDCOLUMNS_H
#include "BitOperations.h"
#include <list>
#include <stdint.h>

class PackedColumns;
class PackedColumnsConstIterator
{
   public:
      friend class PackedColumns;
      PackedColumnsConstIterator& operator++();
      long int operator*() const;
      bool operator==(const PackedColumnsConstIterator& it) const;
      bool operator!=(const PackedColumnsConstIterator& it) const
      {
         return !(*this == it);
      }
      PackedColumnsConstIterator() : first_column_(0), packed_columns_(0), current_bit_(0)
      {}
   private:
      PackedColumnsConstIterator(const PackedColumns& pc, bool end = false);
      long int first_column_;
      uint32_t packed_columns_;
      int current_bit_;
};

class PackedColumns
{
   public:
      friend class PackedColumnsConstIterator;
      PackedColumns() : first_column_(-BitOperations::BITS_IN_WORD), packed_columns_(0)
      {}
      bool xor(size_t column)
      {
         bool removed = false;
         if (first_column_ == -BitOperations::BITS_IN_WORD) first_column_ = static_cast<long int>(column);
         if (in_range(column))
         {
            size_t bit = column - first_column_;
            if (BitOperations::bitSet(bit, packed_columns_))
            {
               BitOperations::clearBit(bit, packed_columns_);
               removed = true;
            }
            else
            {
               BitOperations::setBit(bit, packed_columns_);
            }
         }
         return removed;
      }
      bool in_range(size_t column) const
      {
         if (static_cast<long int>(column) >= first_column_ && static_cast<long int>(column) < first_column_ + BitOperations::BITS_IN_WORD) return true;
         else return false;
      }
      size_t first_column() const
      {
         return first_column_;
      }
      bool operator<(const PackedColumns& pc) const
      {
         return (first_column_ < pc.first_column_);
      }
      long int highest_column() const
      {
         if (first_column_ < 0) return -1;
         int bit = BitOperations::highestBitSet(packed_columns_);
         return (first_column_ + bit);
      }
      uint32_t packed_columns() const
      {
         return packed_columns_;
      }
      size_t size() const
      {
         return BitOperations::bitCount(packed_columns_);
      }
      typedef PackedColumnsConstIterator const_iterator;

      const_iterator begin() const
      {
         return const_iterator(*this);
      }

      const_iterator end() const
      {
         return const_iterator(*this, true);
      }
   private:
      long int first_column_;
      uint32_t packed_columns_;
};

inline PackedColumnsConstIterator::PackedColumnsConstIterator(const PackedColumns& pc, bool end)
      : first_column_(pc.first_column_), packed_columns_(pc.packed_columns_), current_bit_(0)
{
   if (end)
   {
      current_bit_ = BitOperations::BITS_IN_WORD;
   }
}

inline PackedColumnsConstIterator& PackedColumnsConstIterator::operator++()
{
   if (current_bit_ < BitOperations::BITS_IN_WORD)
   {
      ++current_bit_;
      while (current_bit_ < BitOperations::BITS_IN_WORD &&
             !BitOperations::bitSet(current_bit_, packed_columns_))
      {
         ++current_bit_;
      }
   }
   return *this;
}

inline long int PackedColumnsConstIterator::operator*() const
{
   //if (current_bit_ == BitOperations::BITS_IN_WORD) return -1;
   return (current_bit_ + first_column_);
}

inline bool PackedColumnsConstIterator::operator==(const PackedColumnsConstIterator& it) const
{
   return (current_bit_ == it.current_bit_ && first_column_ == it.first_column_ && packed_columns_ == it.packed_columns_);
}

typedef std::list<PackedColumns> packed_list_type;
typedef packed_list_type::const_iterator packed_list_const_iterator;
typedef packed_list_type::iterator packed_list_iterator;
typedef packed_list_type::reverse_iterator packed_list_reverse_iterator;
#endif
