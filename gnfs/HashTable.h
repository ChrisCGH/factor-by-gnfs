#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <string>

template <class K, class D>
struct HashTableEntry
{
   K key_;
   D data_;
   HashTableEntry* next_;
   HashTableEntry() : data_(0), next_(0)
   {}
   ~HashTableEntry()
   {
      delete next_;
   }
};

template <class K>
struct HashTableEntry<K, std::string>
{
   K key_;
   std::string data_;
   HashTableEntry* next_;
   HashTableEntry() : next_(0)
   {}
   ~HashTableEntry()
   {
      delete next_;
   }
};

template <class K, class D, class H, unsigned long int S = 200003L>
class HashTable
{
   public:
      HashTable() : hash_table_(0)
      {}
      ~HashTable()
      {
         clear();
         delete [] hash_table_;
      }
      D& operator[] (const K& key)
      {
         HashTableEntry<K, D>* prev = 0;
         HashTableEntry<K, D>* hte = find_entry(key, &prev);
         if (!hte)
         {
            hte = new HashTableEntry<K, D>;
            hte->key_ = key;
            if (prev) prev->next_ = hte;
            else
            {
               unsigned long int h = (H::hash(key) % S);
               hash_table_[h] = hte;
            }
         }
         return hte->data_;
      }

      bool find(const K& key)
      {
         return (find_entry(key) != 0);
      }

      bool end()
      {
         return false;
      }

      void remove(const K& key)
      {
         HashTableEntry<K, D>* prev = 0;
         HashTableEntry<K, D>* hte = find_entry(key, &prev);
         if (hte)
         {
            if (prev)
            {
               prev->next_ = hte->next_;
            }
            else
            {
               unsigned long int h = (H::hash(key) % S);
               hash_table_[h] = hte->next_;
            }
            hte->next_ = 0;
            delete hte;
         }
      }

      void clear()
      {
         if (!hash_table_) return;
         for (size_t i = 0; i < S; ++i)
         {
            HashTableEntry<K, D>* hte = hash_table_[i];
            delete hte;
            hash_table_[i] = 0;
         }
      }

   private:
      void init()
      {
         if (hash_table_) return;
         hash_table_ = new HashTableEntry<K, D>* [ S ];
         for (size_t i = 0; i < S; ++i)
         {
            hash_table_[i] = 0;
         }
      }
      HashTableEntry<K, D>* find_entry(const K& key, HashTableEntry<K, D>** prev = 0)
      {
         init();
         unsigned long int h = (H::hash(key) % S);
         HashTableEntry<K, D>* hte = hash_table_[h];
         while (hte && hte->key_ != key)
         {
            if (prev) *prev = hte;
            hte = hte->next_;
         }
         return hte;
      }

      HashTableEntry<K, D>** hash_table_;
};
#endif
