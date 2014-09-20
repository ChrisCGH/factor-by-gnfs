#ifndef POINTERHASHTABLE_H
#define POINTERHASHTABLE_H

template <class POINTER, class KEY, int S=1024>
class PointerHashTable
{
   public:
      typedef KEY (*KEYFUNCTION)(POINTER);
      typedef int (*COMPARE)(POINTER, KEY);
      PointerHashTable(KEYFUNCTION key_function, COMPARE comp) : key_function_(key_function), comp_(comp)
      {}
      ~PointerHashTable()
      {}
      POINTER find(KEY key)
      {
         POINTER_LIST& pl = hash_table_[hash(key)].P_list_;
         POINTER_LIST_ITERATOR found = std::lower_bound(pl.begin(), pl.end(), key, comp_);
         if (found != pl.end() && key_function_(*found) == key) return (*found);
         return 0;
      }
      void clear()
      {
         for (int i = 0; i < S; ++i) hash_table_[i].P_list_.clear();
      }
      void load(POINTER first)
      {
         for (POINTER iter = first;
               iter;
               iter = iter->next_)
         {
            KEY key = key_function_(iter);
            hash_table_[hash(key)].P_list_.push_back(iter);
         }
      }

   private:
      typedef std::vector<POINTER> POINTER_LIST;
      typedef typename POINTER_LIST::iterator POINTER_LIST_ITERATOR;
      struct HashTableEntry
      {
         POINTER_LIST P_list_;
      };
      HashTableEntry hash_table_[S];
      long int hash(KEY key)
      {
         return (reinterpret_cast<long int>(key) & (S - 1));
      }
      KEYFUNCTION key_function_;
      COMPARE comp_;
};
#endif
