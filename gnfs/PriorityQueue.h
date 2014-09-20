#ifndef PRIORITYQUEUE_H
#define PRIORITYQUEUE_H

#include <iostream>
#include <utility>
#include <list>
#include <exception>

#if 0
template <class T, int S>
class PriorityQueue
{
   public:
      PriorityQueue() : min_weight_(1e200), next_slot_(0)
      {}
      ~PriorityQueue()
      {}
      void add(const T& data, double weight)
      {
         if (next_slot_ < S)
         {
            pq_[next_slot_].first = data;
            pq_[next_slot_].second = weight;
            if (weight < min_weight_) min_weight_ = weight;
            ++next_slot_;
            if (next_slot_ >= S)
            {
               std::sort(pq_, pq_ + size(), less);
            }
            return;
         }

         // pq_ is full
         // ignore if weight is too small
         if (weight <= min_weight_) return;

         pq_[0].first = data;
         pq_[0].second = weight;
         size_t slot = 1;
         while (slot < S && weight > pq_[slot].second)
         {
            std::swap(pq_[slot - 1], pq_[slot]);
         }
         min_weight_ = pq_[0].second;
      }

   private:
      static bool less(const std::pair<T, double>& e1, const std::pair<T, double>& e2)
      {
         return (e1.second < e2.second);
      }
      static bool more(const std::pair<T, double>& e1, const std::pair<T, double>& e2)
      {
         return (e1.second > e2.second);
      }
   public:
      void sort()
      {
         std::sort(pq_, pq_ + size(), more);
      }

      size_t size()
      {
         return next_slot_;
      }

      void display(std::ostream& os = std::cerr)
      {
         os << "PriorityQueue : allocation = " << S << ", size = " << size() << std::endl;
         for (weighted_list_iterator iter = pq_;
               iter != pq_ + size();
               ++iter)
         {
            os << iter->second << std::endl;
         }
      }

   private:
      typedef typename std::pair<T, double>* weighted_list_iterator;

   public:
      typedef weighted_list_iterator iterator;

      iterator begin()
      {
         return pq_;
      }

      iterator end()
      {
         return pq_ + size();
      }

   private:
      std::pair<T, double> pq_[S];
      double min_weight_;
      size_t next_slot_;
};
#else
template <class T, int S>
class PriorityQueue
{
   public:
      PriorityQueue() : min_weight_(0.0)
      {
         pq_.reserve(2*S);
      }
      ~PriorityQueue()
      {}
      void add(const T& data, double weight)
      {
         if (pq_.size() == 0)
         {
            min_weight_ = weight;
         }
         else if (weight < min_weight_) min_weight_ = weight;
         pq_.push_back(std::make_pair<T, double>(data, weight));
         if (pq_.size() > S*2)
         {
            sort();
         }
      }

   private:
      static bool less(const std::pair<T, double>& e1, const std::pair<T, double>& e2)
      {
         return (e1.second < e2.second);
      }
   public:
      void sort()
      {
         //pq_.sort(less);
         std::sort(pq_.begin(), pq_.end(), less);
         if (pq_.size() > S)
         {
            size_t n = pq_.size() - S;
            pq_.erase(pq_.begin(), pq_.begin() + n - 1);
         }
      }

      void display(std::ostream& os = std::cerr)
      {
         os << "PriorityQueue : allocation = " << S << ", size = " << pq_.size() << std::endl;
         for (weighted_list_iterator iter = pq_.begin();
               iter != pq_.end();
               ++iter)
         {
            os << iter->second << std::endl;
         }
      }

   private:
      //typedef std::list<std::pair<T, double> > weighted_list_type;
      typedef std::vector<std::pair<T, double> > weighted_list_type;
      typedef typename weighted_list_type::iterator weighted_list_iterator;
      typedef typename weighted_list_type::reverse_iterator reverse_weighted_list_iterator;

   public:
      class iterator
      {
         public:
            friend class PriorityQueue;
            iterator& operator++()
            {
               if (curr_ != pq_.rend())
               {
                  ++curr_;
                  ++count_;
                  if (count_ >= S) curr_ = pq_.rend();
               }
               return *this;
            }
            T operator*()
            {
               return curr_->first;
            }
            int operator==(const iterator& it) const
            {
               return (&pq_ == &it.pq_ && curr_ == it.curr_);
            }
            int operator!=(const iterator& it) const
            {
               return !(*this == it);
            }
         private:
            iterator(weighted_list_type& pq, bool first = true) : pq_(pq)
            {
               if (first)
               {
                  curr_ = pq_.rbegin();
                  count_ = 0;
               }
               else curr_ = pq_.rend();
            }
            weighted_list_type& pq_;
            reverse_weighted_list_iterator curr_;
            long int count_;
      };

      iterator begin()
      {
         iterator it(pq_);
         return it;
      }

      iterator end()
      {
         iterator it(pq_, false);
         return it;
      }

      size_t size() const { return pq_.size(); }
   private:
      weighted_list_type pq_;
      double min_weight_;
};
#endif
#endif
