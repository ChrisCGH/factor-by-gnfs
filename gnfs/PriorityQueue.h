#ifndef PRIORITYQUEUE_H
#define PRIORITYQUEUE_H

#include <iostream>
#include <utility>
#include <list>
#include <exception>

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
         pq_.push_back(std::make_pair(data, weight));
         if (pq_.size() > S*2)
         {
            sort();
         }
      }

   private:
   public:
      void sort()
      {
         std::sort(pq_.begin(), pq_.end(), [](const std::pair<T, double>& e1, const std::pair<T, double>& e2){ return (e1.second < e2.second); });
         if (pq_.size() > S)
         {
            size_t n = pq_.size() - S;
            pq_.erase(pq_.begin(), pq_.begin() + n - 1);
         }
      }

      void display(std::ostream& os = std::cerr)
      {
         os << "PriorityQueue : allocation = " << S << ", size = " << pq_.size() << std::endl;
         for (auto& wl: pq_)
         {
            os << wl.second << std::endl;
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
