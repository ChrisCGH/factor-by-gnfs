#ifndef ROOT_H
#define ROOT_H
#if 1
struct PrimeValuation
{
   PrimeValuation(long int p, int v) : first(p), second(v) {}
   PrimeValuation() : first(0L), second(0) {}
   long int first;
   int second;
};
#else
typedef std::pair<long int, int> PrimeValuation;
#endif

class PrimeList
{
   public:
      typedef PrimeValuation* iterator;
      typedef const PrimeValuation* const_iterator;
      PrimeList() : begin_(0), end_(0)
      {
      }
      iterator begin()
      {
         return begin_;
      }
      iterator end()
      {
         return end_;
      }
      const_iterator begin() const
      {
         return begin_;
      }
      const_iterator end() const
      {
         return end_;
      }
      void reserve(size_t s)
      {
         if (last_ - next_ <= static_cast<int>(s))
         {
            newchunk();
         }
      }
      void push_back(const PrimeValuation& pv)
      {
         if (next_ == last_)
         {
            newchunk();
            iterator out_it = next_;
            for (iterator it = begin_; it != end_; ++it)
            {
               *out_it = *it;
               ++out_it;
            }
            begin_ = next_;
            end_ = out_it;
            next_ = out_it;
         }
         if (!begin_)
         {
            begin_ = next_;
            end_ = next_;
         }
         *end_ = pv;
         ++end_;
         ++next_;
      }
      size_t size()
      {
         return (end() - begin());
      }
   private:
      iterator begin_;
      iterator end_;
      static iterator next_;
      static iterator last_;
      static std::vector<iterator> chunks_;
      static void newchunk()
      {
         const size_t chunk_size = 1024L * 1024L;
         next_ = new PrimeValuation[chunk_size + 1];
         last_ = next_ + chunk_size;
         chunks_.push_back(next_);
      }
      static void clear()
      {
         for (auto& chunk: chunks_)
         {
            delete [] chunk;
         }
         chunks_.clear();
      }
};

struct Relation
{
   long long int a;
   long int b;
   //typedef std::vector<std::pair<long int, int> > primelist;
   typedef PrimeList primelist;
   primelist primes_;
   Relation(long long int aa, long int bb) : a(aa), b(bb)
   {}
};

typedef std::vector<Relation*> RelationList;

#ifdef _NUMBERFIELD_H
void squareRoot(const char* filename, std::vector<AlgebraicNumber*>& delta,
                std::vector<int>& s, AlgebraicNumber& root_gamma_L,
                RelationList& relationNumer, RelationList& relationDenom,
                bool debug, const std::string& dump_file);
#endif

#endif
