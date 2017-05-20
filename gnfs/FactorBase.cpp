#include "FactorBase.h"
#include "LongModular.h"
#include <fstream>
#include <string.h>
#include <stdlib.h>

FactorBase::FactorBase(const Polynomial<VeryLong>& f, long int bound)
      : f_(f), B_(bound), factor_base_size_(0), factor_base_capacity_(0), factor_base_(0)
{
   generate();
}

FactorBase::FactorBase(const Polynomial<VeryLong>& f, long int bound, const char* filename)
      : f_(f), B_(bound), factor_base_size_(0), factor_base_capacity_(0), factor_base_(0)
{
   write_as_we_generate(filename);
   read(filename);
}

FactorBase::FactorBase(const char* filename)
      : factor_base_size_(0), factor_base_capacity_(0), factor_base_(0)
{
   // read factor base from file
   try
   {
      read(filename);
   }
   catch (...)
   {
      std::cerr << "Caught exception" << std::endl;
      throw std::string("Caught exception");
   }
}

FactorBase::~FactorBase()
{
   std::free(factor_base_);
   std::free(root_array_);
}

void FactorBase::write(const char* filename)
{
   // write factor base to file
   std::fstream file(filename, std::ios::out);
   if (!file)
   {
      std::cerr << "Failed to open " << filename << " for writing" << std::endl;
      throw std::string("FactorBase::write : failed to open file");
   }
   // write out coefficients of f_
   file << "Factor base" << std::endl;
   file << f_.deg() << std::endl;
   for (int i = 0; i <= f_.deg(); i++)
   {
      file << f_.coefficient(i) << std::endl;
   }
   // write out bound B_
   file << B_ << std::endl;
   int32_t p = 0;
   // write out factor base
   for (auto iter = begin();
         iter != end();
         ++iter)
   {
      p = iter->p;
      file << p;
      for (auto root_iter = begin(iter);
            root_iter != end(iter);
            ++root_iter)
      {
         file << " " << *root_iter;
      }
      file << std::endl;
   }

   // write out inert primes after a separator
   file << "!" << std::endl;
   for (auto& ip: inert_primes_)
   {
      file << ip << std::endl;
   }
}

void FactorBase::read(const char* filename)
{
   // read factor base from file
   std::fstream file(filename, std::ios::in);
   if (!file)
   {
      std::cerr << "Failed to open " << filename << " for reading" << std::endl;
      throw std::string("FactorBase::read : failed to open file");
   }
   std::vector<VeryLong> c;

   std::string str;
   getline(file, str);
   if ("Factor base" != str)
   {
      file.close();
      std::cerr << "FactorBase::read : bad file" << std::endl;
      throw std::string("FactorBase::read : bad file");
   }
   getline(file, str);
   int degree = atoi(str.c_str());
   //if (degree <= 0 || degree > 6)
   if (degree <= 0)
   {
      file.close();
      std::cerr << "FactorBase::read : bad polynomial degree" << std::endl;
      throw std::string("FactorBase::read : bad polynomial degree");
   }
   c.resize(degree + 1);

   // read polynomial coefficients
   for (int i = 0; i <= degree; i++)
   {
      if (getline(file, str))
      {
         c[i] = str;
      }
      else
      {
         file.close();
         std::cerr << "FactorBase::read : premature end of file" << std::endl;
         throw std::string("FactorBase::read : premature end of file");
      }
   }
   f_ = Polynomial<VeryLong>(c);

   // read bound

   getline(file, str);
   B_ = atol(str.c_str());

   // read factor base
   int32_t p = 0;
   long int count = 0;
   int fb_done = 0;
   static char* buffer = 0;
   static std::string::size_type buflen = 0;
   while (!fb_done && getline(file, str))
   {
      // format:
      // p r1 r2 r3 ...
      if (str.size() > buflen)
      {
         delete [] buffer;
         buflen = str.size();
         buffer = new char [ buflen + 1 ];
      }
      strcpy(buffer, str.c_str());
      if (buffer[0] == '!') fb_done = 1;
      else if (buffer[0])
      {
         add(buffer);
         count++;
      }
   }
   // read inert primes
   int icount = 0;
   while (getline(file, str))
   {
      // format:
      // p r1 r2 r3 ...
      if ("" != str)
      {
         p = atol(str.c_str());
         icount++;
         inert_primes_.push_back(p);
      }
   }

   //std::cerr << "FactorBase::read : " << count << " primes read" << std::endl;
   //std::cerr << "                   " << icount << " inert primes read" << std::endl;
   shrink();
}

void FactorBase::add_extra(int32_t p)
{
   if (exists(p)) return;

   std::cerr << "Adding extra prime <" << p << ">" << std::endl;
   std::vector<LongModular> roots;
   find_roots_mod_p<VeryLong, long int, LongModular>(f_, (long int)p, roots);
   add_extra(p, roots);
}

void FactorBase::queue_extra(int32_t p)
{
   LongModular::set_default_modulus(p);
   if (p < B_) return;
   extra_prime_queue_.push_back(p);
}

void FactorBase::add_extra_queue()
{
   if (extra_prime_queue_.empty()) return;
   std::cerr << "Sorting extra_prime_queue_ (size = " << extra_prime_queue_.size() << ") ... " << std::flush;
   std::sort(extra_prime_queue_.begin(), extra_prime_queue_.end());
   std::cerr << "done" << std::endl;

   std::cerr << "Removing duplicates from extra_prime_queue_ ... " << std::flush;
   auto new_end = std::unique(extra_prime_queue_.begin(), extra_prime_queue_.end());
   size_t unique_primes = new_end - extra_prime_queue_.begin();
   std::cerr << "(size = " << unique_primes << ") done" << std::endl;

   std::cerr << "Adding primes to factor base " << std::flush;
   long int count = 0;
   for (auto extra_iter = extra_prime_queue_.begin();
         extra_iter != new_end;
         ++extra_iter)
   {
      int32_t p = *extra_iter;
      std::vector<LongModular> roots;
      find_roots_mod_p<VeryLong, long int, LongModular>(f_, (long int)p, roots);
      //add_extra(p, roots);
      add(p, roots);
      ++count;
      if (count % 1000 == 0) std::cerr << "." << std::flush;
   }
   std::cerr << " done" << std::endl;
   extra_prime_queue_.clear();
}


void FactorBase::write_as_we_generate(const char* filename)
{
   // write factor base to file
   std::fstream file(filename, std::ios::out);
   if (!file)
   {
      std::cerr << "Failed to open " << filename << " for writing" << std::endl;
      throw std::string("FactorBase::write : failed to open file");
   }
   // write out coefficients of f_
   file << "Factor base" << std::endl;
   file << f_.deg() << std::endl;
   for (int i = 0; i <= f_.deg(); i++)
   {
      file << f_.coefficient(i) << std::endl;
   }
   // write out bound B_
   file << B_ << std::endl;

   // generate factor base for f, i.e. pairs (p,r) with p <= B, f(r) = 0 mod p
   const int32_t first_prime = 2L;
   int32_t p = zpnextb(first_prime);
   std::vector<LongModular> roots;
   std::vector<int32_t> factor_base;

   std::cerr << "FactorBase::write_as_we_generate(): Generating factor base up to " << B_ << " ..." << std::endl;
   long int count = 0;
   while (p <= B_)
   {
//      std::cerr << "p = " << p << std::endl;
      roots.clear();
      factor_base.clear();
      find_roots_mod_p<VeryLong, long int, LongModular>(f_, (long int)p, roots);
      if (roots.size() > 0)
      {
         file << p;
         for (auto& root1: roots)
         {
            int32_t root = root1.get_long();
            size_t i = 0;
            while (i < factor_base.size() && factor_base[i] != root)
            {
               i++;
            }
            if (i >= factor_base.size())
            {
               factor_base.push_back(root);
               file << " " << root;
            }
         }
         // Projective roots
         if (f_.coefficient(f_.deg()) % (long)p == 0)
         {
            factor_base.push_back(p);
         }
         file << std::endl;
      }
      // Check for inert prime
      if (roots.size() == 0 && f_.coefficient(f_.deg()) % (long)p != 0)
      {
         // p could be inert, since there are no roots of f(X) = 0 mod p
         // but we need to check that f(X) is actually irreducible mod p
         std::vector<Polynomial<LongModular> > factors;
         factor_over_F_p<VeryLong, long int, LongModular>(f_, p, factors);
         //factor_over_F_p<VeryLong, VeryLong, VeryLong, LongModular>(f_, p, factors);
         if (factors.size() == 1)
         {
            inert_primes_.push_back(p);
         }
      }
      p = zpnext();
      count++;
      if (count % 1000 == 0) std::cerr << "." << std::flush;
      if (count % 800000 == 0) std::cerr << std::endl;
   }
   // write out inert primes after a separator
   file << "!" << std::endl;
   for (size_t i = 0; i < inert_primes_.size(); i++)
   {
      file << inert_primes_[i] << std::endl;
   }
   std::cerr << std::endl << "... Done" << std::endl;
}

void FactorBase::generate()
{
   // generate factor base for f, i.e. pairs (p,r) with p <= B, f(r) = 0 mod p
   const int32_t first_prime = 2L;
   int32_t p = zpnextb(first_prime);
   std::vector<LongModular> roots;

   std::cerr << "FactorBase::generate(): Generating factor base up to " << B_ << " ..." << std::endl;
   long int count = 0;
   while (p <= B_)
   {
//      std::cerr << "p = " << p << std::endl;
      roots.clear();
      find_roots_mod_p<VeryLong, long int, LongModular>(f_, p, roots);
      add(p, roots);
      // Check for inert prime
      if (roots.size() == 0 && f_.coefficient(f_.deg()) % (long)p != 0)
      {
         // p could be inert, since there are no roots of f(X) = 0 mod p
         // but we need to check that f(X) is actually irreducible mod p
         std::vector<Polynomial<LongModular> > factors;
         factor_over_F_p<VeryLong, long int, LongModular>(f_, p, factors);
         if (factors.size() == 1)
         {
            inert_primes_.push_back(p);
         }
      }
      p = zpnext();
      count++;
      if (count % 1000 == 0) std::cerr << "." << std::flush;
      if (count % 800000 == 0) std::cerr << std::endl;
   }
   std::cerr << std::endl << "... Done" << std::endl;
   shrink();
}

namespace
{
int32_t o_p;
size_t o_count;
int32_t o_root_info[10];
}

FactorBase::a_const_root_iterator FactorBase::begin(int32_t p) const
{
   if (o_p == p) return o_root_info;
   R_p r_p(p);
   auto iter = std::lower_bound(begin(), end(), r_p);
   if (iter != end()) return begin(iter);
   // look in factor_base_overflow_
   auto oiter = factor_base_overflow_.find(p);
   if (oiter != factor_base_overflow_.end())
   {
      o_p = p;
      o_count = oiter->second.size();
      for (size_t i = 0; i < o_count; ++i)
      {
         o_root_info[i] = oiter->second[i];
      }
      return o_root_info;
   }
   return 0;
}

FactorBase::a_const_root_iterator FactorBase::end(int32_t p) const
{
   R_p r_p(p);
   auto iter = std::lower_bound(begin(), end(), r_p);
   if (iter != end()) return end(iter);
   auto root_iter = begin(p);
   if (root_iter) return root_iter + o_count;
   return 0;
}

void FactorBase::dump(const char* filename) const
{
   // write factor base to file
   std::fstream file(filename, std::ios::out);
   if (!file)
   {
      std::cerr << "FactorBase::dump() : Failed to open " << filename << " for writing" << std::endl;
      return;
   }
   // write out coefficients of f_
   file << "Factor base" << std::endl;
   file << f_.deg() << std::endl;
   for (int i = 0; i <= f_.deg(); i++)
   {
      file << f_.coefficient(i) << std::endl;
   }
   // write out bound B_
   file << highest_prime() << std::endl;
   int32_t p = 0;
   // write out factor base
   for (auto iter = begin();
         iter != end();
         ++iter)
   {
      p = iter->p;
      file << p;
      for (auto root_iter = begin(iter);
            root_iter != end(iter);
            ++root_iter)
      {
         file << " " << *root_iter;
      }
      file << std::endl;
   }

   // write out inert primes after a separator
   file << "!" << std::endl;
   for (size_t i = 0; i < inert_primes_.size(); i++)
   {
      file << inert_primes_[i] << std::endl;
   }

   // write out overflow primes after a separator
   file << "!Overflow!" << std::endl;

   for (auto& p1: factor_base_overflow_)
   {
      int32_t p = p1.first;
      file << p;
      const std::vector<int32_t>& roots = p1.second;
      for (auto& root: roots)
      {
         file << " " << root;
      }
      file << std::endl;
   }

}

void FactorBase::clear()
{
   free(root_array_);
   free(factor_base_);
}
