#ifndef POLYNOMIAL_INL
#define POLYNOMIAL_INL
template <> 
inline Polynomial<long int>::Polynomial(const std::vector<std::pair<int, std::string> >& terms)
{
   _coefficients.resize(get_max_power(terms) + 1); 
   for (size_t i = 0; i < terms.size(); ++i)
   {
      _coefficients[terms[i].first] += ::atol(terms[i].second.c_str());
   }
   check_degree();
}

template <> 
inline Polynomial<double>::Polynomial(const std::vector<std::pair<int, std::string> >& terms)
{
   _coefficients.resize(get_max_power(terms) + 1); 
   for (size_t i = 0; i < terms.size(); ++i)
   {
      _coefficients[terms[i].first] += ::atof(terms[i].second.c_str());
   }
   check_degree();
}

template <> 
inline Polynomial<Quotient<long int> >::Polynomial(const std::vector<std::pair<int, std::string> >& terms)
{
   _coefficients.resize(get_max_power(terms) + 1); 
   for (size_t i = 0; i < terms.size(); ++i)
   {
      const char* c = terms[i].second.c_str();
      const char* d = c;
      while (c && *c && *c != '/')
      {
         ++c;
      }
      if (*c == '/')
      {
          _coefficients[terms[i].first] += Quotient<long int>(::atoi(d), ::atoi(c + 1));
      }
      else
      {
          _coefficients[terms[i].first] += Quotient<long int>(::atoi(d));
      }
      //std::cout << "_coefficients[" << terms[i].first << "] = [" << _coefficients[terms[i].first] << "]" << std::endl;
   }
   check_degree();
}
template <> 
inline Polynomial<Quotient<VeryLong> >::Polynomial(const std::vector<std::pair<int, std::string> >& terms)
{
   _coefficients.resize(get_max_power(terms) + 1);  
   for (size_t i = 0; i < terms.size(); ++i)
   {
      const char* c = terms[i].second.c_str();
      const char* d = c;

      while (c && *c && *c != '/')
      {
         ++c;
      }
      if (*c == '/')
      {
          _coefficients[terms[i].first] += Quotient<VeryLong>(VeryLong(std::string(d, c - d)), VeryLong(std::string(c + 1)));
      }
      else
      {
          _coefficients[terms[i].first] += Quotient<VeryLong>(VeryLong(std::string(d)));
      }
   }
   check_degree();
}

#endif
