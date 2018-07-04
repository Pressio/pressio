
#ifndef CORE_MATRIX_MATH_BASE_HPP_
#define CORE_MATRIX_MATH_BASE_HPP_

#include "../core_matrix_traits.hpp"

namespace core{
    
template<typename derived_type>
class matrixMathBase
{
private:
  using traits_t = details::traits<derived_type>;
  using der_t = typename traits_t::derived_t;
  using wrap_t = typename traits_t::wrapped_t;
  using ord_t = typename traits_t::ordinal_t;
  using sc_t = typename traits_t::scalar_t;    
public:
  void scale(sc_t factor){
    // this = factor * this
    this->underlying().scaleImpl(factor);
  };  
    
private:  
  friend der_t;
  matrixMathBase() = default;
  ~matrixMathBase() = default; 

private:  
  der_t & underlying(){
    return static_cast<der_t &>(*this);
  };
  der_t const& underlying() const{
    return static_cast<der_t const&>(*this);
  };

};//end class
} // end namespace core
#endif
