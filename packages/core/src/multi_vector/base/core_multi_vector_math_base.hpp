
#ifndef CORE_MULTIVECTOR_BASE_MULTIVECTOR_MATH_BASE_HPP_
#define CORE_MULTIVECTOR_BASE_MULTIVECTOR_MATH_BASE_HPP_

#include "../core_multi_vector_traits.hpp"

namespace core{
    
template<typename derived_type>
class MultiVectorMathBase
  : private core::details::CrtpBase<MultiVectorMathBase<derived_type>>
{

private:
  using sc_t = typename details::traits<derived_type>::scalar_t;

public:
  void scale(sc_t factor){
    // this = factor * this
    this->underlying().scaleImpl(factor);
  };  

private:
  friend derived_type;
  friend core::details::CrtpBase<MultiVectorMathBase<derived_type>>;

  MultiVectorMathBase() = default;
  ~MultiVectorMathBase() = default;

};//end class
  
} // end namespace core
#endif
