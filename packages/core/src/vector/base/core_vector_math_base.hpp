
#ifndef CORE_VECTOR_MATH_BASE_HPP_
#define CORE_VECTOR_MATH_BASE_HPP_

#include "../core_vector_traits.hpp"
#include <cmath>

namespace core{
    
template<typename derived_type>
class vectorMathBase
  : private core::details::crtpBase<vectorMathBase<derived_type>>
{
private:
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using der_t = typename details::traits<derived_type>::derived_t;
  using wrap_t = typename details::traits<derived_type>::wrapped_t;

public:
  template <typename op_t>
  void inPlaceOp(op_t op, sc_t a1, sc_t a2, const der_t & other){
    // this = a1*this op a2*other;
    this->underlying().inPlaceOpImpl(op, a1, a2, other);
  }
  void scale(sc_t & factor){
    // this = factor * this
    this->underlying().scaleImpl(factor);
  };  
  void norm1(sc_t & result) const {
    this->underlying().norm1Impl(result);
  };
  void norm2(sc_t & result) const {
    this->underlying().norm2Impl(result);
  };
  void normInf(sc_t & result) const {
    this->underlying().normInfImpl(result);
  };
  void minValue(sc_t & result) const {
    this->underlying().minValueImpl(result);
  };
  void maxValue(sc_t & result) const {
    this->underlying().maxValueImpl(result);
  };

private:
  friend derived_type;
  friend core::details::crtpBase<vectorMathBase<derived_type>>;

  vectorMathBase() = default;
  ~vectorMathBase() = default;

};//end class

  
} // end namespace core
#endif
