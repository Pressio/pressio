
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
  // template <typename op_t>
  // void inPlaceOp(op_t op, sc_t a1, sc_t a2,
  // 		 const derived_type & other){
  //   // this = a1*this op a2*other;
  //   this->underlying().inPlaceOpImpl(op, a1, a2, other);
  // }

  // void scale(sc_t & factor){
  //   // this = factor * this
  //   this->underlying().scaleImpl(factor);
  // };  

  // void norm1(sc_t * result) const {
  //   this->underlying().norm1Impl(result);
  // };

  // void norm2(sc_t * result) const {
  //   this->underlying().norm2Impl(result);
  // };

  // void normInf(sc_t * result) const {
  //   this->underlying().normInfImpl(result);
  // };

  // void minValue(sc_t * result) const {
  //   this->underlying().minValueImpl(result);
  // };

  // void maxValue(sc_t * result) const {
  //   this->underlying().maxValueImpl(result);
  // };

private:
  friend derived_type;
  friend core::details::CrtpBase<MultiVectorMathBase<derived_type>>;

  MultiVectorMathBase() = default;
  ~MultiVectorMathBase() = default;

};//end class
  
} // end namespace core
#endif
