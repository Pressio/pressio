
#ifndef CORE_VECTOR_BASE_VECTOR_MATH_BASE_HPP_
#define CORE_VECTOR_BASE_VECTOR_MATH_BASE_HPP_

#include "../core_vector_traits.hpp"
#include <cmath>
#include <functional>
  
namespace core{
    
template<typename derived_type>
class VectorMathBase
  : private core::details::CrtpBase<
  VectorMathBase<derived_type>>{

  using sc_t = typename details::traits<derived_type>::scalar_t;

public:
  template <typename op_t,
	    typename T = derived_type>
  void inPlaceOp(sc_t a1, sc_t a2,
		 const T & other){
    // this = a1*this op a2*other;
    this->underlying().template inPlaceOpImpl<op_t,T>(a1, a2, other);
  }

  template <typename op_t,
  	    typename T = derived_type>
  void inPlaceOp(sc_t a1, const T & x1,
  		 sc_t a2, const T & x2){
    // this = a1*x1 op a2*x2;
    this->underlying().template inPlaceOpImpl<op_t,T>(a1, x1, a2, x2);
  }
  
  template <typename op0_t,
	    typename T = derived_type,
	    typename op1_t = op0_t,
  	    typename op2_t = op0_t,
  	    typename op3_t = op0_t>
  void inPlaceOp(sc_t a0,
  		 sc_t a1, const T & x1,
  		 sc_t a2, const T & x2,
  		 sc_t a3, const T & x3,
  		 sc_t a4, const T & x4){
    // this = a0 * this op0 (a1*x1) op1 (a2*x2) op2 (a3*x3) op3 (a4*x4)
    this->underlying().template inPlaceOpImpl<op0_t, T, op1_t,
  					      op2_t, op3_t>(a0, a1, x1,
							    a2, x2,
							    a3, x3,
							    a4, x4);
  }
  
  void putScalar(sc_t value) {
    this->underlying().putScalarImpl(value);}

  void scale(sc_t factor){
    // this = factor * this
    this->underlying().scaleImpl(factor);}

  void norm1(sc_t & result) const {
    this->underlying().norm1Impl(result);}

  void norm2(sc_t & result) const {
    this->underlying().norm2Impl(result);}
  
  void normInf(sc_t & result) const {
    this->underlying().normInfImpl(result);}

  void minValue(sc_t & result) const {
    this->underlying().minValueImpl(result);}

  void maxValue(sc_t & result) const {
    this->underlying().maxValueImpl(result);}

private:
  friend derived_type;
  friend core::details::CrtpBase<VectorMathBase<derived_type>>;

  VectorMathBase() = default;
  ~VectorMathBase() = default;

};//end class

  
} // end namespace core
#endif
