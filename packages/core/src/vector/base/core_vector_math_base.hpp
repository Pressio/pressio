
#ifndef CORE_VECTOR_BASE_VECTOR_MATH_BASE_HPP_
#define CORE_VECTOR_BASE_VECTOR_MATH_BASE_HPP_

#include "../core_vector_traits.hpp"
#include <cmath>
#include <functional>
  
namespace rompp{
namespace core{
    
template<typename derived_type>
class VectorMathBase
  : private core::details::CrtpBase<
  VectorMathBase<derived_type>>{

  using sc_t = typename details::traits<derived_type>::scalar_t;

public:
  
  template <typename T,
  	    core::meta::enable_if_t<
	      std::is_same<T, sc_t>::value> * = nullptr>
  void putScalar(T value) {
    this->underlying().putScalarImpl(value);}

  template <typename T,
  	    core::meta::enable_if_t<
	      std::is_same<T, sc_t>::value> * = nullptr>
  void scale(T factor){
    // this = factor * this
    this->underlying().scaleImpl(factor);}

  template <typename T,
  	    core::meta::enable_if_t<
	      std::is_same<T, sc_t>::value> * = nullptr>
  void norm1(T & result) const {
    this->underlying().norm1Impl(result);}

  template <typename T,
  	    core::meta::enable_if_t<
	      std::is_same<T, sc_t>::value> * = nullptr>
  void norm2(T & result) const {
    this->underlying().norm2Impl(result);}
  
  template <typename T,
  	    core::meta::enable_if_t<
	      std::is_same<T, sc_t>::value> * = nullptr>
  void normInf(T & result) const {
    this->underlying().normInfImpl(result);}

  template <typename T,
  	    core::meta::enable_if_t<
	      std::is_same<T, sc_t>::value> * = nullptr>
  void minValue(T & result) const {
    this->underlying().minValueImpl(result);}

  template <typename T,
  	    core::meta::enable_if_t<
	      std::is_same<T, sc_t>::value> * = nullptr>
  void maxValue(T & result) const {
    this->underlying().maxValueImpl(result);}

private:
  friend derived_type;
  friend core::details::CrtpBase<VectorMathBase<derived_type>>;

  VectorMathBase() = default;
  ~VectorMathBase() = default;

};//end class

  
} // end namespace core
}//end namespace rompp
#endif
