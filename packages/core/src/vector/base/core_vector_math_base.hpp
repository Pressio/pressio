
#ifndef CORE_VECTOR_MATH_BASE_HPP_
#define CORE_VECTOR_MATH_BASE_HPP_

#include "../core_vector_traits.hpp"

namespace core{
    
template<typename derived_type>
class vectorMathBase
{
public:
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using der_t = typename details::traits<derived_type>::derived_t;
  using wrap_t = typename details::traits<derived_type>::wrapped_t;

public:

  // template <typename op_t>
  // void inPlaceOp(op_t op, sc_t a1, sc_t a2, const der_t & vin){
  //   // this = a1*this op a2*vin;
  //   this->underlying().inPlaceOpImpl(op,a1,a2,vin);
  // }

  // sc_t norm2() const {
  //   return this->underlying().norm2Impl();
  // };

private:
  friend derived_type;
  vectorMathBase() = default;
  ~vectorMathBase() = default;

  der_t & underlying(){
    return static_cast<der_t &>(*this);
  };
  der_t const& underlying() const{
    return static_cast<der_t const&>(*this);
  };

};//end class

} // end namespace core
#endif
