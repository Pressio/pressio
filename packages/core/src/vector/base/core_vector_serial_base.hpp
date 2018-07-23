
#ifndef CORE_VECTOR_BASE_VECTOR_SERIAL_BASE_HPP_
#define CORE_VECTOR_BASE_VECTOR_SERIAL_BASE_HPP_

#include "../core_vector_traits.hpp"
#include "../../core_operators_base.hpp"

namespace core{
    
template<typename derived_type>
class VectorSerialBase
  : private core::details::CrtpBase<VectorSerialBase<derived_type>>,
    public SubscriptingOperatorsBase<
	     VectorSerialBase<derived_type>,
             typename details::traits<derived_type>::scalar_t,
	     typename details::traits<derived_type>::ordinal_t>
{

  static_assert(details::traits<derived_type>::isSerial==1,
  "OOPS: distributed concrete vector inheriting from serial base!");
  
private:
  using this_t = VectorSerialBase<derived_type>;
  using ord_t = typename details::traits<derived_type>::ordinal_t;

public:
  ord_t size() const {
    return this->underlying().sizeImpl();
  };
  void resize(ord_t newSize) {
    this->underlying().resizeImpl(newSize);
  };
    
private:
  friend derived_type;
  friend core::details::CrtpBase<this_t>;
  friend SubscriptingOperatorsBase<this_t,
             typename details::traits<derived_type>::scalar_t,
	     typename details::traits<derived_type>::ordinal_t>;

  VectorSerialBase() = default;
  ~VectorSerialBase() = default;
    
};//end class
  
} // end namespace core
#endif
