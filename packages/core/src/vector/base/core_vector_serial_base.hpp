
#ifndef CORE_VECTOR_SERIAL_BASE_HPP_
#define CORE_VECTOR_SERIAL_BASE_HPP_

#include "../core_vector_traits.hpp"
#include "../../core_operators_base.hpp"

namespace core{
    
template<typename derived_type>
class VectorSerialBase
  : private core::details::CrtpBase<VectorSerialBase<derived_type>>
{
private:
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using der_t = typename details::traits<derived_type>::derived_t;
  using wrap_t = typename details::traits<derived_type>::wrapped_t;
  using ord_t = typename details::traits<derived_type>::ordinal_t;

  static_assert(details::traits<derived_type>::isSerial==1,
		"OOPS: non-serial concrete vector inheriting from serial base!");
  
public:
  ord_t size() const {
    return this->underlying().sizeImpl();
  };
  void resize(ord_t newSize) {
    this->underlying().resizeImpl(newSize);
  };
    
private:
  friend derived_type;
  friend core::details::CrtpBase<VectorSerialBase<derived_type>>;

  VectorSerialBase() = default;
  ~VectorSerialBase() = default;
    
};//end class
  
} // end namespace core
#endif
