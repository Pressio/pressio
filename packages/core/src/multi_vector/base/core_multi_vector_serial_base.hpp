
#ifndef CORE_MULTIVECTOR_BASE_MULTIVECTOR_SERIAL_BASE_HPP_
#define CORE_MULTIVECTOR_BASE_MULTIVECTOR_SERIAL_BASE_HPP_

#include "../core_multi_vector_traits.hpp"

namespace core{
    
template<typename derived_type>
class MultiVectorSerialBase
  : private core::details::CrtpBase<MultiVectorSerialBase<derived_type>>
{
  static_assert(details::traits<derived_type>::isSerial==1,
  "OOPS: distributed concrete vector inheriting from serial base!");
 
private:
	using this_t = MultiVectorSerialBase<derived_type>;
public:
    
private:
  friend derived_type;
  friend core::details::CrtpBase<this_t>;

  MultiVectorSerialBase() = default;
  ~MultiVectorSerialBase() = default;
    
};//end class
  
} // end namespace core
#endif
