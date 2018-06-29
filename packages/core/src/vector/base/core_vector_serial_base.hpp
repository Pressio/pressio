
#ifndef CORE_VECTOR_SERIAL_BASE_HPP_
#define CORE_VECTOR_SERIAL_BASE_HPP_

#include "../core_vector_traits.hpp"
#include "../../core_operators_base.hpp"

namespace core{
    
template<typename derived_type>
class vectorSerialBase
{
public:
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using der_t = typename details::traits<derived_type>::derived_t;
  using wrap_t = typename details::traits<derived_type>::wrapped_t;
  using ord_t = typename details::traits<derived_type>::ordinal_t;

  static_assert(details::traits<derived_type>::isSerial==1,
		"OOPS: non-serial concrete vector inheriting from serial base!");
  
public:
  size_t size() const {
    return this->underlying().sizeImpl();
  };
  void resize(size_t newSize) {
    this->underlying().resizeImpl(newSize);
  };
  bool empty() const {
    return this->underlying().emptyImpl();
  };
    
private:
  friend derived_type; 
  vectorSerialBase() = default;
  ~vectorSerialBase() = default;
 
  der_t & underlying(){
    return static_cast<der_t &>(*this);
  };
  der_t const& underlying() const{
    return static_cast<der_t const&>(*this);
  };
   
};//end class
    
} // end namespace core
#endif
