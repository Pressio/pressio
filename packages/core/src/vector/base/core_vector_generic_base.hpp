
#ifndef CORE_VECTOR_GENERIC_BASE_HPP_
#define CORE_VECTOR_GENERIC_BASE_HPP_

#include "../core_vector_traits.hpp"
#include "../../core_operators_base.hpp"

namespace core{
    
template<typename derived_type>
class vectorGenericBase
  : public subscriptingOperatorsBase<
  vectorGenericBase<derived_type>,
    typename details::traits<derived_type>::scalar_t,
    //select ordinal type based on whether the vector is serial or distributed
    typename std::conditional<details::traits<derived_type>::isSerial==1,
			      typename
			      details::traits<derived_type>::ordinal_t,
			      typename
			      details::traits<derived_type>::local_ordinal_t
			      >::type>
{
private:
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using der_t = typename details::traits<derived_type>::derived_t;
  using wrap_t = typename details::traits<derived_type>::wrapped_t;

public:
  wrap_t const * data() const {
    return this->underlying().dataImpl();
  }

  wrap_t * data(){
    return this->underlying().dataImpl();
  }

  void putScalar(sc_t value) {
    this->underlying().putScalarImpl(value);
  }

  void setZero() {
    this->underlying().setZeroImpl();
  }
  
  bool empty() const {
    return this->underlying().emptyImpl();
  };

  // inherits also subscripting operator [] from base (see above)

private:
  friend derived_type;
  vectorGenericBase() = default;
  ~vectorGenericBase() = default;
  
  der_t & underlying(){
    return static_cast<der_t &>(*this);
  }
  der_t const& underlying() const{
    return static_cast<der_t const&>(*this);
  }
    
};//end class    
} // end namespace core
#endif
