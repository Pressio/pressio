
#ifndef CORE_VECTOR_DISTRIBUTED_BASE_HPP_
#define CORE_VECTOR_DISTRIBUTED_BASE_HPP_

#include "../core_vector_traits.hpp"
#include "core_operators_base.hpp"

namespace core{
    
template<typename derived_type>
class vectorDistributedBase
{
private:
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using der_t = typename details::traits<derived_type>::derived_t;
  using wrap_t = typename details::traits<derived_type>::wrapped_t;
  using LO_t = typename details::traits<derived_type>::local_ordinal_t;
  using GO_t = typename details::traits<derived_type>::global_ordinal_t;
  using map_t = typename details::traits<derived_type>::data_map_t;
  using comm_t =  typename details::traits<derived_type>::communicator_t;

  static_assert( details::traits<derived_type>::isDistributed==1,
		 "OOPS: non-distributed concrete vector inheriting from distributed base!");
    
public:
  size_t globalSize() const {
    return this->underlying().globalSizeImpl();
  };

  size_t localSize() const {
    return this->underlying().localSizeImpl();
  };

  map_t const & getDataMap() const{
    return this->underlying().getDataMapImpl();
  }    

  void putScalar(sc_t value) {
    return this->underlying().putScalarImpl(value);
  }    

private:
  friend derived_type;
  vectorDistributedBase() = default;
  ~vectorDistributedBase() = default;

  der_t & underlying(){
    return static_cast<der_t &>(*this);
  };
  der_t const& underlying() const{
    return static_cast<der_t const&>(*this);
  };

};//end class
} // end namespace core
#endif
