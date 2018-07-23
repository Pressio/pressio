
#ifndef CORE_VECTOR_BASE_VECTOR_DISTRIBUTED_BASE_HPP_
#define CORE_VECTOR_BASE_VECTOR_DISTRIBUTED_BASE_HPP_

#include "../core_vector_traits.hpp"
#include "../../core_operators_base.hpp"

namespace core{
    
template<typename derived_type>
class VectorDistributedBase
  : private core::details::CrtpBase<
              VectorDistributedBase<derived_type>>,
    public SubscriptingOperatorsBase<
	     VectorDistributedBase<derived_type>,
             typename details::traits<derived_type>::scalar_t,
             typename details::traits<derived_type>::local_ordinal_t>
{

  static_assert( details::traits<derived_type>::isDistributed==1,
  "OOPS: serial concrete vector inheriting from distributed base!");

private:
  using this_t = VectorDistributedBase<derived_type>;

  using sc_t = typename details::traits<derived_type>::scalar_t;
  using LO_t = typename details::traits<derived_type>::local_ordinal_t;
  using GO_t = typename details::traits<derived_type>::global_ordinal_t;
  using map_t = typename details::traits<derived_type>::data_map_t;
    
public:
  GO_t globalSize() const {
    return this->underlying().globalSizeImpl();
  };

  LO_t localSize() const {
    return this->underlying().localSizeImpl();
  };

  void replaceGlobalValues(GO_t numentries,
			   const GO_t * indices,
			   const sc_t * values){
    this->underlying().replaceGlobalValuesImpl(numentries, indices, values);
  }
  
private:
  friend derived_type;
  friend core::details::CrtpBase<this_t>;
  friend SubscriptingOperatorsBase<this_t,
             typename details::traits<derived_type>::scalar_t,
	     typename details::traits<derived_type>::local_ordinal_t>;

  VectorDistributedBase() = default;
  ~VectorDistributedBase() = default;

};//end class
} // end namespace core
#endif
