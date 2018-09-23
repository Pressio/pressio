
#ifndef CORE_VECTOR_BASE_VECTOR_SHAREDMEM_BASE_HPP_
#define CORE_VECTOR_BASE_VECTOR_SHAREDMEM_BASE_HPP_

#include "../core_vector_traits.hpp"
#include "../../shared_base/core_operators_base.hpp"

namespace rompp{
namespace core{
    
template<typename derived_type>
class VectorSharedMemBase
  : private core::details::CrtpBase<
     VectorSharedMemBase<derived_type>>,
  public Subscripting1DOperatorsBase<
    VectorSharedMemBase<derived_type>, 
    typename details::traits<derived_type>::scalar_t,
    typename details::traits<derived_type>::ordinal_t>
{

  static_assert(details::traits<derived_type>::is_shared_mem==1,
  "OOPS: distributed concrete vector inheriting from sharedMem base!");
  
  using this_t = VectorSharedMemBase<derived_type>;
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using ord_t = typename details::traits<derived_type>::ordinal_t;

public:
  ord_t size() const {
    return this->underlying().sizeImpl();
  };
    
private:
  friend derived_type;
  friend core::details::CrtpBase<this_t>;
  friend Subscripting1DOperatorsBase< this_t, sc_t, ord_t>; 
  VectorSharedMemBase() = default;
  ~VectorSharedMemBase() = default;
    
};//end class
  
} // end namespace core
}//end namespace rompp
#endif
