
#ifndef CORE_SHARED_BASE_CONTAINER_NONRESIZABLE_BASE_HPP_
#define CORE_SHARED_BASE_CONTAINER_NONRESIZABLE_BASE_HPP_

#include "../core_ConfigDefs.hpp"

namespace rompp{
namespace core{
    
template<typename derived_type, int ndim>
class ContainerNonResizableBase
  : private core::details::CrtpBase<
  ContainerNonResizableBase<derived_type,ndim>>{
  
private:
  friend derived_type;
  using this_t = ContainerNonResizableBase<derived_type, ndim>;
  friend core::details::CrtpBase<this_t>;
  ContainerNonResizableBase() = default;
  ~ContainerNonResizableBase() = default;
  
};//end class  
} // end namespace core
}//end namespace rompp
#endif
