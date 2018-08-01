
#ifndef CORE_SHARED_BASE_CONTAINER_DISTRIBUTED_BASE_HPP_
#define CORE_SHARED_BASE_CONTAINER_DISTRIBUTED_BASE_HPP_


#include "core_ConfigDefs.hpp"

namespace core{
    
template<typename derived_type, typename comm_t>
class ContainerDistributedBase
  : private core::details::CrtpBase<
  ContainerDistributedBase<derived_type,comm_t>>
{

public:
  comm_t const & commCRef() const{
    return this->underlying().commCRefImpl();
  }
  
private:
  friend derived_type;
  friend core::details::CrtpBase<
    ContainerDistributedBase<derived_type,comm_t>>;

  ContainerDistributedBase() = default;
  ~ContainerDistributedBase() = default;
  
};//end class  
} // end namespace core
#endif
