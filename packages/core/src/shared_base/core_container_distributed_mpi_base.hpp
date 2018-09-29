
#ifdef HAVE_MPI
#ifndef CORE_SHARED_BASE_CONTAINER_DISTRIBUTED_MPI_BASE_HPP_
#define CORE_SHARED_BASE_CONTAINER_DISTRIBUTED_MPI_BASE_HPP_

#include "../core_ConfigDefs.hpp"

namespace rompp{
namespace core{
    
template<typename derived_type, typename comm_t>
class ContainerDistributedMpiBase
  : private core::details::CrtpBase<
  ContainerDistributedMpiBase<derived_type,comm_t> >{
  
public:
  comm_t const & commCRef() const{
    return this->underlying().commCRefImpl();
  }
  
private:
  friend derived_type;
  friend core::details::CrtpBase<
    ContainerDistributedMpiBase<derived_type,comm_t>>;
  ContainerDistributedMpiBase() = default;
  ~ContainerDistributedMpiBase() = default;
  
};//end class  
} // end namespace core
}//end namespace rompp
#endif
#endif