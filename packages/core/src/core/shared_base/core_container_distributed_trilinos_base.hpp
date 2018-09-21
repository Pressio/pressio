
#ifndef CORE_SHARED_BASE_CONTAINER_DISTRIBUTED_TRILINOS_BASE_HPP_
#define CORE_SHARED_BASE_CONTAINER_DISTRIBUTED_TRILINOS_BASE_HPP_

#include "core_ConfigDefs.hpp"

namespace core{
    
template<typename derived_type, typename map_t>
class ContainerDistributedTrilinosBase
  : private core::details::CrtpBase<
  ContainerDistributedTrilinosBase<derived_type, map_t> >{

public:
  map_t const & getDataMap() const{
    return this->underlying().getDataMapImpl();}

  void replaceDataMap(const map_t & mapObj){
    return this->underlying().replaceDataMapImpl(mapObj);}
  
private:
  friend derived_type;
  friend core::details::CrtpBase<
    ContainerDistributedTrilinosBase<derived_type,map_t>>;

  ContainerDistributedTrilinosBase() = default;
  ~ContainerDistributedTrilinosBase() = default;

};//end class
} // end namespace core
#endif
