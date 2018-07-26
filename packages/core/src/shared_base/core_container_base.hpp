
#ifndef CORE_SHARED_BASE_CONTAINER_BASE_HPP_
#define CORE_SHARED_BASE_CONTAINER_BASE_HPP_

#include "core_ConfigDefs.hpp"

namespace core{
    
template<typename derived_type, 
         typename wrapped_t>
class ContainerBase 
  : private core::details::CrtpBase<
  ContainerBase<derived_type, wrapped_t>>
{

public:
  wrapped_t const * data() const {
    return this->underlying().dataImpl();
  }

  wrapped_t * data(){
    return this->underlying().dataImpl();
  }
  
  bool empty() const {
    return this->underlying().emptyImpl();
  };

  void setZero() {
    this->underlying().setZeroImpl();
  }

private:
  friend derived_type;
  friend core::details::CrtpBase<
    ContainerBase<derived_type, wrapped_t>>;

  ContainerBase() = default;
  ~ContainerBase() = default;
      
};//end class
  
} // end namespace core
#endif
