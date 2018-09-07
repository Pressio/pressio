
#ifndef CORE_SHARED_BASE_CONTAINER_BASE_HPP_
#define CORE_SHARED_BASE_CONTAINER_BASE_HPP_

#include "core_ConfigDefs.hpp"

namespace core{
    
template<typename derived_type, 
         typename wrapped_t>
class ContainerBase 
  : private core::details::CrtpBase<
  ContainerBase<derived_type, wrapped_t>>{

public:
  wrapped_t const * data() const {
    return this->underlying().dataImpl();}

  wrapped_t * data(){
    return this->underlying().dataImpl();}

  wrapped_t dataCp(){
    return this->underlying().dataCpImpl();}
  
  bool empty() const {
    return this->underlying().emptyImpl();}

  void setZero() {
    this->underlying().setZeroImpl();}

  void matchLayoutWith(const derived_type & other){
    this->underlying().matchLayoutWithImpl(other);}

  template <typename T= derived_type,
	    typename std::enable_if<
  	      core::details::traits<T>::isSharedMem==0
  	      ,int>::type = 0
  	    >
  bool isDistributedGlobally() const{
    return this->underlying().isDistributedGloballyImpl();
  }

  template <typename T= derived_type,
	    typename std::enable_if<
  	      core::details::traits<T>::isSharedMem==1
  	      ,int>::type = 0
  	    >
  bool isDistributedGlobally() const{
    return false;
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
