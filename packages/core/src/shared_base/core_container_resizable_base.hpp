
#ifndef CORE_SHARED_BASE_CONTAINER_RESIZABLE_BASE_HPP_
#define CORE_SHARED_BASE_CONTAINER_RESIZABLE_BASE_HPP_

#include "core_ConfigDefs.hpp"

namespace core{
    
template<typename derived_type, int ndim>
class ContainerResizableBase
  : private core::details::CrtpBase<
  ContainerResizableBase<derived_type,ndim>>{
  
public:
  template <int ndim_ = ndim,
	    typename std::enable_if<ndim_==2,int>::type * = nullptr>
  void resize(std::size_t nrows, std::size_t ncols){
    this->underlying().resizeImpl(nrows, ncols);}

  template <int ndim_ = ndim,
	    typename std::enable_if<ndim_==1,int>::type * = nullptr>
  void resize(std::size_t news){
    this->underlying().resizeImpl(news);}
  
private:
  friend derived_type;
  using this_t = ContainerResizableBase<derived_type, ndim>;
  friend core::details::CrtpBase<this_t>;
  ContainerResizableBase() = default;
  ~ContainerResizableBase() = default;
  
};//end class  
} // end namespace core
#endif
