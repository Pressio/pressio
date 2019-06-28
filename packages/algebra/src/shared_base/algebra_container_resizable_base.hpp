
#ifndef ALGEBRA_SHARED_BASE_CONTAINER_RESIZABLE_BASE_HPP_
#define ALGEBRA_SHARED_BASE_CONTAINER_RESIZABLE_BASE_HPP_

#include "../algebra_ConfigDefs.hpp"

namespace rompp{ namespace algebra{

template<typename derived_type, int ndim>
class ContainerResizableBase
  : private utils::details::CrtpBase<
  ContainerResizableBase<derived_type,ndim>>{

  using this_t = ContainerResizableBase<derived_type, ndim>;

public:
  template <int ndim_ = ndim,
	    typename std::enable_if<ndim_==2,int>::type * = nullptr>
  void resize(std::size_t nrows, std::size_t ncols){
    this->underlying().resizeImpl(nrows, ncols);}

  template <int ndim_ = ndim,
	    typename std::enable_if<ndim_==1,int>::type * = nullptr>
  void resize(std::size_t news){
    this->underlying().resizeImpl(news);}

  void matchLayoutWith(const derived_type & other){
    this->underlying().matchLayoutWithImpl(other);}

private:
  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_type>::type;

  friend utils::details::CrtpBase<this_t>;

  ContainerResizableBase() = default;
  ~ContainerResizableBase() = default;

};//end class

}}//end namespace rompp::algebra
#endif
