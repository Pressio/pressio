
#ifndef CORE_SHARED_BASE_CONTAINER_BASE_HPP_
#define CORE_SHARED_BASE_CONTAINER_BASE_HPP_

#include "../core_ConfigDefs.hpp"

namespace rompp{ namespace core{

template<typename derived_type, typename wrapped_t>
class ContainerBase
  : private core::details::CrtpBase<
  ContainerBase<derived_type, wrapped_t>>{

  using this_t 	 = ContainerBase<derived_type, wrapped_t>;
  using scalar_t = typename core::details::traits<derived_type>::scalar_t;

public:
  wrapped_t const * data() const {
    return this->underlying().dataImpl();}

  wrapped_t * data(){
    return this->underlying().dataImpl();}

  wrapped_t dataCp(){
    return this->underlying().dataCpImpl();}

  bool empty() const {
    return this->underlying().emptyImpl();}

  void scale(scalar_t value) {
    this->underlying().scaleImpl(value);
  }

  void setZero() {
    this->underlying().setZeroImpl();}

  template <typename T= derived_type,
	    ::rompp::mpl::enable_if_t<
  	      core::details::traits<T>::is_shared_mem==0,
  	      int> = 0
  	    >
  bool isDistributedGlobally() const{
    return this->underlying().isDistributedGloballyImpl();
  }

  template <typename T= derived_type,
	    ::rompp::mpl::enable_if_t<
  	      core::details::traits<T>::is_shared_mem==1,
	      int> = 0
  	    >
  bool isDistributedGlobally() const{
    return false;
  }

private:
  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_type>::type;

  friend core::details::CrtpBase<this_t>;

  ContainerBase() = default;
  ~ContainerBase() = default;

};//end class

}}//end namespace rompp::core
#endif
