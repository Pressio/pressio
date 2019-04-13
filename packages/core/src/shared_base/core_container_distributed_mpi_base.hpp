
#ifdef HAVE_MPI
#ifndef CORE_SHARED_BASE_CONTAINER_DISTRIBUTED_MPI_BASE_HPP_
#define CORE_SHARED_BASE_CONTAINER_DISTRIBUTED_MPI_BASE_HPP_

#include "../core_ConfigDefs.hpp"
#include "../meta/core_meta_basic.hpp"

namespace rompp{ namespace core{

template<typename derived_type, typename comm_t>
class ContainerDistributedMpiBase
  : private core::details::CrtpBase<
  ContainerDistributedMpiBase<derived_type,comm_t> >{

  using this_t = ContainerDistributedMpiBase<derived_type,comm_t>;

public:

  template <typename T = comm_t,
	    typename std::enable_if<
	      !meta::is_teuchos_rcp<T>::value and 
        !::rompp::mpl::is_std_shared_ptr<T>::value 
	      >::type * = nullptr
	    >
  T const & commCRef() const{
    return this->underlying().commCRefImpl();
  }

  template <typename T= comm_t,
  	    typename std::enable_if<
        meta::is_teuchos_rcp<T>::value or 
        ::rompp::mpl::is_std_shared_ptr<T>::value 
  	      >::type * = nullptr>
  T comm() const{
    return this->underlying().commImpl();
  }

private:
  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_type>::type;

  friend core::details::CrtpBase<this_t>;

  ContainerDistributedMpiBase() = default;
  ~ContainerDistributedMpiBase() = default;

};//end class

}}//end namespace rompp::core
#endif
#endif
