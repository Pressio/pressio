
#ifndef CORE_SHARED_BASE_CONTAINER_SUBSCRIPTABLE_BASE_HPP_
#define CORE_SHARED_BASE_CONTAINER_SUBSCRIPTABLE_BASE_HPP_

#include "../core_ConfigDefs.hpp"

namespace rompp{ namespace core{


template<typename derived_type,
	 typename scalar_t,
	 typename ord_t>
class ContainerSubscriptable1DBase
  : private core::details::CrtpBase<
  ContainerSubscriptable1DBase<derived_type, scalar_t, ord_t>>{

  using this_t = ContainerSubscriptable1DBase<derived_type, scalar_t, ord_t>;

public:
  scalar_t & operator[] (ord_t i){
    return this->underlying()[i];
  }

  scalar_t const & operator[] (ord_t i) const{
    return this->underlying()[i];
  }

  scalar_t & operator() (ord_t i){
    return this->underlying()(i);
  }

  scalar_t const & operator() (ord_t i) const{
    return this->underlying()(i);
  }

private:
  friend derived_type;
  friend core::details::CrtpBase<this_t>;

  ContainerSubscriptable1DBase() = default;
  ~ContainerSubscriptable1DBase() = default;

};//end class
//-------------------------------------------------------



template<typename derived_type,
	 typename scalar_t,
	 typename ord1_t,
	 typename ord2_t = ord1_t>
class ContainerSubscriptable2DBase
  : private core::details::CrtpBase<
  ContainerSubscriptable2DBase<derived_type,
			       scalar_t,
			       ord1_t,
			       ord2_t>>{

  using this_t = ContainerSubscriptable2DBase<derived_type,
					      scalar_t,
					      ord1_t,
					      ord2_t>;

public:
  scalar_t & operator() (ord1_t i, ord2_t j){
    return this->underlying()(i,j);
  }

  scalar_t const & operator() (ord1_t i, ord2_t j) const{
    return this->underlying()(i,j);
  }

private:
  friend derived_type;
  friend core::details::CrtpBase<this_t>;

  ContainerSubscriptable2DBase() = default;
  ~ContainerSubscriptable2DBase() = default;

};//end class


}}//end namespace rompp::core
#endif
