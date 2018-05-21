
#ifndef CORE_VECTOR_DIST_BASE_HPP_
#define CORE_VECTOR_DIST_BASE_HPP_

#include "core_vector_traits.hpp"


namespace core
{
    
template<typename derived_type>
class vectorDistBase
{
public:
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using der_t = typename details::traits<derived_type>::derived_t;
  using wrap_t = typename details::traits<derived_type>::wrapped_t;
  using LO_t = typename details::traits<derived_type>::local_ordinal_t;
  using GO_t = typename details::traits<derived_type>::global_ordinal_t;
  using map_t = typename details::traits<derived_type>::map_t;
  using comm_t = typename details::traits<derived_type>::comm_t;

  der_t & underlying(){
    return static_cast<der_t &>(*this);
  };
  der_t const& underlying() const{
    return static_cast<der_t const&>(*this);
  };
  //--------------------------------------------

  size_t globalSize() const {
    return this->underlying().globalSizeImpl();
  };
  size_t localSize() const {
    return this->underlying().localSizeImpl();
  };

  // template<typename U=map_t,
  // 	   typename std::enable_if< !std::is_same<U,void>::value >::type >
  map_t const & getMap() const{
    return this->underlying().getMapImpl();
  }    
  // template<class U=map_t>
  // const typename std::enable_if< !std::is_same<U,void>::value, U>::type & 
  // getMap() const{
  //   return this->underlying().getMapImpl();
  // }
  
};

    
} // end namespace core

#endif
