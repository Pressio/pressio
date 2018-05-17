
#ifndef CORE_VECTORBASEIMPL_HPP
#define CORE_VECTORBASEIMPL_HPP

#include "core_vectorTraits.hpp"


namespace core
{
    
template<typename derived_type>
class vectorBaseImpl
{
public:
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using der_t = typename details::traits<derived_type>::derived_t;
  using wrap_t = typename details::traits<derived_type>::wrapped_t;

  der_t & underlying(){
    return static_cast<der_t &>(*this);
  };
  der_t const& underlying() const{
    return static_cast<der_t const&>(*this);
  };
  //--------------------------------------------
  
  wrap_t const * view() const {
    return this->underlying().viewImpl();
  };
  wrap_t & getNonConstRefToData() {
    return this->underlying().getNonConstRefToDataImpl();
  };
  
};

    
} // end namespace core

#endif
