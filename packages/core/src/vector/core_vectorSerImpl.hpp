
#ifndef CORE_VECTORSERIMPL_HPP
#define CORE_VECTORSERIMPL_HPP

#include "core_vectorTraits.hpp"

namespace core
{
    
template<typename derived_type>
class vectorSerImpl
{
public:
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using der_t = typename details::traits<derived_type>::derived_t;
  using wrap_t = typename details::traits<derived_type>::wrapped_t;
  using ord_t = typename details::traits<derived_type>::ordinal_t;

  der_t & underlying(){
    return static_cast<der_t &>(*this);
  };
  der_t const& underlying() const{
    return static_cast<der_t const&>(*this);
  };
  //--------------------------------------------
  
  size_t size() const {
    return this->underlying().sizeImpl();
  };
  void resize(size_t newSize) {
    this->underlying().resizeImpl(newSize);
  };
  
};

    
} // end namespace core

#endif
