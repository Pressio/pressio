
#ifndef CORE_MATRIX_DENSE_DISTRIBUTED_BASE_HPP_
#define CORE_MATRIX_DENSE_DISTRIBUTED_BASE_HPP_

#include "../core_matrix_traits.hpp"


namespace core
{
    
template<typename derived_type>
class matrixDenseDistributedBase
{
public:
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using der_t = typename details::traits<derived_type>::derived_t;
  using wrap_t = typename details::traits<derived_type>::wrapped_t;

private:  
  der_t & underlying(){
    return static_cast<der_t &>(*this);
  };
  der_t const& underlying() const{
    return static_cast<der_t const&>(*this);
  };


  

};

    
} // end namespace core

#endif
