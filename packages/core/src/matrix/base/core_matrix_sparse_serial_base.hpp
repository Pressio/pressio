
#ifndef CORE_MATRIX_SPARSE_SERIAL_BASE_HPP_
#define CORE_MATRIX_SPARSE_SERIAL_BASE_HPP_

#include "../core_matrix_traits.hpp"


namespace core
{
    
template<typename derived_type>
class matrixSparseSerialBase
{
public:
  using der_t = typename details::traits<derived_type>::derived_t;
  using wrap_t = typename details::traits<derived_type>::wrapped_t;

private:  
  der_t & underlying(){
    return static_cast<der_t &>(*this); };
  der_t const& underlying() const{
    return static_cast<der_t const&>(*this); };
  

};

    
} // end namespace core

#endif
