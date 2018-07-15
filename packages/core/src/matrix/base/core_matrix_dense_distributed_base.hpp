
#ifndef CORE_MATRIX_DENSE_DISTRIBUTED_BASE_HPP_
#define CORE_MATRIX_DENSE_DISTRIBUTED_BASE_HPP_

#include "../core_matrix_traits.hpp"

namespace core{
    
template<typename derived_type>
class matrixDenseDistributedBase
  : private core::details::crtpBase<matrixDenseDistributedBase<derived_type>>
{
private:
  using der_t = typename details::traits<derived_type>::derived_t;
  using wrap_t = typename details::traits<derived_type>::wrapped_t;
  
private:  
  friend derived_type;
  friend core::details::crtpBase<matrixDenseDistributedBase<derived_type>>;

  matrixDenseDistributedBase() = default;
  ~matrixDenseDistributedBase() = default;
 
};//end class
  
} // end namespace core
#endif
