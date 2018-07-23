
#ifndef CORE_MATRIX_BASE_MATRIX_DENSE_DISTRIBUTED_BASE_HPP_
#define CORE_MATRIX_BASE_MATRIX_DENSE_DISTRIBUTED_BASE_HPP_

#include "../core_matrix_traits.hpp"

namespace core{
    
template<typename derived_type>
class MatrixDenseDistributedBase
  : private core::details::CrtpBase<
  MatrixDenseDistributedBase<derived_type>>
{

  static_assert( details::traits<derived_type>::isDistributed==1,
		 "OOPS: non-distributed matrix inheriting \
from dense distributed base!");
  
private:  
  friend derived_type;
  friend core::details::CrtpBase<
    MatrixDenseDistributedBase<derived_type>>;

  MatrixDenseDistributedBase() = default;
  ~MatrixDenseDistributedBase() = default;
 
};//end class
  
} // end namespace core
#endif
