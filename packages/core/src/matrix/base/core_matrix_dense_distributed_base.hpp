
#ifndef CORE_MATRIX_BASE_MATRIX_DENSE_DISTRIBUTED_BASE_HPP_
#define CORE_MATRIX_BASE_MATRIX_DENSE_DISTRIBUTED_BASE_HPP_

#include "../core_matrix_traits.hpp"

namespace core{
    
template<typename derived_type>
class MatrixDenseDistributedBase
  : private core::details::CrtpBase<
  MatrixDenseDistributedBase<derived_type>>{

  static_assert( details::traits<derived_type>::isSharedMem==0,
  "OOPS: non-distributed matrix inheriting from dense distributed base!");
  static_assert( details::traits<derived_type>::isDense==1,
  "OOPS: non-dense matrix inheriting from dense distributed base!");

  using traits_t = details::traits<derived_type>;
  using sc_t = typename traits_t::scalar_t;
  using LO_t = typename traits_t::local_ordinal_t;
  using GO_t = typename traits_t::global_ordinal_t;
  
public:
  void replaceGlobalValue(GO_t globRow, GO_t globCol, sc_t value){
    this->underlying().replaceGlobalValueImpl(globRow, globCol, value);
  }
    
private:  
  friend derived_type;
  friend core::details::CrtpBase<
    MatrixDenseDistributedBase<derived_type>>;

  MatrixDenseDistributedBase() = default;
  ~MatrixDenseDistributedBase() = default;
 
};//end class
  
} // end namespace core
#endif
