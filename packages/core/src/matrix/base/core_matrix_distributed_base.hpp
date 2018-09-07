
#ifndef CORE_MATRIX_BASE_MATRIX_DISTRIBUTED_BASE_HPP_
#define CORE_MATRIX_BASE_MATRIX_DISTRIBUTED_BASE_HPP_

#include "../core_matrix_traits.hpp"

namespace core{
    
template<typename derived_type>
class MatrixDistributedBase
  : private core::details::CrtpBase<
  MatrixDistributedBase<derived_type>>{

  static_assert( details::traits<derived_type>::isSharedMem==0,
   "OOPS: non-distributed matrix inheriting from distributed base!");

  using traits_t = details::traits<derived_type>;
  using LO_t = typename traits_t::local_ordinal_t;
  using GO_t = typename traits_t::global_ordinal_t;

public:
  LO_t localRows() const{
    return this->underlying().localRowsImpl();}

  LO_t localCols() const{
    return this->underlying().localColsImpl();}

  GO_t globalRows() const{
    return this->underlying().globalRowsImpl();}

  GO_t globalCols() const{
    return this->underlying().globalColsImpl();}
  
private:
  friend derived_type;
  friend core::details::CrtpBase<MatrixDistributedBase<derived_type>>;
  MatrixDistributedBase() = default;
  ~MatrixDistributedBase() = default;
  
};//end class  
} // end namespace core
#endif
