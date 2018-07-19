
#ifndef CORE_MATRIX_BASE_MATRIX_DISTRIBUTED_BASE_HPP_
#define CORE_MATRIX_BASE_MATRIX_DISTRIBUTED_BASE_HPP_

#include "../core_matrix_traits.hpp"

namespace core{
    
template<typename derived_type>
class MatrixDistributedBase
  : private core::details::CrtpBase<MatrixDistributedBase<derived_type>>
{
private:
  using traits_t = details::traits<derived_type>;

  using sc_t = typename traits_t::scalar_t;
  using der_t = typename traits_t::derived_t;
  using LO_t = typename traits_t::local_ordinal_t;
  using GO_t = typename traits_t::global_ordinal_t;
  using wrap_t = typename traits_t::wrapped_t;
  using row_map_t = typename traits_t::row_map_t;
  using col_map_t = typename traits_t::col_map_t;
  using comm_t =  typename traits_t::communicator_t;

  static_assert( details::traits<derived_type>::isDistributed==1,
  "OOPS: non-distributed matrix inheriting from sparse distributed base!");

public:
  LO_t localRows() const{
    return this->underlying().localRowsImpl();
  }

  LO_t localCols() const{
    return this->underlying().localColsImpl();
  }

  GO_t globalRows() const{
    return this->underlying().globalRowsImpl();
  }

  GO_t globalCols() const{
    return this->underlying().globalColsImpl();
  }
  
  row_map_t const & getRowDataMap() const{
    return this->underlying().getRowDataMapImpl();
  }
 
  col_map_t const & getColDataMap() const{
    return this->underlying().getColDataMapImpl();
  }

  bool isFillingCompleted() const{
    return this->underlying().isFillingCompletedImpl();
  }

  void fillingIsCompleted(){
    this->underlying().fillingIsCompletedImpl();
  }
  
private:
  friend derived_type;
  friend core::details::CrtpBase<MatrixDistributedBase<derived_type>>;

  MatrixDistributedBase() = default;
  ~MatrixDistributedBase() = default;
  
};//end class
  
} // end namespace core
#endif
