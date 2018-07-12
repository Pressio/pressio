
#ifndef CORE_MATRIX_SPARSE_DISTRIBUTED_EPETRA_HPP_
#define CORE_MATRIX_SPARSE_DISTRIBUTED_EPETRA_HPP_

#include <Eigen/Core>
#include "../base/core_matrix_generic_base.hpp"
#include "../base/core_matrix_distributed_base.hpp"
#include "../base/core_matrix_sparse_distributed_base.hpp"
#include "../../core_operators_base.hpp"

namespace core{

template <typename wrapped_type>
class matrix<wrapped_type,
	     typename
	     std::enable_if<
	       core::meta::is_matrixSparseDistributedEpetra<
		 wrapped_type>::value
	       >::type
	     >
  : public matrixGenericBase< matrix<wrapped_type> >,
    public matrixDistributedBase< matrix<wrapped_type> >,
    public matrixSparseDistributedBase< matrix<wrapped_type> >
{

private:
  using derived_t = matrix<wrapped_type>;
  using traits_t = details::traits<derived_t>;

  using sc_t = typename traits_t::scalar_t;
  using der_t = typename traits_t::derived_t;
  using LO_t = typename traits_t::local_ordinal_t;
  using GO_t = typename traits_t::global_ordinal_t;
  using row_map_t = typename traits_t::row_map_t;
  using col_map_t = typename traits_t::col_map_t;
  using comm_t =  typename traits_t::communicator_t;
  using wrap_t = typename traits_t::wrapped_t;
  
public:
  matrix() = delete;

  explicit matrix(const row_map_t & rowMap,
		  LO_t NumEntriesPerRow,
		  bool StaticProfile=false)
    : data_(Copy, rowMap, NumEntriesPerRow, StaticProfile){}

  explicit matrix(const wrapped_type & objin)
    : data_(objin){}
  
  ~matrix() = default;

private:
  //----------------
  //from general base
  //----------------
  wrap_t const * dataImpl() const{
    return &data_;
  }
  wrap_t * dataImpl(){
    return &data_;
  }

  //----------------------
  //from distributed base
  //----------------------
  LO_t localRowsImpl() const{
    return data_.NunMyRows();
  }

  LO_t localColsImpl() const{
    return data_.NunMyCols();
  }

  GO_t globalRowsImpl() const{
    return data_.NumGlobalRows();
  }

  GO_t globalColsImpl() const{
    return data_.NumGlobalCols();
  }
  
  row_map_t const & getRowDataMapImpl() const{
    return data_.RowMap();
  }
 
  col_map_t const & getColDataMapImpl() const{
    return data_.ColMap();
  }

  bool isFillingCompleted() const{
    return data_.Filled();
  }

  void fillingIsCompleted(){
    data_.FillComplete();
  }

  //----------------------------
  //from sparse distributed base
  //----------------------------
  void insertGlobalValuesImpl(GO_t targetRow,
			      LO_t numEntries,
			      const sc_t * values,
			      const GO_t * indices)
  {
    data_.InsertGlobalValues(targetRow, numEntries, values, indices);
  }
  
private:
  friend matrixGenericBase< derived_t >;
  friend matrixSerialBase< derived_t >;
  friend matrixDenseSerialBase< derived_t >;

private:
  wrap_t data_;
     
};//end class 
}//end namespace core 
#endif
