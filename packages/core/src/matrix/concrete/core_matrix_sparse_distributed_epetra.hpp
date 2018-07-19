
#ifndef CORE_MATRIX_CONCRETE_MATRIX_SPARSE_DISTRIBUTED_EPETRA_HPP_
#define CORE_MATRIX_CONCRETE_MATRIX_SPARSE_DISTRIBUTED_EPETRA_HPP_

#include <Eigen/Core>
#include "../base/core_matrix_generic_base.hpp"
#include "../base/core_matrix_distributed_base.hpp"
#include "../base/core_matrix_sparse_distributed_base.hpp"
#include "../../core_operators_base.hpp"

namespace core{

template <typename wrapped_type>
class Matrix<wrapped_type,
	     typename
	     std::enable_if<
	       core::meta::is_matrix_sparse_distributed_epetra<
		 wrapped_type>::value
	       >::type
	     >
  : public MatrixGenericBase< Matrix<wrapped_type> >,
    public MatrixDistributedBase< Matrix<wrapped_type> >,
    public MatrixSparseDistributedBase< Matrix<wrapped_type> >
{

private:
  using derived_t = Matrix<wrapped_type>;
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
  Matrix() = delete;

  explicit Matrix(const row_map_t & rowMap,
		  LO_t NumEntriesPerRow,
		  bool StaticProfile=false)
    : data_(Copy, rowMap, NumEntriesPerRow, StaticProfile){}

  explicit Matrix(const wrapped_type & objin)
    : data_(objin){}
  
  ~Matrix() = default;

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

  // void addToDiagonalImpl(sc_t value) {
  //   // check matrix is diagonal
  //   assert(this->globalRows()==this->globalCols());
  // };
  
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
  friend MatrixGenericBase< derived_t >;
  friend MatrixSerialBase< derived_t >;
  friend MatrixDenseSerialBase< derived_t >;

private:
  wrap_t data_;
     
};//end class 
}//end namespace core 
#endif
