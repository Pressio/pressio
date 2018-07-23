
#ifndef CORE_MATRIX_CONCRETE_MATRIX_SPARSE_DISTRIBUTED_EPETRA_HPP_
#define CORE_MATRIX_CONCRETE_MATRIX_SPARSE_DISTRIBUTED_EPETRA_HPP_

#include <Eigen/Core>
#include "../base/core_matrix_generic_base.hpp"
#include "../base/core_matrix_distributed_base.hpp"
#include "../base/core_matrix_sparse_distributed_base.hpp"
#include "../base/core_matrix_sparse_distributed_trilinos.hpp"
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
    public MatrixSparseDistributedBase< Matrix<wrapped_type> >,
    public MatrixSparseDistributedTrilinos< Matrix<wrapped_type> >
{

private:
  using derived_t = Matrix<wrapped_type>;
  using traits_t = details::traits<derived_t>;

  using sc_t = typename traits_t::scalar_t;
  using LO_t = typename traits_t::local_ordinal_t;
  using GO_t = typename traits_t::global_ordinal_t;
  using comm_t =  typename traits_t::communicator_t;
  using wrap_t = typename traits_t::wrapped_t;
  using row_map_t = typename traits_t::row_map_t;
  using col_map_t = typename traits_t::col_map_t;
  using range_map_t = typename traits_t::range_map_t;
  using domain_map_t = typename traits_t::domain_map_t;
  using crs_graph_t = typename traits_t::crs_graph_t;
  
public:
  Matrix() = delete;

  explicit Matrix(const row_map_t & rowMap,
		  LO_t NumEntriesPerRow,
		  bool StaticProfile=false)
    : data_(Copy, rowMap, NumEntriesPerRow, StaticProfile){}

  explicit Matrix(const wrap_t & objin)
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

  comm_t const & commCRef() const{
    return data_.Comm();
  }
  
  //--------------------------------
  //from sparse distributed trilinos
  //--------------------------------
  bool isFillingCompletedImpl() const{
    return data_.Filled();
  }

  void fillingIsCompletedImpl(){
    data_.FillComplete();
  }

  row_map_t const & getRowDataMapImpl() const{
    return data_.RowMap();
  }
 
  col_map_t const & getColDataMapImpl() const{
    return data_.ColMap();
  }

  range_map_t const & getRangeDataMapImpl() const{
    return data_.RangeMap();
  }
 
  domain_map_t const & getDomainDataMapImpl() const{
    return data_.DomainMap();
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
  friend MatrixDistributedBase< derived_t >;
  friend MatrixSparseDistributedBase< derived_t >;
  friend MatrixSparseDistributedTrilinos< derived_t >;

private:
  wrap_t data_;
     
};//end class 
}//end namespace core 
#endif
