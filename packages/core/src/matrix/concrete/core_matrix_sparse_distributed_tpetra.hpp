
#ifdef HAVE_TRILINOS
#ifndef CORE_MATRIX_CONCRETE_MATRIX_SPARSE_DISTRIBUTED_TPETRA_HPP_
#define CORE_MATRIX_CONCRETE_MATRIX_SPARSE_DISTRIBUTED_TPETRA_HPP_

#include "../../shared_base/core_container_base.hpp"
#include "../../shared_base/core_container_distributed_mpi_base.hpp"
#include "../../shared_base/core_container_distributed_trilinos_base.hpp"

#include "../base/core_matrix_base.hpp"
#include "../base/core_matrix_sparse_base.hpp"
#include "../base/core_matrix_distributed_base.hpp"
#include "../base/core_matrix_sparse_distributed_base.hpp"
#include "../base/core_matrix_sparse_distributed_trilinos_base.hpp"

namespace rompp{ namespace core{

template <typename wrapped_type>
class Matrix<wrapped_type,
	     core::meta::enable_if_t<
	       core::meta::is_matrix_sparse_distributed_tpetra<
		 wrapped_type>::value>
	     >
  : public ContainerBase< Matrix<wrapped_type>, wrapped_type >,
    public MatrixDistributedBase< Matrix<wrapped_type> >,
    public MatrixSparseDistributedTrilinosBase< Matrix<wrapped_type> >,
    public ContainerDistributedMpiBase< Matrix<wrapped_type>,
      typename details::traits<Matrix<wrapped_type>>::communicator_t>/*
    public MatrixBase< Matrix<wrapped_type> >,
    public MatrixSparseBase< Matrix<wrapped_type> >,
    public MatrixSparseDistributedBase< Matrix<wrapped_type> >,
    */
{

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

public:
  Matrix() = delete;

  Matrix(Teuchos::RCP<const row_map_t>  rowMap,
	 LO_t NumEntriesPerRow)
    : data_(rowMap, NumEntriesPerRow){}

  Matrix(const derived_t &) = delete;

  // explicit Matrix(const wrap_t & objin)
  //   : data_(objin.getRowMap(),
  // 	    objin.getColMap(),
  // 	    objin.getGlobalMaxNumRowEntries()){
  //   assert(objin.isFillComplete());
  // }

  ~Matrix() = default;

private:
  wrap_t const * dataImpl() const{
    return &data_;
  }

  wrap_t * dataImpl(){
    return &data_;
  }

  wrap_t dataCpImpl(){
    return data_;
  }

  LO_t localRowsImpl() const{
    return data_.getNodeNumRows();
  }

  LO_t localColsImpl() const{
    return data_.getNodeNumCols();
  }

  GO_t globalRowsImpl() const{
    return data_.getGlobalNumRows();
  }

  GO_t globalColsImpl() const{
    return data_.getGlobalNumCols();
  }

  comm_t commImpl() const{
    return data_.getComm();
  }

  bool isFillingCompletedImpl() const{
    return data_.isFillComplete();
  }

  void fillingIsCompletedImpl(){
    data_.fillComplete();
  }

  void fillingIsCompletedImpl(Teuchos::RCP<const domain_map_t> const & dmap,
                              Teuchos::RCP<const range_map_t> const & rmap){
    // this is needed for rectangular matrices
    data_.fillComplete(dmap, rmap);
  }

private:
  friend ContainerBase< derived_t, wrapped_type >;
  friend MatrixDistributedBase< derived_t >;
  friend MatrixSparseDistributedTrilinosBase< derived_t >;
  friend ContainerDistributedMpiBase< derived_t, comm_t >;
  // friend MatrixBase< derived_t >;
  // friend MatrixSparseBase< derived_t >;
  // friend MatrixSparseDistributedBase< derived_t >;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace rompp::core
#endif
#endif
