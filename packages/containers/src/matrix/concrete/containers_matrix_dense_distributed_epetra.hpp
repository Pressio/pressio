
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_MATRIX_CONCRETE_MATRIX_DENSE_DISTRIBUTED_EPETRA_HPP_
#define CONTAINERS_MATRIX_CONCRETE_MATRIX_DENSE_DISTRIBUTED_EPETRA_HPP_

#include "../../shared_base/containers_container_base.hpp"
#include "../../shared_base/containers_container_distributed_mpi_base.hpp"
#include "../../shared_base/containers_container_distributed_trilinos_base.hpp"
#include "../../shared_base/containers_container_subscriptable_base.hpp"

#include "../base/containers_matrix_distributed_base.hpp"
#include "../base/containers_matrix_dense_distributed_base.hpp"
#include "../base/containers_matrix_base.hpp"

namespace rompp{ namespace containers{

template <typename wrapped_type>
class Matrix<wrapped_type,
	     ::rompp::mpl::enable_if_t<
	       containers::meta::is_dense_matrix_epetra<
		 wrapped_type>::value>
	     >
  : public ContainerBase< Matrix<wrapped_type>, wrapped_type >,
    public ContainerDistributedTrilinosBase<
      Matrix<wrapped_type>,
      typename details::traits<Matrix<wrapped_type>>::row_map_t >,
    public ContainerDistributedMpiBase<
      Matrix<wrapped_type>,
      typename details::traits<Matrix<wrapped_type>>::communicator_t >,
    public ContainerSubscriptable2DBase<
     Matrix<wrapped_type>,
     typename details::traits<Vector<wrapped_type>>::scalar_t,
     typename details::traits<Vector<wrapped_type>>::local_ordinal_t,
     typename details::traits<Vector<wrapped_type>>::global_ordinal_t>,
    public MatrixBase< Matrix<wrapped_type> >,
    public MatrixDistributedBase< Matrix<wrapped_type> >,
    public MatrixDenseDistributedBase< Matrix<wrapped_type> >{

  using this_t = Matrix<wrapped_type>;
  using traits_t = details::traits<this_t>;
  using sc_t = typename traits_t::scalar_t;
  using LO_t = typename traits_t::local_ordinal_t;
  using GO_t = typename traits_t::global_ordinal_t;
  using comm_t =  typename traits_t::communicator_t;
  using wrap_t = typename traits_t::wrapped_t;
  using row_map_t = typename traits_t::row_map_t;

public:

  Matrix() = delete;

  Matrix(const row_map_t & rowMap, GO_t ncols)
    : data_(rowMap, ncols){}

  explicit Matrix(const wrap_t & objin)
    : data_(objin){}

  ~Matrix() = default;

public:
  sc_t & operator()(LO_t irow, GO_t icol){
    assert(icol < this->globalCols() );
    assert(irow < this->localRows() );
    return data_[icol][irow];
  }

  sc_t const & operator()(LO_t irow, GO_t icol)const{
    assert(icol < this->globalCols() );
    assert(irow < this->localRows() );
    return data_[icol][irow];
  }

private:
  wrap_t const * dataImpl() const{
    return &data_;
  }

  wrap_t * dataImpl(){
    return &data_;
  }

  bool emptyImpl() const{
    return this->globalCols()==0 &&
      this->globalRows()==0 ? true : false;
  }

  void setZeroImpl() {
    data_.PutScalar(static_cast<sc_t>(0));
  }

  bool isDistributedGloballyImpl() const{
    return data_.DistributedGlobal();
  }

  void matchLayoutWithImpl(const this_t & other){
    data_.ReplaceMap( other.getDataMap() );
  }

  row_map_t const & getDataMapImpl() const{
    return data_.Map();
  }

  void replaceDataMapImpl(const row_map_t & mapObj){
    data_.ReplaceMap(mapObj);
  }

  LO_t localRowsImpl() const{
    return data_.MyLength();
  }

  LO_t localColsImpl() const{
    return data_.NumVectors();
  }

  GO_t globalRowsImpl() const{
    return data_.GlobalLength();
  }

  GO_t globalColsImpl() const{
    return data_.NumVectors();
  }

  comm_t const & commCRefImpl() const{
    return data_.Comm();
  }

  void replaceGlobalValueImpl(GO_t globalRowIndex,
			      GO_t globalColIndex,
			      sc_t value){
    data_.ReplaceGlobalValue(globalRowIndex,globalColIndex, value);
  }


private:
  friend ContainerBase< this_t, wrapped_type >;
  friend MatrixBase< this_t >;
  friend MatrixDistributedBase< this_t >;
  friend MatrixDenseDistributedBase< this_t >;
  friend ContainerDistributedMpiBase< this_t, comm_t >;
  friend ContainerDistributedTrilinosBase< this_t, row_map_t >;
  friend ContainerSubscriptable2DBase< this_t, sc_t, LO_t, GO_t>;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace rompp::containers
#endif
#endif
