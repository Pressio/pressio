
#ifndef CORE_MATRIX_CONCRETE_MATRIX_SPARSE_SHAREDMEM_EIGEN_HPP_
#define CORE_MATRIX_CONCRETE_MATRIX_SPARSE_SHAREDMEM_EIGEN_HPP_

#include "../../shared_base/core_container_base.hpp"
#include "../../shared_base/core_container_resizable_base.hpp"
#include "../../shared_base/core_container_subscriptable_base.hpp"

#include "../base/core_matrix_base.hpp"
#include "../base/core_matrix_sparse_base.hpp"
#include "../base/core_matrix_sharedmem_base.hpp"
#include "../base/core_matrix_sparse_sharedmem_base.hpp"

namespace rompp{ namespace core{

template <typename wrapped_type>
class Matrix<wrapped_type,
	     core::meta::enable_if_t<
	       core::meta::is_matrix_sparse_sharedmem_eigen<
		 wrapped_type >::value>
	     >
  : public ContainerBase< Matrix<wrapped_type>, wrapped_type >,
    public MatrixBase< Matrix<wrapped_type> >,
    public MatrixSparseBase< Matrix<wrapped_type> >,
    public MatrixSharedMemBase< Matrix<wrapped_type> >,
    public MatrixSparseSharedMemBase< Matrix<wrapped_type> >,
    public ContainerResizableBase<Matrix<wrapped_type>, 2>
{

  using derived_t = Matrix<wrapped_type>;
  using mytraits = details::traits<derived_t>;
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename mytraits::ordinal_t;
  using wrap_t = typename mytraits::wrapped_t;
  using der_t = typename mytraits::derived_t;

public:
  Matrix() = delete;

  explicit Matrix(ord_t nrows, ord_t ncols) {
    this->resize(nrows, ncols);
    this->compress();
  }

  // row-major matrix constructor
  template <typename U = ord_t,
	    typename std::enable_if<
	      mytraits::is_row_major==1, U>::type * = nullptr>
  explicit Matrix(U nrows, U ncols, U nonZerosPerRow) {
    this->resize(nrows, ncols);
    if( nonZerosPerRow > ncols )
      throw std::runtime_error(
    "SPARSE MATRIX CNTR: estimated nonzeros larger then size of cols");
    data_.reserve(Eigen::VectorXi::Constant(nrows,nonZerosPerRow));
    this->compress();
  }

  // col-major matrix constructor
  template <typename U = ord_t,
	    typename std::enable_if<
	      mytraits::is_row_major==0, U>::type * = nullptr>
  explicit Matrix(U nrows, U ncols, U nonZerosPerCol) {
    this->resize(nrows, ncols);
    if( nonZerosPerCol > nrows )
      throw std::runtime_error(
    "SPARSE MATRIX CNTR: estimated nonzeros larger then size of rows");
    data_.reserve(Eigen::VectorXi::Constant(ncols,nonZerosPerCol));
    this->compress();
  }

  explicit Matrix(const wrap_t & other) : data_(other){
    this->compress();
  }

  ~Matrix() = default;

public:

  // note here that we return by copy
  sc_t operator() (ord_t row, ord_t col) const{
    // eigen returns 0 if the item is zero
    return data_.coeff(row,col);
  }

  derived_t & operator+=(const derived_t & other) {
    assert(haveCompatibleDimensions(*this, other) );
    this->data_ += *other.data();
    return *this;
  }

  derived_t & operator-=(const derived_t & other) {
    assert(haveCompatibleDimensions(*this, other) );
    this->data_ -= *other.data();
    return *this;
  }

private:
  wrap_t const * dataImpl() const{
    return &data_;
  };

  wrap_t * dataImpl(){
    return &data_;
  };

  void addToDiagonalImpl(sc_t value) {
    auto ide(data_);
    ide.setIdentity();
    ide.coeffs() *= value;
    data_ += ide;
  };

  // from sharedMem base
  ord_t rowsImpl() const{
    return data_.rows();
  }

  ord_t colsImpl() const{
    return data_.cols();
  }

  void resizeImpl(ord_t nrows, ord_t ncols){
    data_.resize(nrows, ncols);
  }

  ord_t nonZerosCountImpl()const{
    return data_.nonZeros();
  }

  void setIdentityImpl(){
    data_.setIdentity();
  }

  void setZeroImpl(){
    data_.setZero();
  }

  // inserting for row major storage
  template <typename U = ord_t,
	    typename std::enable_if<
	      mytraits::is_row_major==1,
	      U>::type * = nullptr>
  void insertValuesImpl(U rowInd, U numEntries,
			const sc_t * values, const U * indices){
    for (ord_t j=0; j<numEntries; ++j)
      data_.insert(rowInd,indices[j]) = values[j];
  }

  // inserting for column major storage
  template <typename U = ord_t,
	    typename std::enable_if<
	      mytraits::is_row_major==0,
	      U>::type * = nullptr>
  void insertValuesImpl(U colInd, U numEntries,
			const sc_t * values, const U * indices){
    for (ord_t i=0; i<numEntries; ++i)
      data_.insert(indices[i],colInd) = values[i];
  }

  void scaleImpl(sc_t & factor){
    data_.coeffs() *= factor;
  };

  void compressImpl(){
    data_.makeCompressed();
  }

  bool isCompressedImpl() const{
    return data_.isCompressed();
  }

private:
  bool haveCompatibleDimensions(const derived_t & m1,
				const derived_t & m2)const{
    return (m1.rows()==m2.rows() &&
	    m1.cols()==m2.cols()) ? true : false;
  }

private:
  friend ContainerBase< derived_t, wrapped_type >;
  friend MatrixBase< derived_t >;
  friend MatrixSparseBase< derived_t >;
  friend MatrixSharedMemBase< derived_t >;
  friend MatrixSparseSharedMemBase< derived_t >;
  friend ContainerResizableBase<derived_t, 2>;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace rompp::core
#endif
