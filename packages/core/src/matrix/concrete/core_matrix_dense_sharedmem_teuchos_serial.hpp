
#ifndef CORE_MATRIX_CONCRETE_MATRIX_DENSE_SHAREDMEM_TEUCHOS_HPP_
#define CORE_MATRIX_CONCRETE_MATRIX_DENSE_SHAREDMEM_TEUCHOS_HPP_

#include "../../shared_base/core_container_base.hpp"
#include "../../shared_base/core_container_resizable_base.hpp"
#include "../../shared_base/core_container_subscriptable_base.hpp"

#include "../base/core_matrix_base.hpp"
#include "../base/core_matrix_sharedmem_base.hpp"
#include "../base/core_matrix_dense_sharedmem_base.hpp"

namespace rompp{ namespace core{

template <typename wrapped_type>
class Matrix<wrapped_type,
	     core::meta::enable_if_t<
	       core::meta::is_teuchos_serial_dense_matrix<
		 wrapped_type>::value>
	     >
  : public ContainerBase< Matrix<wrapped_type>, wrapped_type >,
    public MatrixBase< Matrix<wrapped_type> >,
    public MatrixSharedMemBase< Matrix<wrapped_type> >,
    public MatrixDenseSharedMemBase< Matrix<wrapped_type> >,
    public ContainerSubscriptable2DBase<
     Matrix<wrapped_type>,
     typename details::traits<Matrix<wrapped_type>>::scalar_t,
     typename details::traits<Matrix<wrapped_type>>::ordinal_t>,
    public ContainerResizableBase<Matrix<wrapped_type>, 2>
{

  using derived_t = Matrix<wrapped_type>;
  using mytraits = typename details::traits<derived_t>;
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename mytraits::ordinal_t;
  using wrap_t = typename mytraits::wrapped_t;
  using der_t = typename mytraits::derived_t;

public:
  Matrix() = default;

  Matrix(ord_t nrows, ord_t ncols) {
    this->resize(nrows,ncols);
  }

  explicit Matrix(const wrap_t & other)
    : data_(other){}

  explicit Matrix(const sc_t * datain)
    : data_(datain){}

  ~Matrix() = default;

public:
  sc_t & operator() (ord_t row, ord_t col){
    return data_(row,col);
  }

  sc_t const & operator() (ord_t row, ord_t col) const{
    return data_(row,col);
  }

  derived_t & operator+=(const derived_t & other) {
    assert( other.rows() == this->rows() );
    assert( other.cols() == this->cols() );
    this->data_ += *other.data();
    return *this;
  }

  derived_t & operator-=(const derived_t & other) {
    assert( other.rows() == this->rows() );
    assert( other.cols() == this->cols() );
    this->data_ -= *other.data();
    return *this;
  }

private:
  wrap_t * dataImpl(){
    return &data_;
  };

  wrap_t const * dataImpl() const{
    return &data_;
  };

  void addToDiagonalImpl(sc_t value) {
    // careful with rectangular matrices
    auto nR = this->rows();
    auto nC = this->cols();
    auto loopLimit = (nR > nC) ? nC : nR;

    for (ord_t i=0; i<loopLimit; i++)
      data_(i,i) += value;
  };

  ord_t rowsImpl() const{
    return data_.numRows();
  }

  ord_t colsImpl() const{
    return data_.numCols();
  }

  void setZeroImpl() {
    data_.putScalar(static_cast<sc_t>(0));
  }

  void resizeImpl(ord_t nrows, ord_t ncols){
    //shape does a resize and sets entries to zero
    data_.shape(nrows, ncols);
  }

  void scaleImpl(sc_t value) {
    data_.scale(value);
  }

private:
  friend ContainerBase< derived_t, wrapped_type >;
  friend MatrixBase< derived_t >;
  friend MatrixSharedMemBase< derived_t >;
  friend MatrixDenseSharedMemBase< derived_t >;
  friend ContainerSubscriptable2DBase< derived_t, sc_t, ord_t>;
  friend ContainerResizableBase<derived_t, 2>;

private:
  wrap_t data_;

};//end class

}}//end namespace rompp::core
#endif
