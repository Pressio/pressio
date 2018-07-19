
#ifndef CORE_MATRIX_DENSE_SERIAL_EIGEN_HPP_
#define CORE_MATRIX_DENSE_SERIAL_EIGEN_HPP_

#include <Eigen/Core>
#include "../base/core_matrix_generic_base.hpp"
#include "../base/core_matrix_serial_base.hpp"
#include "../base/core_matrix_dense_serial_base.hpp"
#include "../../core_operators_base.hpp"

namespace core{

template <typename wrapped_type>
class Matrix<wrapped_type,
	     typename
	     std::enable_if<
	       core::meta::is_matrix_dense_serial_eigen<
		 wrapped_type>::value
	       >::type
	     >
  : public MatrixGenericBase< Matrix<wrapped_type> >,
    public MatrixSerialBase< Matrix<wrapped_type> >,
    public MatrixDenseSerialBase< Matrix<wrapped_type> >,
    public ArithmeticOperatorsBase< Matrix<wrapped_type>>,
    public CompoundAssignmentOperatorsBase< Matrix<wrapped_type>>
{

private:
  using derived_t = Matrix<wrapped_type>;
  using mytraits = typename details::traits<derived_t>;  
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename mytraits::ordinal_t;
  using wrap_t = typename mytraits::wrapped_t;
  using der_t = typename mytraits::derived_t;
  
public:
  Matrix() = default;

  explicit Matrix(ord_t nrows, ord_t ncols) {
    this->resize(nrows,ncols);
  }

  explicit Matrix(const wrap_t & other) : data_(other){}

  ~Matrix() = default;

public:
  sc_t & operator() (ord_t row, ord_t col){
    // check if we are withinbound 
    return data_(row,col);
  }

  sc_t const & operator() (ord_t row, ord_t col) const{
    // check if we are withinbound 
    return data_(row,col);
  }

  derived_t operator+(const derived_t & other) const{
    assert( other.rows() == this->rows() );
    assert( other.cols() == this->cols() );
    derived_t res(other.rows(), other.cols());
    *res.data() = this->data_ + *other.data();
    return res;
  }

  derived_t operator-(const derived_t & other) const{
    assert( other.rows() == this->rows() );
    assert( other.cols() == this->cols() );
    derived_t res(other.rows(), other.cols());
    *res.data() = this->data_ - (*other.data());
    return res;
  }
  
  derived_t operator*(const derived_t & other) const{
    assert( other.rows() == this->rows() );
    assert( other.cols() == this->cols() );
    derived_t res(other.rows(), other.cols());
    *res.data() = this->data_ * (*other.data());
    return res;
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
  //from generic base
  wrap_t * dataImpl(){
    return &data_;
  };  

  wrap_t const * dataImpl() const{
    return &data_;
  };

  void addToDiagonalImpl(sc_t value) {
    // check matrix is diagonal
    assert(this->rows()==this->cols());
    for (ord_t ir=0; ir<this->rows(); ir++)
      data_(ir,ir) += value;
  };
  
  // from serial base
  ord_t rowsImpl() const{
    return data_.rows();
  }

  ord_t colsImpl() const{
    return data_.cols();
  }

  void resizeImpl(ord_t nrows, ord_t ncols){
    static_assert(mytraits::isStatic == false,
    "Trying to resize a matrix wrapping a static Eigen matrix!");
    data_.resize(nrows, ncols);
    data_ = wrapped_type::Zero(nrows,ncols);
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
