
#ifndef CORE_MATRIX_DENSE_SERIAL_EIGEN_HPP_
#define CORE_MATRIX_DENSE_SERIAL_EIGEN_HPP_

#include <Eigen/Core>
#include "./base/core_matrix_generic_base.hpp"
#include "./base/core_matrix_dense_serial_base.hpp"
#include "../core_operators_base.hpp"


namespace core{

template <typename wrapped_type>
class matrix<wrapped_type,
	     typename
	     std::enable_if<
	       core::meta::is_matrixDenseSerialEigen<wrapped_type
						     >::value
	       >::type
	     >
  : public matrixGenericBase< matrix<wrapped_type> >,
    public matrixDenseSerialBase< matrix<wrapped_type> >,
    // maybe move operators inheritance to one of the bases
    public arithmeticOperatorsBase<matrix<wrapped_type>>,
    public compoundAssignmentOperatorsBase<matrix<wrapped_type>>
{
public:
  using derived_t = matrix<wrapped_type>;
  using mytraits = typename details::traits<derived_t>;  
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename mytraits::ordinal_t;
  using wrap_t = typename mytraits::wrapped_t;
  using der_t = typename mytraits::derived_t;
  
public:
  matrix() = default;

  matrix(ord_t nrows, ord_t ncols) {    
    // need to check that the wrapped type is NOT a static matrix from Eigen
    // otherwise we cannot resizee a static object.
    static_assert(mytraits::isStatic == false,
		  "You are trying to resize a matrix wrapping a static Eigen matrix!");
    //   this->resize(nrows, ncols);
   this->data_ = wrapped_type::Zero(nrows,ncols);
  }

  matrix(const wrap_t & other) : data_(other){}

  ~matrix() = default;

public:
  sc_t & operator() (ord_t row, ord_t col){
    // check if we are withinbound 
    return data_(row,col);
  }

  sc_t const & operator() (ord_t row, ord_t col) const{
    // check if we are withinbound 
    return data_(row,col);
  }

  derived_t operator+(const derived_t & other) {
    assert( other.rows() == this->rows() );
    assert( other.cols() == this->cols() );
    derived_t res(other.rows(), other.cols());
    *res.data() = this->data_ + *other.data();
    return res;
  }

  derived_t operator-(const derived_t & other) {
    assert( other.rows() == this->rows() );
    assert( other.cols() == this->cols() );
    derived_t res(other.rows(), other.cols());
    *res.data() = this->data_ - (*other.data());
    return res;
  }
  
  derived_t operator*(const derived_t & other) {
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
  wrap_t const * dataImpl() const{
    return &data_;
  };

  wrap_t * dataImpl(){
    return &data_;
  };
  
  ord_t rowsImpl() const{
    return data_.rows();
  }

  ord_t colsImpl() const{
    return data_.cols();
  }

  void resizeImpl(ord_t nrows, ord_t ncols){
    data_.resize(nrows, ncols);
    //need to check if the wrapped matrix is static size,
    //if so, we cannot resize it
    //    data_ = wrap_t::Zero(nrows,ncols);
  }

private:
  friend matrixGenericBase< derived_t >;
  friend matrixDenseSerialBase< derived_t >;

private:
  wrap_t data_;
     
};
 
}//end namespace core 
#endif







// void transposeImpl(derived_t & result) const{
//   result.getNonConstRefToData() = data_.transpose();
// }

// sc_t dotImpl(const der_t & b) const{
//   // what is this?
//   // dot product: <this,b>
//   sc_t res = 0.0;
//   for (size_t i=0; i<this->size(); i++)
//     res += data_[i]*b[i];
//   return res;
// };

// template <typename op_t>
// void applyOpImpl(op_t op, sc_t a1,
// 		   sc_t a2, const der_t & vin){
//   // what is this?
//   // this = a1*this op a2*vin;
//   for (size_t i=0; i<this->size(); i++)
//     data_[i] = op()( a1*data_[i], a2*vin[i] );
// }

// sc_t norm2Impl() const{
//   sc_t result = 0;
//   for (size_t i=0; i<this->size(); i++)
//     result += data_[i]*data_[i];
//   return result;
// };

