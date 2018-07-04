
#ifndef CORE_MATRIX_SPARSE_SERIAL_EIGEN_HPP_
#define CORE_MATRIX_SPARSE_SERIAL_EIGEN_HPP_

#include <Eigen/Core>
#include "../base/core_matrix_generic_base.hpp"
#include "../base/core_matrix_sparse_serial_base.hpp"
#include "../base/core_matrix_math_base.hpp"
#include "../../core_operators_base.hpp"
#include <stdexcept>

namespace core{

template <typename wrapped_type>
class matrix<wrapped_type,
	     typename
	     std::enable_if<
	       core::meta::is_matrixSparseSerialEigen<wrapped_type
					 	      >::value
	       >::type
	     >
  : public matrixGenericBase< matrix<wrapped_type> >,
    public matrixSparseSerialBase< matrix<wrapped_type> >,
    public matrixMathBase< matrix<wrapped_type> >,
    public arithmeticOperatorsBase< matrix<wrapped_type> >,
    public compoundAssignmentOperatorsBase< matrix<wrapped_type> >
{
private:
  using derived_t = matrix<wrapped_type>;
  using mytraits = details::traits<derived_t>;
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename mytraits::ordinal_t;
  using wrap_t = typename mytraits::wrapped_t;
  using der_t = typename mytraits::derived_t;
  
public: 
  matrix() = delete;

  explicit matrix(ord_t nrows, ord_t ncols) {
    this->resize(nrows, ncols);
    this->compress();
  }

  // row-major matrix constructor (U is just a trick to enable sfinae)
  template <typename U = ord_t,
	    typename std::enable_if<mytraits::isRowMajor==1,U>::type * = nullptr>
  explicit matrix(U nrows, U ncols, U nonZerosPerRow) {
    this->resize(nrows, ncols);
    if( nonZerosPerRow > ncols )
      throw std::runtime_error("SPARSE MATRIX CNTR: estimated nonzeros larger then size of cols");
    data_.reserve(Eigen::VectorXi::Constant(nrows,nonZerosPerRow));
    this->compress();
  }

  // col-major matrix constructor
  template <typename U = ord_t, 
	    typename std::enable_if<mytraits::isRowMajor==0,U>::type * = nullptr>
  explicit matrix(U nrows, U ncols, U nonZerosPerCol) {
    this->resize(nrows, ncols);
    if( nonZerosPerCol > nrows )
      throw std::runtime_error("SPARSE MATRIX CNTR: estimated nonzeros larger then size of rows");
    data_.reserve(Eigen::VectorXi::Constant(ncols,nonZerosPerCol));
    this->compress();
  }

  explicit matrix(const wrap_t & other) : data_(other){
    this->compress();
  }

  ~matrix() = default;
  
public: 
  // note here that we return by copy
  sc_t operator() (ord_t row, ord_t col) const{
    // eigen returns 0 if the item is zero
    return data_.coeff(row,col);
  }
  derived_t operator+(const derived_t & other) {
    assert(haveCompatibleDimensions(*this, other) );
    derived_t res(other.rows(), other.cols());
    *res.data() = this->data_ + *other.data();
    return res;
  }
  derived_t operator-(const derived_t & other) {
    assert(haveCompatibleDimensions(*this, other) );
    derived_t res(other.rows(), other.cols());
    *res.data() = this->data_ - (*other.data());
    return res;
  }
  derived_t operator*(const derived_t & other) {
    assert(haveCompatibleDimensions(*this, other) );
    derived_t res(other.rows(), other.cols());
    *res.data() = this->data_ * (*other.data());
    return res;
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
  //from generic base
  wrap_t const * dataImpl() const{
    return &data_;
  };
  wrap_t * dataImpl(){
    return &data_;
  };  

  //from sparse serial base
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
	      mytraits::isRowMajor==1,
	      U>::type * = nullptr>
  void insertValuesImpl(U rowInd, U numEntries,
			const sc_t * values, const U * indices){
    for (ord_t j=0; j<numEntries; ++j)
      data_.insert(rowInd,indices[j]) = values[j];
  }
  // inserting for column major storage
  template <typename U = ord_t, 
	    typename std::enable_if<
	      mytraits::isRowMajor==0,
	      U>::type * = nullptr>
  void insertValuesImpl(U colInd, U numEntries,
			const sc_t * values, const U * indices){
    for (ord_t i=0; i<numEntries; ++i)
      data_.insert(indices[i],colInd) = values[i];
  }

  // from math base
  void scaleImpl(sc_t & factor){
    data_.coeffs() *= factor;
  };  
  
private:
  bool haveCompatibleDimensions(const derived_t & m1,
				const derived_t & m2)const{
    return (m1.rows()==m2.rows() && m1.cols()==m2.cols()) ? true : false;
  }

  void compress(){
    data_.makeCompressed();
  }

private:
  friend matrixGenericBase< derived_t >;
  friend matrixSparseSerialBase< derived_t >;
  friend matrixMathBase< derived_t >;

private:
  wrap_t data_;
     
};//end class 
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
