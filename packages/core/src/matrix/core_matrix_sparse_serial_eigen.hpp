
#ifndef CORE_MATRIX_SPARSE_SERIAL_EIGEN_HPP_
#define CORE_MATRIX_SPARSE_SERIAL_EIGEN_HPP_

#include <Eigen/Core>
#include "./base/core_matrix_generic_base.hpp"
#include "./base/core_matrix_sparse_serial_base.hpp"
#include "../core_operators_base.hpp"


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
    public matrixSparseSerialBase< matrix<wrapped_type> >
    // // maybe move operators inheritance to one of the bases
    // public arithmeticOperatorsBase<matrix<wrapped_type>>,
    // public compoundAssignmentOperatorsBase<matrix<wrapped_type>>
{
public:
  using derived_t = matrix<wrapped_type>;
  using sc_t = typename details::traits<derived_t>::scalar_t;
  using ord_t = typename details::traits<derived_t>::ordinal_t;
  using wrap_t = typename details::traits<derived_t>::wrapped_t;
  using der_t = typename details::traits<derived_t>::derived_t;
  
public:
  matrix() = default;
  matrix(ord_t nrows, ord_t ncols) {
    //   this->resize(nrows, ncols);
  }
  matrix(const wrap_t & other) : data_(other){}
  ~matrix() = default;
  
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
    // data_.resize(nrows, ncols);
    // //need to check if the wrapped matrix is static size,
    // //if so, we cannot resize it
  }

private:
  friend matrixGenericBase< derived_t >;
  friend matrixSparseSerialBase< derived_t >;

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

