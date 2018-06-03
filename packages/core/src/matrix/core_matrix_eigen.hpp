
#ifndef CORE_MATRIX_EIGEN_HPP_
#define CORE_MATRIX_EIGEN_HPP_

#include "core_matrix_generic_base.hpp"
#include <Eigen/Core>


namespace core{


template <typename wrapped_type>
class matrix<wrapped_type,
	     typename std::enable_if< core::meta::is_matrixEigen<wrapped_type>::value
				      >::type
	     >
  : public matrixGenericBase<matrix<wrapped_type> >
{
public:
  using derived_t = matrix<wrapped_type>;
  using sc_t = typename details::traits<derived_t>::scalar_t;
  using ord_t = typename details::traits<derived_t>::ordinal_t;
  using wrap_t = typename details::traits<derived_t>::wrapped_t;
  using der_t = typename details::traits<derived_t>::derived_t;
  static_assert(std::is_same<ord_t,int>::value, "ordinal_type should be int for this class");
  static_assert(details::traits<derived_t>::isEigen==1, "isEigen");
  //static_assert( std::is_same<sc_t, double>::value, "isEigen");
  
private:
  wrap_t data_;

public:
  matrix(const wrap_t & other) : data_(other)
  {}
  matrix(){
    std::cout << "default" << std::endl;
  }
  ~matrix(){}
  
  // sc_t & operator [] (ord_t i){
  //   return data_[i];
  // };
  // sc_t const & operator [] (ord_t i) const{
  //   return data_[i];
  // };  
  //-----------------------------------

  void transposeImpl(derived_t & result) const{
    result.getNonConstRefToData() = data_.transpose();
  }
  
  wrap_t const * viewImpl() const{
    return &data_;
  };

  wrap_t & getNonConstRefToDataImpl(){
    return data_;
  };
  
  size_t rows() const{
    return data_.rows();
  }

  size_t cols() const{
    return data_.cols();
  }

  void resizeImpl(size_t nrows, size_t ncols){
    data_.resize(nrows, ncols);
  }
  

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

  
};
  
}//end namespace core
  
#endif

