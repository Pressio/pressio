
#ifndef CORE_MATRIX_EIGENSTATIC_HPP
#define CORE_MATRIX_EIGENSTATIC_HPP

#include "core_matrixBaseImpl.hpp"
#include <Eigen/Core>


namespace core{


  template <typename scalar_type,
	    int RowsAtCompileTime,
	    int ColsAtCompileTime
	    >
  class matrix< Eigen::Matrix<scalar_type,RowsAtCompileTime,ColsAtCompileTime>,
		scalar_type,
		int,
	        void,
		RowsAtCompileTime,
		ColsAtCompileTime,
		void,void,void>
	       // typename std::enable_if<std::is_same<RowsAtCompileTime,int>::value>::type
	       // >
    // : public matrixBaseImpl<matrix<Eigen::Matrix<scalar_type,RowsAtCompileTime,ColsAtCompileTime>,
    // 				   scalar_type,int,void,RowsAtCompileTime,ColsAtCompileTime
    // 				   >
    // 			    >
{
public:
  // using derived_t = matrix<Eigen::Matrix<scalar_type,RowsAtCompileTime,ColsAtCompileTime>,
  // 			   scalar_type,int,void,RowsAtCompileTime,ColsAtCompileTime,void,void,void>;
  // using sc_t = typename details::traits<derived_t>::scalar_t;
  // using der_t = typename details::traits<derived_t>::derived_t;
  // using wrap_t = typename details::traits<derived_t>::wrapped_t;
  // using ord_t = typename details::traits<derived_t>::ordinal_t;
  // static_assert(std::is_same<ord_t,int>::value, "ordinal_type should be int for this class");

private:
  //  Eigen::Matrix<sc_t,RowsAtCompileTime,ColsAtCompileTime> data_;

public:

  // sc_t & operator [] (ord_t i){
  //   return data_[i];
  // };
  // sc_t const & operator [] (ord_t i) const{
  //   return data_[i];
  // };  
  //-----------------------------------
  
  // wrap_t const * viewImpl() const{
  //   return &data_;
  // };
  // wrap_t & getNonConstRefToDataImpl(){
  //   return data_;
  // };

  // void resizeImpl(size_t val) {
  //   data_.resize(val);
  // };
  // //-----------------------------------

  // size_t sizeImpl() const {
  //   return data_.size();
  // };
  // //-----------------------------------
  
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



  // template <typename scalar_type,
  // 	    typename RowsAtCompileTime,
  // 	    typename ColsAtCompileTime>
  // using matrixStatic = matrix<Eigen::Matrix<scalar_type,RowsAtCompileTime,ColsAtCompileTime>,
  // 			      scalar_type,int,void,RowsAtCompileTime,ColsAtCompileTime,void,void,void
  // 			      >;
  // 					    // std::enable_if<std::is_same<RowsAtCompileTime,int>::value>::type,
  // 					    // std::enable_if<std::is_same<ColsAtCompileTime,int>::value>::type>,

  
}//end namespace core
  
#endif

