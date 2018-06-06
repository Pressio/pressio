
#ifndef CORE_VECTOR_SERIAL_EIGEN_HPP_
#define CORE_VECTOR_SERIAL_EIGEN_HPP_

#include "core_meta.hpp"
#include "./base/core_vector_generic_base.hpp"
#include "./base/core_vector_serial_base.hpp"
#include "./base/core_vector_math_base.hpp"



namespace core{
  
template <typename wrapped_type>
class vector<wrapped_type,
	     typename std::enable_if<core::meta::is_vectorEigen<wrapped_type>::value
				     >::type
	     >
  : public vectorGenericBase< vector<wrapped_type> >,
    public vectorSerialBase< vector<wrapped_type> >,
    public vectorMathBase< vector<wrapped_type> >
{
public:    
  using derived_t = vector<wrapped_type>;
  using sc_t = typename details::traits<derived_t>::scalar_t;
  using ord_t = typename details::traits<derived_t>::ordinal_t;
  using der_t = typename details::traits<derived_t>::derived_t;
  using wrap_t = typename details::traits<derived_t>::wrapped_t;

private:
  wrap_t data_;

public:
  vector(){}
  vector(const wrap_t & src) : data_(src){}
  ~vector(){}

  
  sc_t & operator [] (ord_t i){
    return data_(i);
  };
  sc_t const & operator [] (ord_t i) const{
    return data_(i);
  };  
  //-----------------------------------
  
  wrap_t const * dataImpl() const{
    return &data_;
  };

  wrap_t & dataImpl(){
    return data_;
  };

  void resizeImpl(size_t val){
    data_.resize(val);
  };

  size_t sizeImpl() const {
    return (data_.rows()==1) ? data_.cols() : data_.rows();
  };
  //-----------------------------------

  // template <typename wrapped_matrix_type, typename U>
  // typename std::enable_if<std::is_same<U,der_t>::value >::type
  // matMultiplyImpl(const matrix<wrapped_matrix_type> & matin, U& result) const
  // {
  //   auto & resVec = result.getNonConstRefToData();
  //   auto * matPtr = matin.view();
  //   assert( matPtr->cols() == this->size() );
  //   //assert( result.size() == this->size() );
  //   resVec = (*matPtr)*data_;
  // };

  // sc_t dotImpl(const der_t & b) const{
  //   // // what is this?
  //   // // dot product: <this,b>
  //   // sc_t res = 0.0;
  //   // for (size_t i=0; i<this->size(); i++)
  //   //   res += data_[i]*b[i];
  //   // return res;
  // };

  // template <typename op_t>
  // void applyOpImpl(op_t op, sc_t a1,
  // 		   sc_t a2, const der_t & vin){
  //   // // what is this?
  //   // // this = a1*this op a2*vin;
  //   // for (size_t i=0; i<this->size(); i++)
  //   //   data_[i] = op()( a1*data_[i], a2*vin[i] );
  // }
  
  // sc_t norm2Impl() const{
  //   // sc_t result = 0;
  //   // for (size_t i=0; i<this->size(); i++)
  //   //   result += data_[i]*data_[i];
  //   // return result;
  // };

    
};

    
}//end namespace core
#endif

