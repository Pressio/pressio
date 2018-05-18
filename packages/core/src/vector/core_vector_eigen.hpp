
#ifndef CORE_VECTOR_EIGEN_HPP
#define CORE_VECTOR_EIGEN_HPP

#include "core_vectorBaseImpl.hpp"
#include "core_vectorMathImpl.hpp"
#include "core_vectorSerImpl.hpp"
#include "meta.hpp"

#include <Eigen/Dense>


//*******************************
// eigen vector wrapper
//*******************************  

namespace core{
  
template <typename wrapped_type>
class vector<wrapped_type,
	     typename std::enable_if< core::meta::is_vectorEigen<wrapped_type>::value
				      >::type
	     >
  : public vectorBaseImpl<vector<wrapped_type> >,
    public vectorSerImpl<vector<wrapped_type> >,
    public vectorMathImpl<vector<wrapped_type> >
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

  sc_t & operator [] (ord_t i){
    return data_[i];
  };
  sc_t const & operator [] (ord_t i) const{
    return data_[i];
  };  
  //-----------------------------------
  
  wrap_t const * viewImpl() const{
    return &data_;
  };
  wrap_t & getNonConstRefToDataImpl(){
    return data_;
  };
  void resizeImpl(size_t val) {
    //data_.resize(val);
  };
  //-----------------------------------

  size_t sizeImpl() const {
    return 1;//data_.size();
  };
  //-----------------------------------
  
  sc_t dotImpl(const der_t & b) const{
    // // what is this?
    // // dot product: <this,b>
    // sc_t res = 0.0;
    // for (size_t i=0; i<this->size(); i++)
    //   res += data_[i]*b[i];
    // return res;
  };

  template <typename op_t>
  void applyOpImpl(op_t op, sc_t a1,
  		   sc_t a2, const der_t & vin){
    // // what is this?
    // // this = a1*this op a2*vin;
    // for (size_t i=0; i<this->size(); i++)
    //   data_[i] = op()( a1*data_[i], a2*vin[i] );
  }
  
  sc_t norm2Impl() const{
    // sc_t result = 0;
    // for (size_t i=0; i<this->size(); i++)
    //   result += data_[i]*data_[i];
    // return result;
  };

    
};

    
}//end namespace core
#endif

