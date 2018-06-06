
#ifndef CORE_VECTOR_SERIAL_USERDEFINED_HPP_
#define CORE_VECTOR_SERIAL_USERDEFINED_HPP_

#include "./base/core_vector_generic_base.hpp"
#include "./base/core_vector_serial_base.hpp"
#include "./base/core_vector_math_base.hpp"


namespace core{

  
template <typename wrapped_type>
class vector<wrapped_type,
	     typename std::enable_if< !core::meta::is_stdlibVector<wrapped_type>::value &&
				      !core::meta::is_vectorEigen<wrapped_type>::value
				    >::type
	     >
  : public vectorGenericBase< vector<wrapped_type> >,
    public vectorSerialBase< vector<wrapped_type> >,
    public vectorMathBase< vector<wrapped_type> >
{
public:
  using derived_t = vector<wrapped_type>;
  using sc_t = typename details::traits<derived_t>::scalar_t;
  using der_t = typename details::traits<derived_t>::derived_t;
  using wrap_t = typename details::traits<derived_t>::wrapped_t;
  using ord_t = typename details::traits<derived_t>::ordinal_t;

private:
   wrap_t data_;

public:
  vector(){};
  vector(const wrap_t & obj) : data_(obj){};
  ~vector(){};


  sc_t & operator [] (ord_t i){
    return data_[i];
  };
  sc_t const & operator [] (ord_t i) const{
    return data_[i];
  };  
  //-----------------------------------

  wrap_t const * dataImpl() const{
    return &data_;
  };

  wrap_t * dataImpl(){
    return data_;
  };

  void resizeImpl(size_t val) {
    data_.resize(val);
  };

  size_t sizeImpl() const {
    return data_.size();
  };

};
  
  
}//end namespace core
#endif





  // sc_t dotImpl(const der_t & b) const{
  //   return data_.dot(b);
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
  //   return data_.norm2();
  // };
