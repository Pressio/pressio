
#ifndef CORE_VECTOR_STDLIB_HPP_
#define CORE_VECTOR_STDLIB_HPP_

#include "../base/core_vector_generic_base.hpp"
#include "../base/core_vector_serial_base.hpp"
#include "../base/core_vector_math_base.hpp"

namespace core{

template <typename wrapped_type>
class vector<wrapped_type,
	     typename std::enable_if<
	       core::meta::is_vectorStdLib<wrapped_type>::value
	       >::type
	     >
  : public vectorGenericBase< vector<wrapped_type> >,
    public vectorSerialBase< vector<wrapped_type> >,
    public vectorMathBase< vector<wrapped_type> >,
    public arithmeticOperatorsBase< vector<wrapped_type> >,
    public compoundAssignmentOperatorsBase< vector<wrapped_type> >
{
private:
  using derived_t = vector<wrapped_type>;
  using sc_t = typename details::traits<derived_t>::scalar_t;
  using der_t = typename details::traits<derived_t>::derived_t;
  using wrap_t = typename details::traits<derived_t>::wrapped_t;
  using ord_t = typename details::traits<derived_t>::ordinal_t;

public:
  vector() = default;
  vector(ord_t insize,
	 sc_t value = static_cast<sc_t>(0) ){
    this->resize(insize, value);
  }
  vector(const std::vector<sc_t> & src) : data_(src){}
  ~vector(){}

public:
  sc_t & operator [] (ord_t i){
    return data_[i];
  };

  sc_t const & operator [] (ord_t i) const{
    return data_[i];
  };  

  derived_t operator+(const derived_t & other) {
    derived_t res(other.size());
    std::transform(this->data_.begin(), this->data_.end(),
		   other.data()->begin(), res.data()->begin(),
		   std::plus<sc_t>());
    return res;
  }

  derived_t operator-(const derived_t & other) {
    derived_t res(other.size()); 
    std::transform(this->data_.begin(), this->data_.end(),
		   other.data()->begin(), res.data()->begin(),
		   std::minus<sc_t>());
    return res;
  }

  derived_t operator*(const derived_t & other) {
    derived_t res(other.size()); 
    std::transform(this->data_.begin(), this->data_.end(),
		   other.data()->begin(), res.data()->begin(),
		   std::multiplies<sc_t>());
    return res;
  }
  
  derived_t & operator+=(const derived_t & other) {
    std::transform(this->data_.begin(), this->data_->end(),
		   other.data()->begin(), this->data_.begin(),
		   std::plus<sc_t>());
    return *this;
  }
  
  derived_t & operator-=(const derived_t & other) {
    std::transform(this->data_.begin(), this->data_->end(),
		   other.data()->begin(), this->data_.begin(),
		   std::minus<sc_t>());
    return *this;
  }


private:
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
  bool emptyImpl() const{
    return data_.empty();
  };

private:
  friend vectorGenericBase< derived_t >;
  friend vectorSerialBase< derived_t >;
  friend vectorMathBase< derived_t >;
  
private:
  std::vector<sc_t> data_;
  
};//end class
  
}//end namespace core  
#endif


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

