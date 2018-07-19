
#ifndef CORE_VECTOR_CONCRETE_VECTOR_STDLIB_HPP_
#define CORE_VECTOR_CONCRETE_VECTOR_STDLIB_HPP_

#include "../base/core_vector_generic_base.hpp"
#include "../base/core_vector_serial_base.hpp"
#include "../base/core_vector_math_base.hpp"

namespace core{

template <typename wrapped_type>
class Vector<wrapped_type,
	     typename std::enable_if<
	       core::meta::is_vector_stdlib<wrapped_type>::value
	       >::type
	     >
  : public VectorGenericBase< Vector<wrapped_type> >,
    public VectorSerialBase< Vector<wrapped_type> >,
    public VectorMathBase< Vector<wrapped_type> >,
    public ArithmeticOperatorsBase< Vector<wrapped_type> >,
    public CompoundAssignmentOperatorsBase< Vector<wrapped_type> >
{
private:
  using derived_t = Vector<wrapped_type>;
  using sc_t = typename details::traits<derived_t>::scalar_t;
  using der_t = typename details::traits<derived_t>::derived_t;
  using wrap_t = typename details::traits<derived_t>::wrapped_t;
  using ord_t = typename details::traits<derived_t>::ordinal_t;

public:
  Vector() = default;
  explicit Vector(ord_t insize,
		  sc_t value = static_cast<sc_t>(0) ){
    this->resize(insize, value);
  }
  explicit Vector(const std::vector<sc_t> & src) : data_(src){}
  ~Vector(){}

public:
  sc_t & operator [] (ord_t i){
    return data_[i];
  };

  sc_t const & operator [] (ord_t i) const{
    return data_[i];
  };  

  derived_t operator+(const derived_t & other) const{
    derived_t res(other.size());
    std::transform(this->data_.begin(), this->data_.end(),
		   other.data()->begin(), res.data()->begin(),
		   std::plus<sc_t>());
    return res;
  }

  derived_t operator-(const derived_t & other) const{
    derived_t res(other.size()); 
    std::transform(this->data_.begin(), this->data_.end(),
		   other.data()->begin(), res.data()->begin(),
		   std::minus<sc_t>());
    return res;
  }

  derived_t operator*(const derived_t & other) const{
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
  //----------------
  //from generic base
  //----------------
  wrap_t const * dataImpl() const{
    return &data_;
  }
  wrap_t * dataImpl(){
    return data_;
  }
  void putScalarImpl(sc_t value) {
    for (ord_t i=0; i<this->size(); i++)
      data_[i] = value;
  }    
  void setZeroImpl() {
    this->putScalarImpl( static_cast<sc_t>(0) );
  }
  bool emptyImpl() const{
    return data_.empty();
  }

  //----------------
  //from serial base
  //----------------
  ord_t sizeImpl() const {
    return data_.size();
  }
  void resizeImpl(ord_t val) {
    data_.resize(val);
  }

  //----------------
  //from math base
  //----------------
  template<typename op_t>
  void inPlaceOpImpl(op_t op, sc_t a1, sc_t a2, const der_t & other){
    // this = a1*this op a2*other;
    for (ord_t i=0; i<this->size(); i++)
      data_[i] = op()( a1*data_[i], a2*other[i] );
  }
  void scaleImpl(sc_t & factor){
    for (ord_t i=0; i<this->size(); i++)
      data_[i] *= factor;
  }
  void norm1Impl(sc_t & result) const {
    result = static_cast<sc_t>(0);
    for (decltype(this->size()) i=0; i<this->size(); i++)
      result += std::abs(data_[i]);
  }
  void norm2Impl(sc_t & res) const {
    res = static_cast<sc_t>(0);
    for (decltype(this->size()) i=0; i<this->size(); i++)
      res += data_[i]*data_[i];
    res = std::sqrt(res);
  }
  void normInfImpl(sc_t & res) const {
    res = std::abs(data_[0]);
    for (decltype(this->size()) i=1; i<this->size(); i++){
      sc_t currVal = std::abs(data_[i]);
      if(currVal>res)
	res = currVal;
    }
  }
  void minValueImpl(sc_t & result) const {
    result = std::max_element(data_.begin(), data_.end());
  }
  void maxValueImpl(sc_t & result) const {
    result = std::min_element(data_.begin(), data_.end());
  }

private:
  friend VectorGenericBase< derived_t >;
  friend VectorSerialBase< derived_t >;
  friend VectorMathBase< derived_t >;
  
private:
  std::vector<sc_t> data_;
  
};//end class  
}//end namespace core  
#endif
