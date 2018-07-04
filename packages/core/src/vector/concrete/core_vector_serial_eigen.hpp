
#ifndef CORE_VECTOR_SERIAL_EIGEN_HPP_
#define CORE_VECTOR_SERIAL_EIGEN_HPP_

#include "../../meta/core_meta_basic.hpp"
#include "../../meta/core_meta_detect_operators.hpp"
#include "../../meta/core_meta_detect_typedefs.hpp"
#include "../base/core_vector_generic_base.hpp"
#include "../base/core_vector_serial_base.hpp"
#include "../base/core_vector_math_base.hpp"

namespace core{
  
template <typename wrapped_type>
class vector<wrapped_type,
	     typename std::enable_if<
	       core::meta::is_vectorEigen<wrapped_type>::value
	       >::type
	     >
  : public vectorGenericBase< vector<wrapped_type> >,
    public vectorSerialBase< vector<wrapped_type> >,
    public vectorMathBase< vector<wrapped_type> >,
    public arithmeticOperatorsBase<vector<wrapped_type>>,
    public compoundAssignmentOperatorsBase<vector<wrapped_type>>
{
private:
  using derived_t = vector<wrapped_type>;
  using mytraits = typename details::traits<derived_t>;  
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename  mytraits::ordinal_t;
  using der_t = typename mytraits::derived_t;
  using wrap_t = typename mytraits::wrapped_t;

public:
  vector() = default;

  explicit vector(ord_t insize){
    this->resize(insize);
  }

  explicit vector(const wrap_t & src) : data_(src){}

  ~vector(){}
  
public:
  sc_t & operator [] (ord_t i){
    //assert(!this->empty());
    return data_(i);
  };

  sc_t const & operator [] (ord_t i) const{
    //assert(!this->empty());
    return data_(i);
  };  

  derived_t operator+(const derived_t & other) {
    assert( other.size() == this->size() );
    derived_t res(other.size());
    *res.data() = this->data_ + *other.data();
    return res;
  }

  derived_t operator-(const derived_t & other) {
    assert( other.size() == this->size() );
    derived_t res(other.size());
    *res.data() = this->data_ - *other.data();
    return res;
  }
  
  derived_t operator*(const derived_t & other) {
    assert( other.size() == this->size() );
    derived_t res(other.size());
    for (size_t i=0; i<this->size(); i++)
      res[i] = this->data_(i) * other[i];
    return res;
  }
  
  derived_t & operator+=(const derived_t & other) {
    assert( other.size() == this->size() );
    this->data_ += *other.data();
    return *this;
  }
  
  derived_t & operator-=(const derived_t & other) {
    assert( other.size() == this->size() );
    this->data_ -= *other.data();
    return *this;
  }

private:
  //-----------------
  //from generic base
  //-----------------
  wrap_t const * dataImpl() const{
    return &data_;
  };
  wrap_t * dataImpl(){
    return &data_;
  };
  void putScalarImpl(sc_t value) {
    data_.setConstant(value);
  }    

  //-----------------
  //from serial base
  //-----------------
  size_t sizeImpl() const {
    return (data_.rows()==1) ? data_.cols() : data_.rows();
  };
  void resizeImpl(size_t val){
    // check that the wrapped type is NOT a static vector from Eigen
    // otherwise we cannot resizee a static vector.
    static_assert(mytraits::isStatic == false,
		  "You cannot resize a vector wrapping a STATIC Eigen vector!");
    data_.resize(val);
  };
  bool emptyImpl() const{
    return this->size()==0 ? true : false;
  };

  //----------------
  //from math base
  //----------------
  template<typename op_t>
  void inPlaceOpImpl(op_t op, sc_t a1, sc_t a2, const der_t & other){
    // this = a1*this op a2*other;
    for (decltype(this->size()) i=0; i<this->size(); i++)
      data_(i) = op()( a1*data_(i), a2*other[i] );
  }
  void scaleImpl(sc_t & factor){
    // this = factor * this;
    for (decltype(this->size()) i=0; i<this->size(); i++)
      data_(i) *= factor;
  };
  void norm1Impl(sc_t & result) const {
    result = static_cast<sc_t>(0);
    for (decltype(this->size()) i=0; i<this->size(); i++)
      result += std::abs(data_(i));
  };
  void norm2Impl(sc_t & res) const {
    res = static_cast<sc_t>(0);
    for (decltype(this->size()) i=0; i<this->size(); i++)
      res += data_(i)*data_(i);
    res = std::sqrt(res);
  };
  void normInfImpl(sc_t & res) const {
    res = std::abs(data_(0));
    for (decltype(this->size()) i=1; i<this->size(); i++){
      sc_t currVal = std::abs(data(i));
      if(currVal>res)
	res = currVal;
    }
  };
  void minValueImpl(sc_t & result) const {
    result = data_.minCoeff();
  };
  void maxValueImpl(sc_t & result) const {
    result = data_.maxCoeff();
  };
  
private:
  friend vectorGenericBase< derived_t >;
  friend vectorSerialBase< derived_t >;
  friend vectorMathBase< derived_t >;

private:
  wrap_t data_;
 
};//end class    
}//end namespace core
#endif
