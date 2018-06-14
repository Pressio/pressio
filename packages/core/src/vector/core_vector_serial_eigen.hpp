
#ifndef CORE_VECTOR_SERIAL_EIGEN_HPP_
#define CORE_VECTOR_SERIAL_EIGEN_HPP_

#include "../meta/core_meta_basic.hpp"
#include "../meta/core_meta_detect_operators.hpp"
#include "../meta/core_meta_detect_typedefs.hpp"
#include "./base/core_vector_generic_base.hpp"
#include "./base/core_vector_serial_base.hpp"
#include "./base/core_vector_math_base.hpp"

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
    // maybe move operators inheritance to serial/generic base
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
  vector(ord_t insize){
    // need to check that the wrapped type is NOT a static vector from Eigen
    // otherwise we cannot resizee a static vector.
    static_assert(mytraits::isStatic == false,
		  "You are trying to resize a vector wrapping a static Eigen vector!");
    this->resize(insize);
  }
  vector(const wrap_t & src) : data_(src){}
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
  wrap_t const * dataImpl() const{
    return &data_;
  };
  wrap_t * dataImpl(){
    return &data_;
  };
  size_t sizeImpl() const {
    return (data_.rows()==1) ? data_.cols() : data_.rows();
  };
  void resizeImpl(size_t val){
    data_.resize(val);
  };
  bool emptyImpl() const{
    return this->size()==0 ? true : false;
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



