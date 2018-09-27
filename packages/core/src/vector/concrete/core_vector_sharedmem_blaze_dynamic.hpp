
#ifndef CORE_VECTOR_CONCRETE_VECTOR_SHAREDMEM_BLAZE_DYNAMIC_HPP_
#define CORE_VECTOR_CONCRETE_VECTOR_SHAREDMEM_BLAZE_DYNAMIC_HPP_

#include "../../shared_base/core_container_base.hpp"
#include "../../shared_base/core_operators_base.hpp"
#include "../../shared_base/core_container_resizable_base.hpp"

#include "../base/core_vector_sharedmem_base.hpp"
#include "../base/core_vector_math_base.hpp"

namespace rompp{
namespace core{
  
template <typename wrapped_type>
class Vector<wrapped_type,
    core::meta::enable_if_t<
      core::meta::is_dynamic_vector_blaze<wrapped_type>::value
      >
    >
  : public ContainerBase< Vector<wrapped_type>, wrapped_type >,
    public VectorSharedMemBase< Vector<wrapped_type> >,
    public VectorMathBase< Vector<wrapped_type> >,
    public CompoundAssignmentOperatorsBase<Vector<wrapped_type>>,
    public ContainerResizableBase<Vector<wrapped_type>, 1>{

  using this_t = Vector<wrapped_type>;
  using mytraits = typename details::traits<this_t>;  
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename  mytraits::ordinal_t;
  using wrap_t = typename mytraits::wrapped_t;

public:  
  Vector() = default;
  ~Vector() = default;
 
  explicit Vector(ord_t insize){
    this->resize(insize);
  }

  explicit Vector(const sc_t * src)
    : data_(src){}

  explicit Vector(const wrap_t & src)
    : data_(src){}

  Vector(this_t const & other)
    : data_(*other.data()){
  }

  // constructor from any expression, force evaluation
  template <typename T,
  	    core::meta::enable_if_t<
  	      T::is_vector_expression> * = nullptr>
  explicit Vector(const T & expr){
    this->resize(expr.size());
    for (ord_t i = 0; i != expr.size(); ++i)
      data_[i] = expr[i];
  }

  // assignment from any expression, force evaluation
  template <typename T,
  	    core::meta::enable_if_t<
  	      T::is_vector_expression> * = nullptr>
  this_t & operator=(const T & expr){
    if(this->size() != expr.size())
      this->resize(expr.size());
    for (ord_t i = 0; i != expr.size(); ++i)
      data_[i] = expr[i];
    return *this;
  }

  // assignment 
  template <typename T,
  	    core::meta::enable_if_t<
  	      std::is_same<T,this_t>::value> * = nullptr>
  this_t & operator=(const T & other){
    data_ = *other.data();
    return *this;
  }
  
  
public:
  sc_t & operator [] (ord_t i){
    return data_[i];
  };

  sc_t const & operator [] (ord_t i) const{
    return data_[i];
  };  

  this_t & operator+=(const this_t & other) {
    assert( other.size() == this->size() );
    this->data_ += *other.data();
    return *this;
  }
  
  this_t & operator-=(const this_t & other) {
    assert( other.size() == this->size() );
    this->data_ -= *other.data();
    return *this;
  }

private:

  wrap_t const * dataImpl() const{
    return &data_;
  }

  wrap_t * dataImpl(){
    return &data_;
  }

  void putScalarImpl(sc_t value) {
    for( typename
	  blaze::DynamicVector<ord_t>::Iterator it=data_.begin();
	 it!=data_.end(); ++it ) {
      it->value() = value;
    }    
  }

  void setZeroImpl() {
    this->putScalarImpl( static_cast<sc_t>(0) );
  }

  // bool emptyImpl() const{
  //   return this->size()==0 ? true : false;
  // }

  ord_t sizeImpl() const {
    return data_.size();
  }

  void resizeImpl(ord_t val){
    data_.resize(val);
  }
  
  void scaleImpl(sc_t & factor){
    // this = factor * this;
    for(typename
	 blaze::DynamicVector<ord_t>::Iterator it=data_.begin();
	 it!=data_.end(); ++it ) {
      it->value() *= factor;
    }
  }

  void norm1Impl(sc_t & result) const {
    result = blaze::l1Norm(data_);
  }

  void norm2Impl(sc_t & res) const {
    res = blaze::l2Norm(data_);
  }

  void normInfImpl(sc_t & res) const {
    res = blaze::maxNorm(data_);
  }

  void minValueImpl(sc_t & result) const {
    result = blaze::min(data_);
  }

  void maxValueImpl(sc_t & result) const {
    result = blaze::max(data_);
  }
  
private:
  friend ContainerBase< this_t, wrapped_type >;
  friend VectorSharedMemBase< this_t >;
  friend VectorMathBase< this_t >;  
  friend CompoundAssignmentOperatorsBase< this_t >;  
  friend ContainerResizableBase<this_t, 1>;

private:
  wrap_t data_;
 
};//end class    
}//end namespace core
}//end namespace rompp
#endif
