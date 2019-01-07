
#ifdef HAVE_TRILINOS
#ifndef CORE_VECTOR_CONCRETE_VECTOR_SHAREDMEM_TEUCHOS_SERIAL_DENSE_HPP_
#define CORE_VECTOR_CONCRETE_VECTOR_SHAREDMEM_TEUCHOS_SERIAL_DENSE_HPP_

#include "../../shared_base/core_container_base.hpp"
#include "../../shared_base/core_container_resizable_base.hpp"
#include "../../shared_base/core_container_subscriptable_base.hpp"
#include "../base/core_vector_sharedmem_base.hpp"

namespace rompp{ namespace core{

template <typename wrapped_type>
class Vector<wrapped_type,
	     core::meta::enable_if_t<
	       core::meta::is_teuchos_serial_dense_vector<wrapped_type>::value
	       >
	     >
  : public ContainerBase< Vector<wrapped_type>, wrapped_type >,
    public VectorSharedMemBase< Vector<wrapped_type> >,
    public ContainerResizableBase<Vector<wrapped_type>, 1>,
    public ContainerSubscriptable1DBase< Vector<wrapped_type>,
     typename details::traits<Vector<wrapped_type>>::scalar_t,
     typename details::traits<Vector<wrapped_type>>::ordinal_t>{

  using this_t = Vector<wrapped_type>;
  using mytraits = typename details::traits<this_t>;
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename  mytraits::ordinal_t;
  using wrap_t = typename mytraits::wrapped_t;

public:
  Vector() = default;

  ~Vector() = default;

  explicit Vector(ord_t insize)
    : data_(insize){}

  explicit Vector(const wrap_t & src)
    : data_(src){}

  Vector(this_t const & other)
    : data_(*other.data()){
  }

  // constructor from any expression, force evaluation
  template <typename T,
  	    core::meta::enable_if_t<
  	      T::is_vector_expression> * = nullptr>
  explicit Vector(const T & expr) : data_(expr.size()){
    //this->resize(expr.size());
    for (ord_t i = 0; i != expr.size(); ++i)
      data_[i] = expr(i);
  }

public:

  // assignment from any expression, force evaluation
  template <typename T,
  	    core::meta::enable_if_t<
  	      T::is_vector_expression> * = nullptr>
  this_t & operator=(const T & expr){
    if(this->size() != expr.size())
      this->resize(expr.size());
    for (ord_t i = 0; i != expr.size(); ++i)
      data_[i] = expr(i);
    return *this;
  }

  // assignment with other vector of same type
  template <typename T,
  	    core::meta::enable_if_t<
  	      std::is_same<T,this_t>::value> * = nullptr>
  this_t & operator=(const T & other){
    data_ = *other.data();
    return *this;
  }

  // assignment with value
  template <typename T,
  	    core::meta::enable_if_t<
  	      std::is_same<T, sc_t>::value> * = nullptr>
  this_t & operator=(const T value){
    for (ord_t i = 0; i != this->size(); ++i)
      data_[i] = value;
    return *this;
  }

  // compound assignment from expression template
  // this += expr
  template <typename T,
  	    core::meta::enable_if_t<
  	      T::is_vector_expression> * = nullptr>
  this_t & operator+=(const T & expr) {
    assert( expr.size() == this->size() );
    for (ord_t i = 0; i != expr.size(); ++i)
      data_[i] += expr(i);
    return *this;
  }

  // compound assignment when type(b) = type(this)
  // this += b
  template <typename T,
  	    core::meta::enable_if_t<
  	      std::is_same<T,this_t>::value> * = nullptr>
  this_t & operator+=(const T & other) {
    assert( other.size() == this->size() );
    this->data_ += *other.data();
    return *this;
  }


  // compound assignment from expression template
  // this -= expr
  template <typename T,
  	    core::meta::enable_if_t<
  	      T::is_vector_expression> * = nullptr>
  this_t & operator-=(const T & expr) {
    assert( expr.size() == this->size() );
    for (ord_t i = 0; i != expr.size(); ++i)
      data_[i] -= expr(i);
    return *this;
  }

  // compound assignment when type(b) = type(this)
  // this -= b
  template <typename T,
  	    core::meta::enable_if_t<
  	      std::is_same<T,this_t>::value> * = nullptr>
  this_t & operator-=(const T & other) {
    assert( other.size() == this->size() );
    this->data_ -= *other.data();
    return *this;
  }

public:
  sc_t & operator [] (ord_t i){
    return data_(i);
  };

  sc_t const & operator [] (ord_t i) const{
    return data_(i);
  };

  sc_t & operator()(ord_t i){
    return data_[i];
  };
  sc_t const & operator()(ord_t i) const{
    return data_[i];
  };

private:

  void matchLayoutWithImpl(const this_t & other){
    this->resize( other.size() );
  }

  wrap_t const * dataImpl() const{
    return &data_;
  }

  wrap_t * dataImpl(){
    return &data_;
  }

  void putScalarImpl(sc_t value) {
    data_.putScalar(value);
  }

  void setZeroImpl() {
    this->putScalar( static_cast<sc_t>(0) );
  }

  bool emptyImpl() const{
    return this->size()==0 ? true : false;
  }

  ord_t sizeImpl() const {
    return data_.length();
  }

  void resizeImpl(ord_t newsize){
    //this will resize and reset to zero
    data_.size(newsize);
  }

private:
  friend ContainerBase< this_t, wrapped_type >;
  friend VectorSharedMemBase< this_t >;
  friend ContainerResizableBase<this_t, 1>;
  friend ContainerSubscriptable1DBase<this_t, sc_t, ord_t>;

private:
  wrap_t data_;

};//end class

}}//end namespace rompp::core
#endif
#endif
