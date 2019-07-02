
#ifdef HAVE_BLAZE
#ifndef CONTAINERS_VECTOR_CONCRETE_VECTOR_SHAREDMEM_BLAZE_DYNAMIC_HPP_
#define CONTAINERS_VECTOR_CONCRETE_VECTOR_SHAREDMEM_BLAZE_DYNAMIC_HPP_

#include "../../shared_base/containers_container_base.hpp"
#include "../../shared_base/containers_container_resizable_base.hpp"
#include "../../shared_base/containers_container_subscriptable_base.hpp"
#include "../base/containers_vector_sharedmem_base.hpp"

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Vector<wrapped_type,
    ::pressio::mpl::enable_if_t<
      containers::meta::is_dynamic_vector_blaze<wrapped_type>::value
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
  	    ::pressio::mpl::enable_if_t<
  	      T::is_vector_expression> * = nullptr>
  explicit Vector(const T & expr){
    this->resize(expr.size());
    for (ord_t i = 0; i != expr.size(); ++i)
      data_[i] = expr(i);
  }

  // assignment from any expression, force evaluation
  template <typename T,
  	    ::pressio::mpl::enable_if_t<
  	      T::is_vector_expression> * = nullptr>
  this_t & operator=(const T & expr){
    if(this->size() != expr.size())
      this->resize(expr.size());
    for (ord_t i = 0; i != expr.size(); ++i)
      data_[i] = expr(i);
    return *this;
  }

  // assignment
  template <typename T,
  	    ::pressio::mpl::enable_if_t<
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

  sc_t & operator()(ord_t i){
    return data_[i];
  };
  sc_t const & operator()(ord_t i) const{
    return data_[i];
  };


  // compound assignment from expression template
  // this += expr
  template <typename T,
  	    ::pressio::mpl::enable_if_t<
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
  	    ::pressio::mpl::enable_if_t<
  	      std::is_same<T,this_t>::value> * = nullptr>
  this_t & operator+=(const T & other) {
    assert( other.size() == this->size() );
    this->data_ += *other.data();
    return *this;
  }


  // compound assignment from expression template
  // this -= expr
  template <typename T,
  	    ::pressio::mpl::enable_if_t<
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
  	    ::pressio::mpl::enable_if_t<
  	      std::is_same<T,this_t>::value> * = nullptr>
  this_t & operator-=(const T & other) {
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

  bool emptyImpl() const{
    return this->size()==0 ? true : false;
  }

  ord_t sizeImpl() const {
    return data_.size();
  }

  void resizeImpl(ord_t val){
    data_.resize(val);
  }


private:
  friend ContainerBase< this_t, wrapped_type >;
  friend VectorSharedMemBase< this_t >;
  friend ContainerResizableBase<this_t, 1>;
  friend ContainerSubscriptable1DBase<this_t, sc_t, ord_t>;

private:
  wrap_t data_ = {};

};//end class
}}//end namespace pressio::containers

#endif /* CONTAINERS_VECTOR_CONCRETE_VECTOR_SHAREDMEM_BLAZE_DYNAMIC_HPP_ */
#endif //HAVE_BLAZE
