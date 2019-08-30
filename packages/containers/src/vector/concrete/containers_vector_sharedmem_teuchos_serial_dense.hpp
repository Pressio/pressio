
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_VECTOR_CONCRETE_VECTOR_SHAREDMEM_TEUCHOS_SERIAL_DENSE_HPP_
#define CONTAINERS_VECTOR_CONCRETE_VECTOR_SHAREDMEM_TEUCHOS_SERIAL_DENSE_HPP_

#include "../../shared_base/containers_container_base.hpp"
#include "../../shared_base/containers_container_resizable_base.hpp"
#include "../../shared_base/containers_container_subscriptable_base.hpp"
#include "../base/containers_vector_sharedmem_base.hpp"

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Vector<wrapped_type,
	     ::pressio::mpl::enable_if_t<
	       containers::meta::is_dense_vector_teuchos<wrapped_type>::value
	       >
	     >
  : public ContainerBase< Vector<wrapped_type>, wrapped_type >,
    public VectorSharedMemBase< Vector<wrapped_type> >,
    public ContainerResizableBase<Vector<wrapped_type>, 1>,
    public ContainerSubscriptable1DBase< Vector<wrapped_type>,
     typename details::traits<Vector<wrapped_type>>::scalar_t,
     typename details::traits<Vector<wrapped_type>>::ordinal_t>,
    public ContainerPrintable1DBase<
	       Vector<wrapped_type>,
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
  	    ::pressio::mpl::enable_if_t<
  	      T::is_vector_expression> * = nullptr>
  explicit Vector(const T & expr) : data_(expr.size()){
    //this->resize(expr.size());
    for (ord_t i = 0; i != expr.size(); ++i)
      data_[i] = expr(i);
  }

public:

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

  // assignment with other vector of same type
  this_t & operator=(const this_t & other){
    data_ = *other.data();
    return *this;
  }

  // assignment with value
  template <typename T,
  	    ::pressio::mpl::enable_if_t<
  	      std::is_same<T, sc_t>::value> * = nullptr>
  this_t & operator=(const T value){
    for (ord_t i = 0; i != this->size(); ++i)
      data_[i] = value;
    return *this;
  }

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
  this_t & operator+=(const this_t & other) {
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
  this_t & operator-=(const this_t & other) {
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

  template <typename stream_t>
  void printImpl(stream_t & os, char c, ord_t nIn) const{
    assert(nIn <= this->size());
    auto nToPrint = (nIn==-1) ? this->size() : nIn;
    ::pressio::utils::impl::setStreamPrecision<stream_t, sc_t>(os);
    if (c=='d') this->printVertically(os, nToPrint);
    if (c=='f') this->printFlatten(os, nToPrint);
  }

  template <typename stream_t>
  void printVertically(stream_t & os, ord_t nToPrint) const{
    for (auto i=0; i<nToPrint; i++)
      os << data_[i] << "\n";
    os << std::endl;
  }

  template <typename stream_t>
  void printFlatten(stream_t & os, ord_t nToPrint) const{
    for (auto i=0; i<nToPrint; i++)
      os << data_[i] << " ";
    os << std::endl;
  }

  void matchLayoutWithImpl(const this_t & other){
    this->resize( other.size() );
  }

  void scaleImpl(sc_t value) {
    data_.scale(value);
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
  friend ContainerPrintable1DBase<this_t, ord_t>;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers
#endif
#endif
