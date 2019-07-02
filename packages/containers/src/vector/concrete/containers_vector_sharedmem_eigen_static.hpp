
#ifndef CONTAINERS_VECTOR_CONCRETE_VECTOR_SHAREDMEM_EIGEN_STATIC_HPP_
#define CONTAINERS_VECTOR_CONCRETE_VECTOR_SHAREDMEM_EIGEN_STATIC_HPP_

#include "../../shared_base/containers_container_base.hpp"
#include "../../shared_base/containers_container_nonresizable_base.hpp"
#include "../../shared_base/containers_container_subscriptable_base.hpp"
#include "../../shared_base/containers_container_printable_base.hpp"
#include "../base/containers_vector_sharedmem_base.hpp"

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Vector<wrapped_type,
	     ::pressio::mpl::enable_if_t<
	       containers::meta::is_static_vector_eigen<wrapped_type>::value
	       >
	     >
  : public ContainerBase< Vector<wrapped_type>, wrapped_type >,
    public VectorSharedMemBase< Vector<wrapped_type> >,
    public ContainerNonResizableBase<Vector<wrapped_type>, 1>,
    public ContainerSubscriptable1DBase<
	      Vector<wrapped_type>,
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

  explicit Vector(const sc_t * src)
    : data_(src){}

  explicit Vector(const wrap_t & src)
    : data_(src){}

  Vector(this_t const & other)
    : data_(*other.data()){}

  // assignment from any expression, force evaluation
  template <typename T,
	    ::pressio::mpl::enable_if_t<
	      T::is_vector_expression> * = nullptr>
  this_t & operator=(const T & expr){
    assert(this->size() == expr.size());
    for (decltype(expr.size()) i = 0; i != expr.size(); ++i)
      data_[i] = expr(i);
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
    return data_(i);
  };
  sc_t const & operator()(ord_t i) const{
    return data_(i);
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

  template <typename stream_t>
  void printImpl(stream_t & os, char c, ord_t nIn) const{
    assert(nIn <= this->size());
    auto nToPrint = (nIn==-1) ? this->size() : nIn;
    ::pressio::utils::impl::setStreamPrecision<stream_t, sc_t>(os);

    if (c=='d')
      this->printVertically(os, nToPrint);
    if (c=='f')
      this->printFlatten(os, nToPrint);
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

  wrap_t const * dataImpl() const{
    return &data_;
  }

  void scaleImpl(sc_t value) {
    data_ *= value;
  }

  wrap_t * dataImpl(){
    return &data_;
  }

  void putScalarImpl(sc_t value) {
    data_.setConstant(value);
  }

  void setZeroImpl() {
    this->putScalarImpl( static_cast<sc_t>(0) );
  }

  bool emptyImpl() const{
    return this->size()==0 ? true : false;
  }

  ord_t sizeImpl() const {
    return (data_.rows()==1) ? data_.cols() : data_.rows();
  }

private:
  friend ContainerBase< this_t, wrapped_type >;
  friend VectorSharedMemBase< this_t >;
  friend ContainerNonResizableBase<this_t, 1>;
  friend ContainerSubscriptable1DBase<this_t, sc_t, ord_t>;
  friend ContainerPrintable1DBase<this_t, ord_t>;

private:
  wrap_t data_ = {};

};//end class
}}//end namespace pressio::containers
#endif







  // template<typename op_t, typename T,
  // 	   ::pressio::mpl::enable_if_t<
  // 	     std::is_same<T,this_t>::value
  // 	     > * = nullptr
  // 	   >
  // void inPlaceOpImpl(sc_t a1, sc_t a2, const T & other){
  //   // this = a1*this op a2*other;
  //   for (decltype(this->size()) i=0; i<this->size(); i++)
  //     data_(i) = op_t()( a1*data_(i), a2*other[i] );
  // }

  // template<typename op_t, typename T,
  // 	   ::pressio::mpl::enable_if_t<
  // 	     std::is_same<T,this_t>::value
  // 	     > * = nullptr
  // 	   >
  // void inPlaceOpImpl(sc_t a1, const T & x1,
  // 		     sc_t a2, const T & x2){
  //   // this = a1*x1 op a2*x2;
  //   assert(this->size() == x1.size());
  //   assert(this->size() == x2.size());
  //   for (ord_t i=0; i<this->size(); i++)
  //     data_(i) = op_t()( a1*x1[i], a2*x2[i] );
  // }

  // template<typename T,
  // 	   typename op1_t, typename op2_t, typename op3_t,
  // 	   ::pressio::mpl::enable_if_t<
  // 	     std::is_same<T,this_t>::value &&
  // 	     std::is_same<op1_t, std::plus<sc_t>>::value &&
  // 	     std::is_same<op2_t,op1_t>::value &&
  // 	     std::is_same<op3_t,op1_t>::value
  // 	     > * = nullptr
  // 	   >
  // void inPlaceOpImpl(sc_t a1, const T & x1,
  // 		     sc_t a2, const T & x2,
  // 		     sc_t a3, const T & x3,
  // 		     sc_t a4, const T & x4){
  //   // this = (a1*x1) + (a2*x2) + (a3*x3) + (a4*x4)
  //   assert(this->size() == x1.size());
  //   assert(this->size() == x2.size());
  //   assert(this->size() == x3.size());
  //   assert(this->size() == x4.size());

  //   for (ord_t i=0; i<this->size(); i++)
  //     data_(i) = op1_t()( op1_t()( op1_t()(a1*x1[i],a2*x2[i]), a3*x3[i] ), a4*x4[i] );
  // }

  // template<typename op0_t, typename T,
  // 	   typename op1_t, typename op2_t, typename op3_t,
  // 	   ::pressio::mpl::enable_if_t<
  // 	     std::is_same<T,this_t>::value &&
  // 	     std::is_same<op0_t, std::plus<sc_t>>::value &&
  // 	     std::is_same<op1_t, op0_t>::value &&
  // 	     std::is_same<op2_t, op1_t>::value &&
  // 	     std::is_same<op3_t, op1_t>::value
  // 	     > * = nullptr
  // 	   >
  // void inPlaceOpImpl(sc_t a0, sc_t a1, const T & x1,
  // 		     sc_t a2, const T & x2,
  // 		     sc_t a3, const T & x3,
  // 		     sc_t a4, const T & x4){
  //   // this = a0 * this + (a1*x1) + (a2*x2) + (a3*x3) + (a4*x4)
  //   assert(this->size() == x1.size());
  //   assert(this->size() == x2.size());
  //   assert(this->size() == x3.size());
  //   assert(this->size() == x4.size());
  //   for (ord_t i=0; i<this->size(); i++)
  //     data_(i) = a0*data_(i) + a1*x1[i] + a2*x2[i]
  // 	+ a3*x3[i] + a4*x4[i];
  // }
