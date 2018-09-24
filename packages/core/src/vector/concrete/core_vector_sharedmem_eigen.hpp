
#ifndef CORE_VECTOR_CONCRETE_VECTOR_SHAREDMEM_EIGEN_HPP_
#define CORE_VECTOR_CONCRETE_VECTOR_SHAREDMEM_EIGEN_HPP_

#include "../../shared_base/core_container_base.hpp"
#include "../../shared_base/core_operators_base.hpp"
#include "../../shared_base/core_container_resizable_base.hpp"
#include "../../shared_base/core_container_nonresizable_base.hpp"

#include "../base/core_vector_sharedmem_base.hpp"
#include "../base/core_vector_math_base.hpp"

namespace rompp{
namespace core{
  
template <typename wrapped_type>
class Vector<wrapped_type,
	     core::meta::enable_if_t<
	       core::meta::is_vector_eigen<wrapped_type>::value
	       >
	     >
  : public ContainerBase< Vector<wrapped_type>, wrapped_type >,
    public VectorSharedMemBase< Vector<wrapped_type> >,
    public VectorMathBase< Vector<wrapped_type> >,
    public ArithmeticOperatorsBase<Vector<wrapped_type>>,
    public CompoundAssignmentOperatorsBase<Vector<wrapped_type>>,
    public std::conditional<
  details::traits<Vector<wrapped_type>>::is_static == true,
  ContainerNonResizableBase<Vector<wrapped_type>, 1>,
  ContainerResizableBase<Vector<wrapped_type>, 1>
  >::type
  
{

  using this_t = Vector<wrapped_type>;
  using mytraits = typename details::traits<this_t>;  
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename  mytraits::ordinal_t;
  using wrap_t = typename mytraits::wrapped_t;

public:

  Vector() = default;

  template <typename T,
	    typename std::enable_if<
	      std::is_same<T,ord_t>::value &&
	      !mytraits::is_static
	      >::type * = nullptr>
  explicit Vector(T insize){
    this->resize(insize);
  }

  explicit Vector(const sc_t * src)
    : data_(src){}

  explicit Vector(const wrap_t & src)
    : data_(src){}

  Vector(this_t const & other)
    : data_(*other.data()){}
  
  ~Vector(){}
  
public:
  sc_t & operator [] (ord_t i){
    //assert(!this->empty());
    return data_(i);
  };

  sc_t const & operator [] (ord_t i) const{
    //assert(!this->empty());
    return data_(i);
  };  

  this_t operator+(const this_t & other) const{
    assert( other.size() == this->size() );
    this_t res(other.size());
    *res.data() = this->data_ + *other.data();
    return res;
  }

  this_t operator-(const this_t & other) const{
    assert( other.size() == this->size() );
    this_t res(other.size());
    *res.data() = this->data_ - *other.data();
    return res;
  }
  
  this_t operator*(const this_t & other) const{
    assert( other.size() == this->size() );
    this_t res(other.size());
    for (decltype(this->size()) i=0; i<this->size(); i++)
      res[i] = this->data_(i) * other[i];
    return res;
  }
  
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

  template <typename T = ord_t,
  	    typename std::enable_if<
  	      !mytraits::is_static, T
  	      >::type * = nullptr>
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

  template <typename T = ord_t,
	    typename std::enable_if<
	      !mytraits::is_static, T
	      >::type * = nullptr>
  void resizeImpl(T val){
    data_.resize(val);
  }

  template<typename op_t, typename T,
	   core::meta::enable_if_t<
	     std::is_same<T,this_t>::value
	     > * = nullptr
	   >
  void inPlaceOpImpl(sc_t a1, sc_t a2, const T & other){
    // this = a1*this op a2*other;
    for (decltype(this->size()) i=0; i<this->size(); i++)
      data_(i) = op_t()( a1*data_(i), a2*other[i] );
  }

  template<typename op_t, typename T,
  	   core::meta::enable_if_t<
  	     std::is_same<T,this_t>::value
  	     > * = nullptr
  	   >
  void inPlaceOpImpl(sc_t a1, const T & x1,
  		     sc_t a2, const T & x2){
    // this = a1*x1 op a2*x2;
    assert(this->size() == x1.size());
    assert(this->size() == x2.size());
    for (ord_t i=0; i<this->size(); i++)
      data_(i) = op_t()( a1*x1[i], a2*x2[i] );
  }

  template<typename T,
  	   typename op1_t, typename op2_t, typename op3_t,
  	   core::meta::enable_if_t<
  	     std::is_same<T,this_t>::value &&
  	     std::is_same<op1_t, std::plus<sc_t>>::value &&
  	     std::is_same<op2_t,op1_t>::value &&
  	     std::is_same<op3_t,op1_t>::value
  	     > * = nullptr
  	   >
  void inPlaceOpImpl(sc_t a1, const T & x1,
  		     sc_t a2, const T & x2,
  		     sc_t a3, const T & x3,
  		     sc_t a4, const T & x4){
    // this = (a1*x1) + (a2*x2) + (a3*x3) + (a4*x4)
    assert(this->size() == x1.size());
    assert(this->size() == x2.size());
    assert(this->size() == x3.size());
    assert(this->size() == x4.size());

    for (ord_t i=0; i<this->size(); i++)
      data_(i) = op1_t()( op1_t()( op1_t()(a1*x1[i],a2*x2[i]), a3*x3[i] ), a4*x4[i] );
  }
  
  template<typename op0_t, typename T,
	   typename op1_t, typename op2_t, typename op3_t,
  	   core::meta::enable_if_t<
  	     std::is_same<T,this_t>::value &&
  	     std::is_same<op0_t, std::plus<sc_t>>::value &&
  	     std::is_same<op1_t, op0_t>::value &&
  	     std::is_same<op2_t, op1_t>::value &&
  	     std::is_same<op3_t, op1_t>::value
  	     > * = nullptr
  	   >
  void inPlaceOpImpl(sc_t a0, sc_t a1, const T & x1,
  		     sc_t a2, const T & x2,
  		     sc_t a3, const T & x3,
  		     sc_t a4, const T & x4){
    // this = a0 * this + (a1*x1) + (a2*x2) + (a3*x3) + (a4*x4)
    assert(this->size() == x1.size());
    assert(this->size() == x2.size());
    assert(this->size() == x3.size());
    assert(this->size() == x4.size());
    for (ord_t i=0; i<this->size(); i++)
      data_(i) = a0*data_(i) + a1*x1[i] + a2*x2[i]
	+ a3*x3[i] + a4*x4[i];
  }

  
  void scaleImpl(sc_t & factor){
    // this = factor * this;
    for (decltype(this->size()) i=0; i<this->size(); i++)
      data_(i) *= factor;
  }

  void norm1Impl(sc_t & result) const {
    result = static_cast<sc_t>(0);
    for (decltype(this->size()) i=0; i<this->size(); i++)
      result += std::abs(data_(i));
  }

  void norm2Impl(sc_t & res) const {
    res = static_cast<sc_t>(0);
    for (decltype(this->size()) i=0; i<this->size(); i++)
      res += data_(i)*data_(i);
    res = std::sqrt(res);
  }

  void normInfImpl(sc_t & res) const {
    res = std::abs(data_(0));
    for (decltype(this->size()) i=1; i<this->size(); i++){
      sc_t currVal = std::abs(data(i));
      if(currVal>res)
	res = currVal;
    }
  }
  void minValueImpl(sc_t & result) const {
    result = data_.minCoeff();
  }
  void maxValueImpl(sc_t & result) const {
    result = data_.maxCoeff();
  }
  
private:
  friend ContainerBase< this_t, wrapped_type >;
  friend VectorSharedMemBase< this_t >;
  friend VectorMathBase< this_t >;
  friend ArithmeticOperatorsBase< this_t >;
  friend CompoundAssignmentOperatorsBase< this_t >;  
  friend typename std::conditional<
    details::traits<this_t>::is_static == true,
    ContainerNonResizableBase<this_t, 1>,
    ContainerResizableBase<this_t, 1>
    >::type;

private:
  wrap_t data_;
 
};//end class    
}//end namespace core
}//end namespace rompp
#endif
