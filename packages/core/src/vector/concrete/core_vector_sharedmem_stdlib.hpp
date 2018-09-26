
#ifndef CORE_VECTOR_CONCRETE_VECTOR_STDLIB_HPP_
#define CORE_VECTOR_CONCRETE_VECTOR_STDLIB_HPP_

#include "../../shared_base/core_container_base.hpp"
#include "../../shared_base/core_operators_base.hpp"
#include "../base/core_vector_sharedmem_base.hpp"
#include "../base/core_vector_math_base.hpp"


namespace rompp{
namespace core{

template <typename wrapped_type>
class Vector<wrapped_type,
	     typename std::enable_if<
	       core::meta::is_vector_stdlib<wrapped_type>::value
	       >::type
	     >
  : public ContainerBase< Vector<wrapped_type>, wrapped_type >,
    public VectorSharedMemBase< Vector<wrapped_type> >,
    public VectorMathBase< Vector<wrapped_type> >,
    public CompoundAssignmentOperatorsBase< Vector<wrapped_type> >,
    public Subscripting1DOperatorsBase< Vector<wrapped_type>, 
              typename details::traits<Vector<wrapped_type>>::scalar_t,
              typename details::traits<Vector<wrapped_type>>::ordinal_t>
{
private:
  using this_t = Vector<wrapped_type>;
  using sc_t = typename details::traits<this_t>::scalar_t;
  using der_t = typename details::traits<this_t>::derived_t;
  using wrap_t = typename details::traits<this_t>::wrapped_t;
  using ord_t = typename details::traits<this_t>::ordinal_t;

public:
  Vector() = default;
  ~Vector(){}
  
  explicit Vector(ord_t insize,
                  sc_t value = static_cast<sc_t>(0) ){
    this->resize(insize, value);
  }

  explicit Vector(const std::vector<sc_t> & src)
    : data_(src){}

  // constructor from any expression, force evaluation
  template <typename T,
	    core::meta::enable_if_t<
	      T::is_vector_expression> * = nullptr>
  Vector(const T & expr){
    this->resize(expr.size());
    for (size_t i = 0; i != expr.size(); ++i)
      data_[i] = expr[i];
  }

  // assignment from any expression, force evaluation
  template <typename T,
	    core::meta::enable_if_t<
	      T::is_vector_expression> * = nullptr>
  this_t & operator=(const T & expr){
    if(this->size() != expr.size())
      this->resize(expr.size());
    for (size_t i = 0; i != expr.size(); ++i)
      data_[i] = expr[i];
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
    std::transform(this->data_.begin(), this->data_->end(),
		   other.data()->begin(), this->data_.begin(),
		   std::plus<sc_t>());
    return *this;
  }
  
  this_t & operator-=(const this_t & other) {
    std::transform(this->data_.begin(), this->data_->end(),
		   other.data()->begin(), this->data_.begin(),
		   std::minus<sc_t>());
    return *this;
  }


private:
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
  //from sharedMem base
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
  friend ContainerBase< this_t, wrapped_type >;
  friend VectorSharedMemBase< this_t >;
  friend VectorMathBase< this_t >;
  friend CompoundAssignmentOperatorsBase< this_t >;  
  friend Subscripting1DOperatorsBase< this_t, sc_t, ord_t>;

private:
  std::vector<sc_t> data_;
  
};//end class  
}//end namespace core  
}//end namespace rompp
#endif
