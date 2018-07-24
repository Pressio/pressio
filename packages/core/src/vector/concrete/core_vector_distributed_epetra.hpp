
#ifndef CORE_VECTOR_CONCRETE_VECTOR_DISTRIBUTED_EPETRA_HPP_
#define CORE_VECTOR_CONCRETE_VECTOR_DISTRIBUTED_EPETRA_HPP_

#include "../../shared_base/core_container_base.hpp"
#include "../base/core_vector_distributed_base.hpp"
#include "../base/core_vector_distributed_trilinos.hpp"
#include "../base/core_vector_math_base.hpp"
#include "../../shared_base/core_container_distributed_base.hpp"
#include "../../shared_base/core_operators_base.hpp"

namespace core{
  
template <typename wrapped_type>
class Vector<wrapped_type,
	     typename
	     std::enable_if<
	       meta::is_vector_epetra<
		 wrapped_type>::value
	       >::type
	     >
  : public ContainerBase< Vector<wrapped_type>, wrapped_type >,
    public VectorDistributedBase< Vector<wrapped_type> >,
    public VectorDistributedTrilinos< Vector<wrapped_type> >,
    public VectorMathBase< Vector<wrapped_type> >, 
    public ContainerDistributedBase< Vector<wrapped_type>, 
              typename details::traits<Vector<wrapped_type>>::communicator_t >, 
    public Subscripting1DOperatorsBase< Vector<wrapped_type>, 
              typename details::traits<Vector<wrapped_type>>::scalar_t,
              typename details::traits<Vector<wrapped_type>>::local_ordinal_t>
{
private:
  using derived_t = Vector<wrapped_type>;
  using sc_t = typename details::traits<derived_t>::scalar_t;
  using LO_t = typename details::traits<derived_t>::local_ordinal_t;
  using GO_t = typename details::traits<derived_t>::global_ordinal_t;
  using der_t = typename details::traits<derived_t>::derived_t;
  using wrap_t = typename details::traits<derived_t>::wrapped_t;
  using map_t = typename details::traits<derived_t>::data_map_t;
  using mpicomm_t = typename details::traits<derived_t>::communicator_t;

public:
  Vector() = delete;
  explicit Vector(const map_t & mapobj) : data_(mapobj){}
  explicit Vector(const wrap_t & vecobj) : data_(vecobj){}
  ~Vector() = default;

public:
  sc_t & operator [] (LO_t i){
    assert(i < this->localSize());
    return data_[i];
  };
  sc_t const & operator [] (LO_t i) const{
    assert(i < this->localSize());
    return data_[i];
  };  

  der_t operator+(const der_t & other) const{
    der_t res( other.getDataMap() );
    for (LO_t i=0; i<this->localSize(); i++)
      res[i] = (*this)[i] + other[i];
    return res;
  }

  der_t operator-(const der_t & other) const{
    der_t res( other.getDataMap() );
    for (LO_t i=0; i<this->localSize(); i++)
      res[i] = (*this)[i] - other[i];
    return res;
  }

  der_t operator*(const der_t & other) const{
    der_t res( other.getDataMap() );
    for (LO_t i=0; i<this->localSize(); i++)
      res[i] = (*this)[i] * other[i];
    return res;
  }

  der_t & operator+=(const der_t & other) {
    this->data_.Update(1.0, *other.data(), 1.0 );
    return *this;
  }

  der_t & operator-=(const der_t & other) {
    this->data_.Update(-1.0, *other.data(), 1.0 );
    return *this;
  }
    
private:

  mpicomm_t const & commCRef() const{
    return data_.Comm();
  }

  wrap_t const * dataImpl() const{
    return &data_;
  }

  wrap_t * dataImpl(){
    return &data_;
  }

  void putScalarImpl(sc_t value) {
    data_.PutScalar(value);
  }

  void setZeroImpl() {
    this->putScalarImpl( static_cast<sc_t>(0) );
  }

  bool emptyImpl() const{
    return this->globalSize()==0 ? true : false;
  }
  
  //----------------
  //from distributed base
  //----------------
  GO_t globalSizeImpl() const {
    return data_.GlobalLength();
  }

  LO_t localSizeImpl() const {
    return data_.MyLength();
  }

  void replaceGlobalValuesImpl(GO_t numentries,
			       const GO_t * indices,
			       const sc_t * values){
    data_.ReplaceGlobalValues(numentries, values, indices);
  }

  //-------------------
  //from trilinos base
  //-------------------
  map_t const & getDataMapImpl() const{
    return data_.Map();
  }

  void replaceDataMapImpl(const map_t & mapObj){
    data_.ReplaceMap(mapObj);
  }
    
  //----------------
  //from math base
  //----------------
  template<typename op_t>
  void inPlaceOpImpl(op_t op, sc_t a1, sc_t a2, const der_t & other){
    // this = a1*this op a2*other;
    for (LO_t i=0; i<this->localSize(); i++)
      data_[i] = op()( a1*data_[i], a2*other[i] );
  }

  void scaleImpl(sc_t & factor){
    data_.Scale(factor);
  }

  void norm1Impl(sc_t & result) const {
    data_.Norm1(&result);
  }

  void norm2Impl(sc_t & result) const {
    data_.Norm2(&result);
  }

  void normInfImpl(sc_t & result) const {
    data_.NormInf(&result);
  }

  void minValueImpl(sc_t & result) const {
    data_.MinValue(&result);
  }

  void maxValueImpl(sc_t & result) const {
    data_.MaxValue(&result);
  }
  
private:
  friend ContainerBase< derived_t, wrapped_type >;
  friend VectorDistributedBase< derived_t >;
  friend VectorDistributedTrilinos< derived_t >;
  friend VectorMathBase< derived_t >;
  friend ContainerDistributedBase< derived_t, mpicomm_t >;
  friend Subscripting1DOperatorsBase< derived_t, sc_t, LO_t>;

private:
  wrap_t data_;

};//end class
}//end namespace core
#endif
