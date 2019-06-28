
#ifdef HAVE_TRILINOS
#ifndef ALGEBRA_VECTOR_CONCRETE_VECTOR_DISTRIBUTED_EPETRA_HPP_
#define ALGEBRA_VECTOR_CONCRETE_VECTOR_DISTRIBUTED_EPETRA_HPP_

#include "../../shared_base/algebra_container_base.hpp"
#include "../../shared_base/algebra_container_distributed_mpi_base.hpp"
#include "../../shared_base/algebra_container_distributed_trilinos_base.hpp"
#include "../../shared_base/algebra_container_subscriptable_base.hpp"
#include "../../shared_base/algebra_container_resizable_base.hpp"
#include "../base/algebra_vector_distributed_base.hpp"

namespace rompp{ namespace algebra{

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
    public ContainerDistributedMpiBase< Vector<wrapped_type>,
     typename details::traits<Vector<wrapped_type>>::communicator_t >,
    public ContainerDistributedTrilinosBase< Vector<wrapped_type>,
     typename details::traits<Vector<wrapped_type>>::data_map_t >,
    public ContainerResizableBase< Vector<wrapped_type>, 1>,
    public ContainerSubscriptable1DBase< Vector<wrapped_type>,
     typename details::traits<Vector<wrapped_type>>::scalar_t,
     typename details::traits<Vector<wrapped_type>>::local_ordinal_t>{

  using this_t = Vector<wrapped_type>;
  using sc_t = typename details::traits<this_t>::scalar_t;
  using LO_t = typename details::traits<this_t>::local_ordinal_t;
  using GO_t = typename details::traits<this_t>::global_ordinal_t;
  using der_t = this_t;
  using wrap_t = typename details::traits<this_t>::wrapped_t;
  using map_t = typename details::traits<this_t>::data_map_t;
  using mpicomm_t = typename details::traits<this_t>::communicator_t;

public:
  Vector() = delete;
  ~Vector() = default;

  explicit Vector(const map_t & mapobj)
    : data_(mapobj){}

  explicit Vector(const wrap_t & vecobj)
    : data_(vecobj){}

  Vector(this_t const & other)
    : data_(*other.data()){}

  // assignment from any expression, force evaluation
  template <typename T,
	    ::rompp::mpl::enable_if_t<
	      T::is_vector_expression> * = nullptr>
  this_t & operator=(const T & expr){
    assert(this->localSize() == expr.localSize());
    for (LO_t i = 0; i != expr.localSize(); ++i)
      data_[i] = expr(i);
    return *this;
  }

public:
  sc_t & operator [] (LO_t i){
    assert(i < this->localSize());
    return data_[i];
  };
  sc_t const & operator [] (LO_t i) const{
    assert(i < this->localSize());
    return data_[i];
  };

  sc_t & operator()(LO_t i){
    assert(i < this->localSize());
    return data_[i];
  };
  sc_t const & operator()(LO_t i) const{
    assert(i < this->localSize());
    return data_[i];
  };


  // compound assignment from expression template
  // this += expr
  template <typename T,
  	    ::rompp::mpl::enable_if_t<
  	      T::is_vector_expression> * = nullptr>
  this_t & operator+=(const T & expr) {
    assert(this->localSize() == expr.localSize());
    for (LO_t i = 0; i != expr.localSize(); ++i)
      data_[i] += expr(i);
    return *this;
  }

  // compound assignment when type(b) = type(this)
  // this += b
  template <typename T,
  	    ::rompp::mpl::enable_if_t<
  	      std::is_same<T,this_t>::value> * = nullptr>
  this_t & operator+=(const T & other) {
    this->data_.Update(1.0, *other.data(), 1.0 );
    return *this;
  }


  // compound assignment from expression template
  // this -= expr
  template <typename T,
  	    ::rompp::mpl::enable_if_t<
  	      T::is_vector_expression> * = nullptr>
  this_t & operator-=(const T & expr) {
    assert(this->localSize() == expr.localSize());
    for (LO_t i = 0; i != expr.localSize(); ++i)
      data_[i] -= expr(i);
    return *this;
  }

  // compound assignment when type(b) = type(this)
  // this -= b
  template <typename T,
  	    ::rompp::mpl::enable_if_t<
  	      std::is_same<T,this_t>::value> * = nullptr>
  this_t & operator-=(const T & other) {
    this->data_.Update(-1.0, *other.data(), 1.0 );
    return *this;
  }


  void print(std::string tag) const{
    ::rompp::utils::io::print_stdout(tag);
    data_.Print(std::cout);
  }

private:

  void matchLayoutWithImpl(const der_t & other){
    data_.ReplaceMap( other.getDataMap() );
  }

  mpicomm_t const & commCRefImpl() const{
    return data_.Comm();
  }

  wrap_t const * dataImpl() const{
    return &data_;
  }

  wrap_t * dataImpl(){
    return &data_;
  }

  bool isDistributedGloballyImpl() const{
    return data_.DistributedGlobal();
  }

  void putScalarImpl(sc_t value) {
    data_.PutScalar(value);
  }

  void setZeroImpl(){
    data_.PutScalar(static_cast<sc_t>(0));
  }

  bool emptyImpl() const{
    return this->globalSize()==0 ? true : false;
  }

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

  map_t const & getDataMapImpl() const{
    return data_.Map();
  }

  void replaceDataMapImpl(const map_t & mapObj){
    data_.ReplaceMap(mapObj);
  }

private:
  friend ContainerBase< this_t, wrapped_type >;
  friend VectorDistributedBase< this_t >;
  friend ContainerDistributedMpiBase< this_t, mpicomm_t >;
  friend ContainerDistributedTrilinosBase< this_t, map_t >;
  friend ContainerResizableBase< this_t, 1>;
  friend ContainerSubscriptable1DBase< this_t, sc_t, LO_t>;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace rompp::algebra
#endif
#endif






  // template<typename op_t, typename T,
  // 	   ::rompp::mpl::enable_if_t<
  // 	     std::is_same<T,this_t>::value
  // 	     > * = nullptr
  // 	   >
  // void inPlaceOpImpl(sc_t a1, sc_t a2, const T & other){
  //   // this = a1*this op a2*other;
  //   for (LO_t i=0; i<this->localSize(); i++)
  //     data_[i] = op_t()( a1*data_[i], a2*other[i] );
  // }


  // template<typename op_t, typename T,
  // 	   ::rompp::mpl::enable_if_t<
  // 	     std::is_same<T,this_t>::value
  // 	     > * = nullptr
  // 	   >
  // void inPlaceOpImpl(sc_t a1, const T & x1,
  // 		     sc_t a2, const T & x2){
  //   // this = a1*x1 op a2*x2;
  //   assert(this->globalSizeImpl() == x1.globalSizeImpl());
  //   assert(this->globalSizeImpl() == x2.globalSizeImpl());

  //   for (LO_t i=0; i<this->localSize(); i++)
  //     data_[i] = op_t()( a1*x1[i], a2*x2[i] );
  // }

  // template<typename op0_t, typename T,
  // 	   typename op1_t, typename op2_t, typename op3_t,
  // 	   ::rompp::mpl::enable_if_t<
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
  //   assert(this->globalSizeImpl() == x1.globalSizeImpl());
  //   assert(this->globalSizeImpl() == x2.globalSizeImpl());
  //   assert(this->globalSizeImpl() == x3.globalSizeImpl());
  //   assert(this->globalSizeImpl() == x4.globalSizeImpl());

  //   for (LO_t i=0; i<this->localSize(); i++)
  //     data_[i] = a0*data_[i] + a1*x1[i] + a2*x2[i]
  // 	+ a3*x3[i] + a4*x4[i];
  // }

  // void minValueImpl(sc_t & result) const {
  //   data_.MinValue(&result);
  // }

  // void maxValueImpl(sc_t & result) const {
  //   data_.MaxValue(&result);
  // }
