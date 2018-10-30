
#ifdef HAVE_TRILINOS
#ifndef CORE_VECTOR_CONCRETE_VECTOR_DISTRIBUTED_TPETRA_HPP_
#define CORE_VECTOR_CONCRETE_VECTOR_DISTRIBUTED_TPETRA_HPP_

#include "../../shared_base/core_container_base.hpp"
#include "../../shared_base/core_container_distributed_mpi_base.hpp"
#include "../../shared_base/core_container_distributed_trilinos_base.hpp"
#include "../../shared_base/core_container_subscriptable_base.hpp"
#include "../../shared_base/core_container_resizable_base.hpp"
#include "../base/core_vector_distributed_base.hpp"

namespace rompp{ namespace core{
  
template <typename wrapped_type>
class Vector<wrapped_type,
	     typename
	     std::enable_if<
	       meta::is_vector_tpetra<
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

  explicit Vector(const wrap_t & vecobj)
    : data_(vecobj){
    myLocDataView_ = data_.getDataNonConst();
    myLocDataCView_= data_.getData();
  }

  Vector(this_t const & other)
    : data_(*other.data()){
    myLocDataView_ = data_.getDataNonConst();
    myLocDataCView_ = data_.getData();
  }
  
  ~Vector() = default;

public:
  sc_t & operator [] (LO_t i){
    assert(!myLocDataView_.is_null());
    assert(i < this->localSize());
    return myLocDataView_[i];
  };
  sc_t const & operator [] (LO_t i) const{
    assert(!myLocDataCView_.is_null());
    assert(i < this->localSize());
    return myLocDataCView_[i];
  };  

  sc_t & operator()(LO_t i){
    assert(!myLocDataView_.is_null());
    assert(i < this->localSize());
    return myLocDataView_[i];
  };
  sc_t const & operator()(LO_t i) const{
    assert(!myLocDataCView_.is_null());
    assert(i < this->localSize());
    return myLocDataCView_[i];
  };  

  // // compound assignment when type(b) = type(this)
  // This gives a linking issue for daxpy in kokkos:blas
  // // this += b 
  // template <typename T,
  // 	    core::meta::enable_if_t<
  // 	      std::is_same<T,this_t>::value> * = nullptr>
  // this_t & operator+=(const T & other) {
  //   this->data_.update(1.0, *other.data(), 1.0 );
  //   return *this;
  // }

//   // compound assignment when type(b) = type(this)
//   // this -= b 
//   template <typename T,
//   	    core::meta::enable_if_t<
//   	      std::is_same<T,this_t>::value> * = nullptr>
//   this_t & operator-=(const T & other) {
//     this->data_.update(-1.0, *other.data(), 1.0 );
//     return *this;
//   }
      
private:

//   void matchLayoutWithImpl(const der_t & other){
//     data_.ReplaceMap( other.getDataMap() );
//   }
  
//   mpicomm_t const & commCRefImpl() const{
//     return data_.Comm();
//   }

  wrap_t const * dataImpl() const{
    return &data_;
  }

  wrap_t * dataImpl(){
    return &data_;
  }

  bool isDistributedGloballyImpl() const{
    return data_.isDistributed();
  }
  
  void putScalarImpl(sc_t value) {
    data_.putScalar(value);
  }

  void setZeroImpl(){
    data_.putScalar(static_cast<sc_t>(0));
  }

  bool emptyImpl() const{
    return this->globalSize()==0 ? true : false;
  }
  
  GO_t globalSizeImpl() const {
    return data_.getGlobalLength();
  }

  LO_t localSizeImpl() const {
    return data_.getLocalLength();
  }

//   void replaceGlobalValuesImpl(GO_t numentries,
// 			       const GO_t * indices,
// 			       const sc_t * values){
//     data_.ReplaceGlobalValues(numentries, values, indices);
//   }

//   map_t const & getDataMapImpl() const{
//     return data_.Map();
//   }

//   void replaceDataMapImpl(const map_t & mapObj){
//     data_.ReplaceMap(mapObj);
//   }
  
private:
  friend ContainerBase< this_t, wrapped_type >;
  friend VectorDistributedBase< this_t >;
  friend ContainerDistributedMpiBase< this_t, mpicomm_t >;
  friend ContainerDistributedTrilinosBase< this_t, map_t >;
  friend ContainerResizableBase< this_t, 1>;
  friend ContainerSubscriptable1DBase< this_t, sc_t, LO_t>;

private:
  wrap_t data_;

  // myLocDataView_: persistent non-const view of my local data
  Teuchos::ArrayRCP<sc_t> myLocDataView_;
  // myLocDataView_: persistent const view of my local data
  Teuchos::ArrayRCP<const sc_t> myLocDataCView_;

};//end class

}}//end namespace rompp::core
#endif
#endif
