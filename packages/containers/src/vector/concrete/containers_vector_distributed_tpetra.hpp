
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_VECTOR_CONCRETE_VECTOR_DISTRIBUTED_TPETRA_HPP_
#define CONTAINERS_VECTOR_CONCRETE_VECTOR_DISTRIBUTED_TPETRA_HPP_

#include "../../shared_base/containers_container_base.hpp"
#include "../../shared_base/containers_container_distributed_mpi_base.hpp"
#include "../../shared_base/containers_container_distributed_trilinos_base.hpp"
#include "../../shared_base/containers_container_subscriptable_base.hpp"
#include "../../shared_base/containers_container_resizable_base.hpp"
#include "../base/containers_vector_distributed_base.hpp"
#include <MatrixMarket_Tpetra.hpp>

namespace pressio{ namespace containers{

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
    /*public ContainerDistributedMpiBase< Vector<wrapped_type>,
     typename details::traits<Vector<wrapped_type>>::communicator_t >,*/
    public ContainerDistributedTrilinosBase< Vector<wrapped_type>,
     typename details::traits<Vector<wrapped_type>>::data_map_t >/*,
    public ContainerResizableBase< Vector<wrapped_type>, 1>,
    public ContainerSubscriptable1DBase< Vector<wrapped_type>,
     typename details::traits<Vector<wrapped_type>>::scalar_t,
     typename details::traits<Vector<wrapped_type>>::local_ordinal_t>*/{

  using this_t = Vector<wrapped_type>;
  using sc_t = typename details::traits<this_t>::scalar_t;
  using LO_t = typename details::traits<this_t>::local_ordinal_t;
  using GO_t = typename details::traits<this_t>::global_ordinal_t;
  using device_t = typename details::traits<this_t>::device_t;
  using der_t = this_t;
  using wrap_t = typename details::traits<this_t>::wrapped_t;
  using map_t = typename details::traits<this_t>::data_map_t;
  using mpicomm_t = typename details::traits<this_t>::communicator_t;

public:
  Vector() = delete;

  explicit Vector(const wrap_t & vecobj)
    // use the deep_copy constructor
    : data_(vecobj, Teuchos::Copy){}

  explicit Vector(Teuchos::RCP<const map_t> mapO)
    : data_(mapO){}

  Vector(this_t const & other)
    : data_(*other.data(), Teuchos::Copy){}

  ~Vector() = default;

public:

  // copy assignment
  this_t & operator=(const this_t & other){
    assert(this->localSize() == other.localSize());
    data_.assign( *other.data() );
    return *this;
  }

  // compound assignment when type(b) = type(this)
  // this += b
  this_t & operator+=(const this_t & other) {
    this->data_.update(1.0, *other.data(), 1.0 );
    return *this;
  }

  // compound assignment when type(b) = type(this)
  // this -= b
  this_t & operator-=(const this_t & other) {
    this->data_.update(-1.0, *other.data(), 1.0 );
    return *this;
  }

  void print(std::string tag) const{
    Tpetra::MatrixMarket::Writer<wrap_t>::writeDense
      (std::cout << std::setprecision(15), data_, tag, tag);
  }

private:

//   void matchLayoutWithImpl(const der_t & other){
//     data_.ReplaceMap( other.getDataMap() );
//   }

//   mpicomm_t const & commCRefImpl() const{
//     return data_.getMap()->getComm();
//   }

  map_t const & getDataMapImpl() const{
    return *data_.getMap();
  }

  wrap_t const * dataImpl() const{
    return &data_;
  }

  wrap_t * dataImpl(){
    return &data_;
  }

  wrap_t dataCpImpl(){
    return data_;
  }

  bool isDistributedGloballyImpl() const{
    return data_.isDistributed();
  }

  void setZeroImpl(){
    data_.putScalar(static_cast<sc_t>(0));
    // putScalar doesn't sync afterwards, so we have to sync manually.
    this->needSync();
  }

  bool emptyImpl() const{
    return this->globalSize()==0 ? true : false;
  }

  void putScalarImpl(sc_t value) {
    data_.putScalar(value);
    // putScalar doesn't sync afterwards, so we have to sync manually.
    this->needSync();
  }

  GO_t globalSizeImpl() const {
    return data_.getGlobalLength();
  }

  LO_t localSizeImpl() const {
    return data_.getLocalLength();
  }

private:
  void needSync(){
    if (data_.template need_sync<Kokkos::HostSpace>())
      data_.template sync<Kokkos::HostSpace> ();
    else if (data_.template need_sync<device_t>())
      data_.template sync<device_t> ();
  }

private:
  friend ContainerBase< this_t, wrapped_type >;
  friend VectorDistributedBase< this_t >;
  // friend ContainerDistributedMpiBase< this_t, mpicomm_t >;
  friend ContainerDistributedTrilinosBase< this_t, map_t >;
  // friend ContainerResizableBase< this_t, 1>;
  // friend ContainerSubscriptable1DBase< this_t, sc_t, LO_t>;

private:
  wrap_t data_ = {};

  // // myLocDataView_: persistent non-const view of my local data
  // Teuchos::ArrayRCP<sc_t> myLocDataView_;
  // // myLocDataView_: persistent const view of my local data
  // Teuchos::ArrayRCP<const sc_t> myLocDataCView_;

};//end class

}}//end namespace pressio::containers
#endif
#endif
