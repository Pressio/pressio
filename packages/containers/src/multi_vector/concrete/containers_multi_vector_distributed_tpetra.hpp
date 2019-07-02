
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_MULTIVECTOR_CONCRETE_MULTIVECTOR_DISTRIBUTED_TPETRA_HPP_
#define CONTAINERS_MULTIVECTOR_CONCRETE_MULTIVECTOR_DISTRIBUTED_TPETRA_HPP_

#include "../../shared_base/containers_container_base.hpp"
#include "../../shared_base/containers_container_distributed_mpi_base.hpp"
#include "../../shared_base/containers_container_distributed_trilinos_base.hpp"
#include "../base/containers_multi_vector_distributed_base.hpp"
#include <MatrixMarket_Tpetra.hpp>

namespace pressio{ namespace containers{

template <typename wrapped_type>
class MultiVector<wrapped_type,
     typename
     std::enable_if<
       meta::is_multi_vector_tpetra<
	    wrapped_type>::value
       >::type
     >
  : public ContainerBase< MultiVector<wrapped_type>, wrapped_type >,
    public MultiVectorDistributedBase< MultiVector<wrapped_type> >,
    public ContainerDistributedTrilinosBase< MultiVector<wrapped_type>,
     typename details::traits<MultiVector<wrapped_type>>::data_map_t>,
  public ContainerDistributedMpiBase< MultiVector<wrapped_type>,
     typename details::traits<MultiVector<wrapped_type>>::communicator_t>
{

private:
  using this_t = MultiVector<wrapped_type>;
  using sc_t = typename details::traits<this_t>::scalar_t;
  using LO_t = typename details::traits<this_t>::local_ordinal_t;
  using GO_t = typename details::traits<this_t>::global_ordinal_t;
  using device_t = typename details::traits<this_t>::device_t;
  using wrap_t = typename details::traits<this_t>::wrapped_t;
  using map_t = typename details::traits<this_t>::data_map_t;
  using mpicomm_t = typename details::traits<this_t>::communicator_t;

public:
  MultiVector() = delete;

  MultiVector(Teuchos::RCP<const map_t> mapobj, GO_t numVectors)
    : data_(mapobj, numVectors){}

  MultiVector(const map_t & mapobj, GO_t numVectors)
    : data_( Teuchos::rcpFromRef(mapobj), numVectors ){}

  explicit MultiVector(const wrap_t & other)
    // use the deep_copy constructor
    : data_(other, Teuchos::Copy){}

  MultiVector(const this_t & other)
    : data_(*other.data(), Teuchos::Copy){}

  ~MultiVector() = default;

  // compound assignment when type(b) = type(this)
  // this += b
  template <typename T,
  	    ::pressio::mpl::enable_if_t<
  	      std::is_same<T,this_t>::value> * = nullptr>
  this_t & operator+=(const T & other) {
    this->data_.update(1.0, *other.data(), 1.0 );
    return *this;
  }

  // copy assignment
  template <typename T,
  	    ::pressio::mpl::enable_if_t<
  	      std::is_same<T,this_t>::value> * = nullptr>
  this_t & operator=(const T & other){
    assert(this->localSize() == other.localSize());
    data_.assign( *other.data() );
    return *this;
  }

  void print(std::string tag) const{
    Tpetra::MatrixMarket::Writer<wrap_t>::writeDense
      (std::cout << std::setprecision(15), data_, tag, tag);
  }

 private:
  wrap_t const * dataImpl() const{
    return &data_;
  }

  wrap_t * dataImpl(){
    return &data_;
  }

  map_t const & getDataMapImpl() const{
    return *data_.getMap();
  }

  bool hasRowMapEqualToImpl(map_t const &othermap) const{
    return data_.getMap()->isSameAs(othermap);
  }

  Teuchos::RCP<const map_t> getRCPDataMapImpl() const{
    return data_.getMap();
  }

  bool emptyImpl() const{
    if (this->globalNumVectors()==0)
      return true;
    else
      return this->globalLength()==0  ? true : false;
  }

  void setZeroImpl() {
    data_.putScalar(static_cast<sc_t>(0));
    // putScalar doesn't sync afterwards, so we have to sync manually.
    this->needSync();
  }

  bool isDistributedGloballyImpl() const{
    return data_.isDistributed();
  }

  GO_t globalNumVectorsImpl() const{
    return data_.getNumVectors();
  }

  LO_t localNumVectorsImpl() const{
    return data_.getNumVectors();
  }

  GO_t globalLengthImpl() const {
    return data_.getGlobalLength();
  };

  LO_t localLengthImpl() const {
    return data_.getLocalLength();
  };

  mpicomm_t commImpl() const{
    return data_.getMap()->getComm();
  }

  void scaleImpl(sc_t factor){
    data_.scale(factor);
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
  friend MultiVectorDistributedBase< this_t >;
  friend ContainerDistributedTrilinosBase< this_t, map_t >;
  friend ContainerDistributedMpiBase< this_t, mpicomm_t >;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers

#endif
#endif /* HAVE_TRILINOS */
