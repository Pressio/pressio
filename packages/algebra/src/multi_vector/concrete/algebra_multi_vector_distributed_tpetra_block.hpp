
#ifdef HAVE_TRILINOS
#ifndef ALGEBRA_MULTIVECTOR_CONCRETE_MULTIVECTOR_DISTRIBUTED_TPETRA_BLOCK_HPP_
#define ALGEBRA_MULTIVECTOR_CONCRETE_MULTIVECTOR_DISTRIBUTED_TPETRA_BLOCK_HPP_

#include "../../shared_base/algebra_container_base.hpp"
#include "../../shared_base/algebra_container_distributed_mpi_base.hpp"
#include "../../shared_base/algebra_container_distributed_trilinos_base.hpp"
#include "../base/algebra_multi_vector_distributed_base.hpp"

namespace rompp{ namespace algebra{

template <typename wrapped_type>
class MultiVector<
  wrapped_type,
  ::rompp::mpl::enable_if_t<
    meta::is_multi_vector_tpetra_block<
      wrapped_type
      >::value
    >
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

  /* Block MV/V still missing a copy constructor,
   * see https://github.com/trilinos/Trilinos/issues/4627
   * so for now we construct this one using other's data */

  explicit MultiVector(const wrap_t & other)
    : data_( *other.getMap(),
	     other.getBlockSize(),
	     other.getNumVectors()){
    // just a trick to copy data
    data_.update(::rompp::utils::constants::one<sc_t>(),
		 other,
		 ::rompp::utils::constants::zero<sc_t>());
  }

  // delegate (for now) to the one above
  MultiVector(const this_t & other)
    : MultiVector(*other.data()){}

  MultiVector(const map_t & map, LO_t blockSize, LO_t numVectors)
    : data_(map, blockSize, numVectors){}

  ~MultiVector() = default;

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

  void setZeroImpl() {
    data_.putScalar( ::rompp::utils::constants::zero<sc_t>() );
    // putScalar doesn't sync afterwards, so we have to sync manually.
    this->needSync();
  }

  GO_t globalNumVectorsImpl() const{
    return data_.getNumVectors();
  }

  LO_t localNumVectorsImpl() const{
    return data_.getNumVectors();
  }

  GO_t globalLengthImpl() const {
    return this->getDataMap().getGlobalNumElements();
  };

  LO_t localLengthImpl() const {
    return this->getDataMap().getNodeNumElements();
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
}}//end namespace rompp::algebra

#endif
#endif /* HAVE_TRILINOS */
