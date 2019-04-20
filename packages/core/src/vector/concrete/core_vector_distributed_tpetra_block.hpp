
#ifdef HAVE_TRILINOS
#ifndef CORE_VECTOR_CONCRETE_VECTOR_DISTRIBUTED_TPETRA_BLOCK_HPP_
#define CORE_VECTOR_CONCRETE_VECTOR_DISTRIBUTED_TPETRA_BLOCK_HPP_

#include "../../shared_base/core_container_base.hpp"
#include "../../shared_base/core_container_distributed_mpi_base.hpp"
#include "../../shared_base/core_container_distributed_trilinos_base.hpp"
#include "../../shared_base/core_container_resizable_base.hpp"
#include "../base/core_vector_distributed_base.hpp"

namespace rompp{ namespace core{

template <typename wrapped_type>
class Vector<
  wrapped_type,
  ::rompp::mpl::enable_if_t<
    meta::is_vector_tpetra_block<
      wrapped_type
      >::value
    >
  >
  : public ContainerBase< Vector<wrapped_type>, wrapped_type >,
    public VectorDistributedBase< Vector<wrapped_type> >,
    public ContainerDistributedTrilinosBase< Vector<wrapped_type>,
     typename details::traits<Vector<wrapped_type>>::data_map_t >
{

  using this_t		= Vector<wrapped_type>;
  using der_t		= this_t;
  using sc_t		= typename details::traits<this_t>::scalar_t;
  using LO_t		= typename details::traits<this_t>::local_ordinal_t;
  using GO_t		= typename details::traits<this_t>::global_ordinal_t;
  using device_t	= typename details::traits<this_t>::device_t;
  using wrap_t		= typename details::traits<this_t>::wrapped_t;
  using map_t		= typename details::traits<this_t>::data_map_t;
  using mpicomm_t	= typename details::traits<this_t>::communicator_t;

public:
  Vector() = delete;

  // explicit Vector(const map_t & mapO,
  // 		  const LO_t blockSize)
  //   : data_(mapO, blockSize){}


  /* Block MV/V still missing a copy constructor,
   * see https://github.com/trilinos/Trilinos/issues/4627
   * so for now we construct this one using other's data */

  explicit Vector(const wrap_t & vecobj)
    : data_( *vecobj.getMap(),
  	     vecobj.getBlockSize()){
    // just a trick to copy data
    data_.update(constants::one<sc_t>(),
		 vecobj,
		 constants::zero<sc_t>());
  }

  // delegate (for now) to the one above
  Vector(this_t const & other)
    : Vector(*other.data()){}

  ~Vector() = default;

private:

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

  void setZeroImpl(){
    data_.putScalar( ::rompp::core::constants::zero<sc_t>() );
    // putScalar doesn't sync afterwards, so we have to sync manually.
    this->needSync();
  }

  void putScalarImpl(sc_t value) {
    data_.putScalar(value);
    // putScalar doesn't sync afterwards, so we have to sync manually.
    this->needSync();
  }

  GO_t globalSizeImpl() const {
    return this->getDataMap().getGlobalNumElements();
  }

  LO_t localSizeImpl() const {
    return this->getDataMap().getNodeNumElements();
  }

  bool emptyImpl() const{
    // TODO: not sure this is great way to check
    return this->globalSize()==0 ? true : false;
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
  friend ContainerDistributedTrilinosBase< this_t, map_t >;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace rompp::core
#endif
#endif
