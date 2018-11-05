
#ifdef HAVE_TRILINOS
#ifndef CORE_MULTIVECTOR_CONCRETE_MULTIVECTOR_DISTRIBUTED_TPETRA_HPP_
#define CORE_MULTIVECTOR_CONCRETE_MULTIVECTOR_DISTRIBUTED_TPETRA_HPP_

#include "../../shared_base/core_container_base.hpp"
#include "../../shared_base/core_container_distributed_mpi_base.hpp"
#include "../../shared_base/core_container_distributed_trilinos_base.hpp"
#include "../base/core_multi_vector_distributed_base.hpp"

namespace rompp{ namespace core{
  
template <typename wrapped_type>
class MultiVector<wrapped_type,
     typename
     std::enable_if<
       meta::is_multi_vector_tpetra<
	    wrapped_type>::value
       >::type
     >
  : public ContainerBase< MultiVector<wrapped_type>, wrapped_type >,
    public MultiVectorDistributedBase< MultiVector<wrapped_type> >/*,
    public ContainerDistributedTrilinosBase< MultiVector<wrapped_type>, 
              typename details::traits<MultiVector<wrapped_type>>::data_map_t >, 
    public ContainerDistributedMpiBase< MultiVector<wrapped_type>, 
    typename details::traits<MultiVector<wrapped_type>>::communicator_t >*/{

private:
  using this_t = MultiVector<wrapped_type>;
  using sc_t = typename details::traits<this_t>::scalar_t;
  using LO_t = typename details::traits<this_t>::local_ordinal_t;
  using GO_t = typename details::traits<this_t>::global_ordinal_t;
  using wrap_t = typename details::traits<this_t>::wrapped_t;
  using map_t = typename details::traits<this_t>::data_map_t;
  using mpicomm_t = typename details::traits<this_t>::communicator_t;

public:
  MultiVector() = delete;

  explicit MultiVector(Teuchos::RCP<const map_t> mapobj, GO_t numVectors)
    : data_(mapobj, numVectors){}

  explicit MultiVector(const wrap_t & other)
    : data_(other){}

  ~MultiVector() = default;
  
private:

  wrap_t const * dataImpl() const{
    return &data_;
  }

  wrap_t * dataImpl(){
    return &data_;
  }

  bool emptyImpl() const{
    if (this->globalNumVectors()==0)
      return true;
    else
      return this->globalLength()==0  ? true : false;
  }

  void setZeroImpl() {
    data_.putScalar(static_cast<sc_t>(0));
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
    
  // void replaceGlobalValueImpl(GO_t globalRowIndex,
  // 			      GO_t vectorIndex,
  // 			      sc_t value){
  //   data_.replaceGlobalValue(globalRowIndex, vectorIndex, value);
  // }
  
  void scaleImpl(sc_t factor){
    data_.scale(factor);
  }
  
private:
  friend ContainerBase< this_t, wrapped_type >;
  friend MultiVectorDistributedBase< this_t >;
  // friend ContainerDistributedMpiBase< this_t, mpicomm_t >;
  // friend ContainerDistributedTrilinosBase< this_t, map_t >;
  
private:
  wrap_t data_;

};//end class

}}//end namespace rompp::core

#endif
#endif

