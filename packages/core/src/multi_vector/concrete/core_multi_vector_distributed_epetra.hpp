
#ifndef CORE_MULTIVECTOR_CONCRETE_MULTIVECTOR_DISTRIBUTED_EPETRA_HPP_
#define CORE_MULTIVECTOR_CONCRETE_MULTIVECTOR_DISTRIBUTED_EPETRA_HPP_

#include "../../shared_base/core_container_base.hpp"
#include "../base/core_multi_vector_distributed_base.hpp"
#include "../base/core_multi_vector_math_base.hpp"
#include "../../shared_base/core_container_distributed_mpi_base.hpp"
#include "../../shared_base/core_operators_base.hpp"
#include "../../shared_base/core_container_distributed_trilinos_base.hpp"

namespace core{
  
template <typename wrapped_type>
class MultiVector<wrapped_type,
     typename
     std::enable_if<
       meta::is_multi_vector_epetra<
	    wrapped_type>::value
       >::type
     >
  : public ContainerBase< MultiVector<wrapped_type>, wrapped_type >,
    public MultiVectorDistributedBase< MultiVector<wrapped_type> >,
    public MultiVectorMathBase< MultiVector<wrapped_type> >,
    public Subscripting2DOperatorsBase< MultiVector<wrapped_type>,
      typename details::traits<MultiVector<wrapped_type>>::scalar_t,
      typename details::traits<MultiVector<wrapped_type>>::local_ordinal_t,
      typename details::traits<MultiVector<wrapped_type>>::global_ordinal_t>,
    public ContainerDistributedTrilinosBase< MultiVector<wrapped_type>, 
              typename details::traits<MultiVector<wrapped_type>>::data_map_t >, 
    public ContainerDistributedMpiBase< MultiVector<wrapped_type>, 
      typename details::traits<MultiVector<wrapped_type>>::communicator_t >
{

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

  explicit MultiVector(const map_t & mapobj, GO_t numVectors)
    : data_(mapobj, numVectors){}

  explicit MultiVector(const wrap_t & other) : data_(other){}

  ~MultiVector() = default;

public:
  sc_t & operator()(LO_t irow, GO_t icol){
    assert(icol < this->globalNumVectors() );
    assert(irow < this->localLength() );
    return data_[icol][irow];
  }

  sc_t const & operator()(LO_t irow, GO_t icol)const{
    assert(icol < this->globalNumVectors() );
    assert(irow < this->localLength() );
    return data_[icol][irow];
  }
  
private:

  void matchLayoutWithImpl(const this_t & other){
    data_.ReplaceMap( other.getDataMap() );
  }

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
    data_.PutScalar(static_cast<sc_t>(0));
  }

  bool isDistributedGloballyImpl() const{
    data_.DistributedGlobal();
  }
  
  mpicomm_t const & commCRefImpl() const{
    return data_.Comm();
  }

  map_t const & getDataMapImpl() const{
    return data_.Map();
  }

  void replaceDataMapImpl(const map_t & mapObj){
    data_.ReplaceMap(mapObj);
  }

  GO_t globalNumVectorsImpl() const{
    return data_.NumVectors();
  }

  LO_t localNumVectorsImpl() const{
    // it is the same because epetra multivectors
    // are distributed on data, but each process owns
    // a part of each vector
    return data_.NumVectors();
  }
  
  GO_t globalLengthImpl() const {
    return data_.GlobalLength();
  };

  LO_t localLengthImpl() const {
    return data_.MyLength();
  };
    
  void replaceGlobalValueImpl(GO_t globalRowIndex,
			       GO_t vectorIndex,
			       sc_t value){
    data_.ReplaceGlobalValue(globalRowIndex, vectorIndex, value);
  }

  void scaleImpl(sc_t factor){
    data_.Scale(factor);
  }
  
private:
  friend ContainerBase< this_t, wrapped_type >;
  friend MultiVectorDistributedBase< this_t >;
  friend MultiVectorMathBase< this_t >;
  friend ContainerDistributedMpiBase< this_t, mpicomm_t >;
  friend ContainerDistributedTrilinosBase< this_t, map_t >;
  friend Subscripting2DOperatorsBase< this_t, sc_t, LO_t, GO_t>;
  
private:
  wrap_t data_;

};//end class

}//end namespace core
#endif
