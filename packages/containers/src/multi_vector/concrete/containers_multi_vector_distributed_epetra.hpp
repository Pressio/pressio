/*
//@HEADER
// ************************************************************************
//
// containers_multi_vector_distributed_epetra.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the 
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions 
// are met:
//
// 1. Redistributions of source code must retain the above copyright 
// notice, this list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright 
// notice, this list of conditions and the following disclaimer in the 
// documentation and/or other materials provided with the distribution.
// 
// 3. Neither the name of the copyright holder nor the names of its 
// contributors may be used to endorse or promote products derived 
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_MULTIVECTOR_CONCRETE_MULTIVECTOR_DISTRIBUTED_EPETRA_HPP_
#define CONTAINERS_MULTIVECTOR_CONCRETE_MULTIVECTOR_DISTRIBUTED_EPETRA_HPP_

#include "../../shared_base/containers_container_base.hpp"
#include "../../shared_base/containers_container_distributed_mpi_base.hpp"
#include "../../shared_base/containers_container_distributed_trilinos_base.hpp"
#include "../base/containers_multi_vector_distributed_base.hpp"

namespace pressio{ namespace containers{

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

  MultiVector(const map_t & mapobj, GO_t numVectors)
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

  // compound assignment when type(b) = type(this)
  // this += b
  this_t & operator+=(const this_t & other) {
    this->data_.Update(1.0, *other.data(), 1.0 );
    return *this;
  }

  void print(std::string tag) const{
    ::pressio::utils::io::print_stdout(tag, utils::io::reset());
    data_.Print(std::cout);
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
    return data_.DistributedGlobal();
  }

  mpicomm_t const & commCRefImpl() const{
    return data_.Comm();
  }

  map_t const & getDataMapImpl() const{
    return data_.Map();
  }

  bool hasRowMapEqualToImpl(map_t const &othermap) const{
    return data_.Map().SameAs(othermap);
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
  friend ContainerDistributedMpiBase< this_t, mpicomm_t >;
  friend ContainerDistributedTrilinosBase< this_t, map_t >;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers

#endif /* CONTAINERS_MULTIVECTOR_CONCRETE_MULTIVECTOR_DISTRIBUTED_EPETRA_HPP_ */
#endif /* HAVE_TRILINOS */
