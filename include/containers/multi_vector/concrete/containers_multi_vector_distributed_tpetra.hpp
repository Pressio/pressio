/*
//@HEADER
// ************************************************************************
//
// containers_multi_vector_distributed_tpetra.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
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

#ifndef CONTAINERS_MULTI_VECTOR_CONCRETE_CONTAINERS_MULTI_VECTOR_DISTRIBUTED_TPETRA_HPP_
#define CONTAINERS_MULTI_VECTOR_CONCRETE_CONTAINERS_MULTI_VECTOR_DISTRIBUTED_TPETRA_HPP_

#include <MatrixMarket_Tpetra.hpp>

namespace pressio{ namespace containers{

template <typename wrapped_type>
class MultiVector<
  wrapped_type,
  mpl::enable_if_t<
    ::pressio::containers::predicates::is_multi_vector_tpetra<
      wrapped_type>::value
    >
  >
{

public:
  using this_t = MultiVector<wrapped_type>;
  using traits = details::traits<this_t>;
  using sc_t = typename details::traits<this_t>::scalar_t;
  using LO_t = typename details::traits<this_t>::local_ordinal_t;
  using GO_t = typename details::traits<this_t>::global_ordinal_t;
  using device_t = typename details::traits<this_t>::device_t;
  using wrap_t = typename details::traits<this_t>::wrapped_t;
  using map_t = typename details::traits<this_t>::data_map_t;
  using mpicomm_t = typename details::traits<this_t>::communicator_t;

public:
  MultiVector() = default;

  explicit MultiVector(const wrap_t & other)
    // use the deep_copy constructor
    : data_(other, Teuchos::Copy){}

  MultiVector(Teuchos::RCP<const map_t> mapobj, GO_t numVectors)
    : data_(mapobj, numVectors){}

  MultiVector(const map_t & mapobj, GO_t numVectors)
    : data_( Teuchos::rcpFromRef(mapobj), numVectors ){}

  explicit MultiVector(wrap_t && other)
    : data_(std::move(other)){}

  // copy constr
  MultiVector(const MultiVector & other)
    : data_(*other.data(), Teuchos::Copy){}

  // delete copy assign to force usage of ops::deep_copy
  MultiVector & operator=(const MultiVector & other) = delete;
  //   if(&other != this){data_.assign( *other.data() );}
  //   return *this;
  // }

  // move cnstr
  MultiVector(MultiVector && other) = default;
  // move assignment
  MultiVector & operator=(MultiVector && other) = default;

  ~MultiVector() = default;

public:
  wrap_t const * data() const{
    return &data_;
  }

  wrap_t * data(){
    return &data_;
  }

  GO_t numVectors() const{
    return data_.getNumVectors();
  }

  GO_t numVectorsGlobal() const{
    return data_.getNumVectors();
  }

  LO_t numVectorsLocal() const{
    // it is the same because epetra multivectors
    // are distributed on data, but each process owns
    // a part of each vector
    return data_.getNumVectors();
  }

  // for distributed objects, extent return the global extent
  GO_t extent(std::size_t i) const{
    assert(i<=1);
    return (i==0) ? data_.getGlobalLength() : data_.getNumVectors();
  }

  LO_t extentLocal(std::size_t i) const{
    // each process owns all cols
    assert(i<=1);
    return (i==0) ? data_.getLocalLength() : data_.getNumVectors();
  }

private:
  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers

#endif  // CONTAINERS_MULTI_VECTOR_CONCRETE_CONTAINERS_MULTI_VECTOR_DISTRIBUTED_TPETRA_HPP_
