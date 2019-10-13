/*
//@HEADER
// ************************************************************************
//
// containers_multi_vector_sharedmem_kokkos.hpp
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

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#ifndef CONTAINERS_MULTI_VECTOR_CONCRETE_VECTOR_SHAREDMEM_KOKKOS_HPP_
#define CONTAINERS_MULTI_VECTOR_CONCRETE_VECTOR_SHAREDMEM_KOKKOS_HPP_

#include "../../shared_base/containers_container_base.hpp"
#include "../base/containers_multi_vector_sharedmem_base.hpp"

namespace pressio{ namespace containers{

template <typename wrapped_type>
class MultiVector<
  wrapped_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_kokkos<wrapped_type>::value
    >
  >
  : public ContainerBase< MultiVector<wrapped_type>, wrapped_type >,
    public MultiVectorSharedMemBase<MultiVector<wrapped_type>>
{

  using this_t = MultiVector<wrapped_type>;
  using mytraits = typename details::traits<this_t>;
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename  mytraits::ordinal_t;
  using wrap_t = typename mytraits::wrapped_t;

  // Views have "view semantics." copy constructor and
  // operator= only do shallow copies.
  // Here, for the time being, we construct wrapper
  // of a view WITHOUT doing shallow copy.
  // We create a new object and deep_copy original.

public:
  MultiVector() = default;

  explicit MultiVector(const wrap_t src)
    : data_{src.label(), src.extent(0), src.extent(1)}{
    Kokkos::deep_copy(data_, src);
  }

  MultiVector(const std::string & label, size_t e1, size_t e2)
    : data_{label, e1, e2}
  {}

  MultiVector(const this_t & other)
    : data_{other.data_.label(),
	    other.data_.extent(0),
	    other.data_.extent(1)}{
    Kokkos::deep_copy(data_, other.data_);
  }

  ~MultiVector(){}

public:
  // copy assign implments copy semantics not view (for time being)
  this_t & operator=(const this_t & other){
    assert(this->length() == other.length());
    assert(this->numVectors() == other.numVectors());
    Kokkos::deep_copy(data_, *other.data());
    return *this;
  }

private:
  wrap_t const * dataImpl() const{
    return &data_;
  }
  wrap_t * dataImpl(){
    return &data_;
  }

  wrap_t dataCpImpl(){
    return data_;
  }

  ord_t numVectorsImpl() const{
    return data_.extent(1);
  }

  ord_t lengthImpl() const{
    return data_.extent(0);
  }

private:
  friend ContainerBase< this_t, wrapped_type >;
  friend MultiVectorSharedMemBase<this_t>;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers
#endif
#endif // PRESSIO_ENABLE_TPL_TRILINOS
