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

#ifndef CONTAINERS_MULTI_VECTOR_CONCRETE_CONTAINERS_MULTI_VECTOR_SHAREDMEM_KOKKOS_HPP_
#define CONTAINERS_MULTI_VECTOR_CONCRETE_CONTAINERS_MULTI_VECTOR_SHAREDMEM_KOKKOS_HPP_

namespace pressio{ namespace containers{

template <typename wrapped_type>
class MultiVector<
  wrapped_type,
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_admissible_as_multi_vector_kokkos<wrapped_type>::value
    >
  >
{
public:
  using this_t = MultiVector<wrapped_type>;
  using traits = details::traits<this_t>;
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

  explicit MultiVector(const wrap_t & src)
    : data_{src.label(), src.extent(0), src.extent(1)}{
    Kokkos::deep_copy(data_, src);
  }

  MultiVector(wrap_t && src)
    : data_(std::move(src)){}

  MultiVector(const std::string & label, size_t e1, size_t e2)
    : data_{label, e1, e2}{}

  // copy constr
  MultiVector(const MultiVector & other)
    : data_{other.data_.label(),
	    other.data_.extent(0),
	    other.data_.extent(1)}
  {
    Kokkos::deep_copy(data_, other.data_);
  }

  // delete copy assign to force usage of ops::deep_copy
  MultiVector & operator=(const MultiVector & other) = delete;
  //   if (&other != this){
  //     assert(this->length() == other.length());
  //     assert(this->numVectors() == other.numVectors());
  //     Kokkos::deep_copy(data_, *other.data());
  //   }
  //   return *this;
  // }

  // move cnstr and assign
  MultiVector(MultiVector && other) = default;
  MultiVector & operator=(MultiVector && other) = default;

  ~MultiVector() = default;

public:
  wrap_t const * data() const{
    return &data_;
  }
  wrap_t * data(){
    return &data_;
  }

  wrap_t dataCp(){
    return data_;
  }

  template< typename _wrapped_type = wrapped_type >
  mpl::enable_if_t<
    // todo: this is not entirely correct because this would work also
    // for UMV space, needs to be fixed
    std::is_same<typename mytraits::memory_space, Kokkos::HostSpace>::value,
    sc_t &>
  operator () (ord_t i, ord_t j){
    assert(i < this->extent(0));
    assert(j < this->extent(1));
    return data_(i,j);
  };

  template< typename _wrapped_type = wrapped_type >
  mpl::enable_if_t<
    // todo: not entirely correct because this would work also
    // for UMV space, needs to be fixed
    std::is_same<typename mytraits::memory_space, Kokkos::HostSpace>::value,
    sc_t const &>
  operator () (ord_t i, ord_t j) const{
    assert(i < this->extent(0));
    assert(j < this->extent(1));
    return data_(i, j);
  };

  ord_t extent(ord_t i) const {
    return data_.extent(i);
  }

  ord_t numVectors() const{
    return data_.extent(1);
  }

  ord_t length() const{
    return data_.extent(0);
  }

private:
  wrap_t data_ = {};
};

}}//end namespace pressio::containers
#endif  // CONTAINERS_MULTI_VECTOR_CONCRETE_CONTAINERS_MULTI_VECTOR_SHAREDMEM_KOKKOS_HPP_
