/*
//@HEADER
// ************************************************************************
//
// containers_static_collection_impl.hpp
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

#ifndef CONTAINERS_COLLECTION_CONTAINERS_STATIC_COLLECTION_IMPL_HPP_
#define CONTAINERS_COLLECTION_CONTAINERS_STATIC_COLLECTION_IMPL_HPP_

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
#include <pybind11/pybind11.h>
#endif

namespace pressio{ namespace containers{ namespace impl{

/*
 * we need to be careful here for types with view semantics.
 * we create things by deep copying
 */

template<typename T, std::size_t n>
class IndexableStaticCollection
{
  static_assert
  (::pressio::containers::predicates::is_wrapper<T>::value,
    "You can only create a IndexableStaticCollection of a pressio container type.");

  static_assert
  (!::pressio::containers::predicates::is_expression<T>::value,
    "You cannot create a IndexableStaticCollection of pressio expressions.");

public:
  using value_type    = T;
  using data_type     = std::array<T, n>;
  using size_type     = std::size_t;

  static constexpr std::size_t size() {
    return n;
  }

  T & operator[](std::size_t i){
    assert( i<n );
    return data_[i];
  }

  T const & operator[](std::size_t i) const{
    assert( i<n );
    return data_[i];
  }

  T & operator()(std::size_t i){
    assert( i<n );
    return data_[i];
  }

  T const & operator()(std::size_t i) const{
    assert( i<n );
    return data_[i];
  }

public:
  template <
  typename _T = T,
  mpl::enable_if_t<std::is_default_constructible<_T>::value, int> = 0
  >
  IndexableStaticCollection(){};

  // constructor for n == 1
  template <std::size_t _n = n, mpl::enable_if_t<_n == 1, int> = 0>
  IndexableStaticCollection(T const & y)
    : data_{{y}}{}

  // constructor for n == 2
  template <std::size_t _n = n, mpl::enable_if_t<_n == 2, int> = 0>
  IndexableStaticCollection(T const & y)
    : data_{{y,y}}{}

  // constructor for n == 3
  template <std::size_t _n = n, mpl::enable_if_t<_n == 3, int> = 0>
  IndexableStaticCollection(T const & y)
    : data_{{y,y,y}}{}

  // constructor for n == 4
  template <std::size_t _n = n, mpl::enable_if_t<_n == 4, int> = 0>
  IndexableStaticCollection(T const & y)
    : data_{{y,y,y,y}}{}

  // constructor for n == 5
  template <std::size_t _n = n, mpl::enable_if_t<_n == 5, int> = 0>
  IndexableStaticCollection(T const & y)
    : data_{{y,y,y,y,y}}{}

  // constructor for n == 6
  template <std::size_t _n = n, mpl::enable_if_t<_n == 6, int> = 0>
  IndexableStaticCollection(T const & y)
    : data_{{y,y,y,y,y,y}}{}

  // constructor for n == 7
  template <std::size_t _n = n, mpl::enable_if_t<_n == 7, int> = 0>
  IndexableStaticCollection(T const & y)
    : data_{{y,y,y,y,y,y,y}}{}

  // copy cnstr
  IndexableStaticCollection(IndexableStaticCollection const & other) = default;
  // copy assignment
  IndexableStaticCollection & operator=(IndexableStaticCollection const & other) = default;
  // move cnstr
  IndexableStaticCollection(IndexableStaticCollection && other) = default;
  // move assignment
  IndexableStaticCollection & operator=(IndexableStaticCollection && other) = default;
  // destructor
  ~IndexableStaticCollection() = default;

private:
  data_type data_;
};

}}}//end namespace pressio::containers::impl
#endif  // CONTAINERS_COLLECTION_CONTAINERS_STATIC_COLLECTION_IMPL_HPP_
