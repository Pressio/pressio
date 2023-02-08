/*
//@HEADER
// ************************************************************************
//
// ops_min_max.hpp
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

#ifndef OPS_EPETRA_OPS_MIN_MAX_HPP_
#define OPS_EPETRA_OPS_MIN_MAX_HPP_

namespace pressio{ namespace ops{

template <typename T>
::pressio::mpl::enable_if_t<
  // TPL/container specific
    (::pressio::is_vector_epetra<T>::value
  || ::pressio::is_multi_vector_epetra<T>::value)
  // scalar compatibility
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value),
  typename ::pressio::Traits<T>::scalar_type
  >
max(T & obj)
{
  assert(::pressio::ops::extent(obj, 0) > 0);
  assert(!::pressio::is_multi_vector_epetra<T>::value
       || ::pressio::ops::extent(obj, 1) > 0);

  using sc_t = typename ::pressio::Traits<T>::scalar_type;
  std::vector<sc_t> x(obj.NumVectors());
  int ret = obj.MaxValue(x.data());
  const auto i = std::max_element(x.begin(), x.end());
  if (i == x.end()) {
    throw std::runtime_error("Can't take max() of no elements");
  }
  return *i;
}

template <typename T>
::pressio::mpl::enable_if_t<
  // TPL/container specific
    (::pressio::is_vector_epetra<T>::value
  || ::pressio::is_multi_vector_epetra<T>::value)
  // scalar compatibility
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value),
  typename ::pressio::Traits<T>::scalar_type
  >
min(const T & obj)
{
  assert(::pressio::ops::extent(obj, 0) > 0);
  assert(!::pressio::is_multi_vector_epetra<T>::value
       || ::pressio::ops::extent(obj, 1) > 0);

  using sc_t = typename ::pressio::Traits<T>::scalar_type;
  std::vector<sc_t> x(obj.NumVectors());
  int ret = obj.MinValue(x.data());
  const auto i = std::min_element(x.begin(), x.end());
  if (i == x.end()) {
    throw std::runtime_error("Can't take min() of no elements");
  }
  return *i;
}

}}//end namespace pressio::ops
#endif  // OPS_EPETRA_OPS_MIN_MAX_HPP_
