/*
//@HEADER
// ************************************************************************
//
// ops_norms.hpp
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

#ifndef OPS_TPETRA_BLOCK_OPS_NORMS_HPP_
#define OPS_TPETRA_BLOCK_OPS_NORMS_HPP_

namespace pressio{ namespace ops{

template <typename T>
std::enable_if_t<
  ::pressio::is_vector_tpetra_block<T>::value
  // scalar compatibility
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value),
  typename decltype(std::declval<T>().getVectorView())::mag_type
  >
norm2(const T & a)
{
  static_assert(
    std::is_same<typename ::pressio::Traits<T>::scalar_type,
    typename decltype(std::declval<T>().getVectorView())::mag_type>::value,
    "Scalar and mag not same");

  /* workaround the non-constness of getVectorView,
   * which is supposed to be const but it is not */
  using vec_t = typename std::remove_const<T>::type;
  const auto a_v = const_cast<vec_t &>(a).getVectorView();
  return a_v.norm2();
}

template <typename T>
std::enable_if_t<
  ::pressio::is_vector_tpetra_block<T>::value
  // scalar compatibility
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value),
  typename decltype(std::declval<T>().getVectorView())::mag_type
  >
norm1(const T & a)
{
  static_assert(
    std::is_same<typename ::pressio::Traits<T>::scalar_type,
    typename decltype(std::declval<T>().getVectorView())::mag_type>::value,
    "Scalar and mag not same");

  /* workaround the non-constness of getVectorView,
   * which is supposed to be const but it is not */
  using vec_t = typename std::remove_const<T>::type;
  const auto a_v = const_cast<vec_t &>(a).getVectorView();
  return a_v.norm1();
}

template <
  typename T,
std::enable_if_t<
  ::pressio::is_expression_column_acting_on_tpetra_block<T>::value
  // scalar compatibility
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value),
  int > = 0
  >
auto norm1(const T & o)
{
  auto o_tp = impl::get_underlying_tpetra_object(o);
  return norm1(o_tp);
}

template <
  typename T,
std::enable_if_t<
  ::pressio::is_expression_column_acting_on_tpetra_block<T>::value
  // scalar compatibility
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value),
  int > = 0
  >
auto norm2(const T & o)
{
  auto o_tp = impl::get_underlying_tpetra_object(o);
  return norm2(o_tp);
}


}}//end namespace pressio::ops
#endif  // OPS_TPETRA_BLOCK_OPS_NORMS_HPP_
