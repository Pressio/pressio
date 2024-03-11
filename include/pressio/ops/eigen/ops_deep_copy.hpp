/*
//@HEADER
// ************************************************************************
//
// ops_deep_copy.hpp
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

#ifndef OPS_EIGEN_OPS_DEEP_COPY_HPP_
#define OPS_EIGEN_OPS_DEEP_COPY_HPP_

namespace pressio{ namespace ops{

template<typename T1, typename T2>
std::enable_if_t<
  // common deep_copy constraints
  ::pressio::Traits<T1>::rank == ::pressio::Traits<T2>::rank
  // TPL/container specific
  && (::pressio::is_native_container_eigen<T1>::value
   || ::pressio::is_expression_acting_on_eigen<T1>::value)
  && (::pressio::is_native_container_eigen<T2>::value
   || ::pressio::is_expression_acting_on_eigen<T2>::value)
  >
deep_copy(T2 & dest, const T1 & src)
{
  assert((matching_extents<T2, T1>::compare(dest, src)));

  const auto src_n = impl::get_native(src);
  auto && dest_n = impl::get_native(dest);
  dest_n = src_n;
}

}}//end namespace pressio::ops
#endif  // OPS_EIGEN_OPS_DEEP_COPY_HPP_
