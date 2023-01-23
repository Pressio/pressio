/*
//@HEADER
// ************************************************************************
//
// ops_rank2_update.hpp
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

#ifndef OPS_KOKKOS_OPS_RANK2_UPDATE_HPP_
#define OPS_KOKKOS_OPS_RANK2_UPDATE_HPP_

#include<KokkosBlas1_axpby.hpp>

namespace pressio{ namespace ops{

//----------------------------------------------------------------------
// overloads for computing: MV = a * MV + b * MV1
// where MV is an kokkos multivector wrapper
//----------------------------------------------------------------------
template<typename T, typename T1, typename alpha_t, typename beta_t>
::pressio::mpl::enable_if_t<
  // rank-1 update common constraints
     ::pressio::Traits<T>::rank == 2
  && ::pressio::Traits<T1>::rank == 2
  // TPL/container specific
  && ::pressio::is_native_container_kokkos<T>::value
  && ::pressio::is_native_container_kokkos<T1>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<T, T1>::value
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value)
  && std::is_convertible<alpha_t, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<T>::scalar_type>::value
  >
update(T & mv, const alpha_t &a,
       const T1 & mv1, const beta_t &b)
{
  /* make sure we don't pass const objects to be modified.
     In kokkos it is legal to modify const views, not for pressio wrappers. */
  static_assert
    (!std::is_const<T1>::value,
     "cannot modify a const-qualified wrapper of a Kokkos view");

  assert(::pressio::ops::extent(mv, 0) == ::pressio::ops::extent(mv1, 0));
  assert(::pressio::ops::extent(mv, 1) == ::pressio::ops::extent(mv1, 1));

  using sc_t = typename ::pressio::Traits<T>::scalar_type;
  sc_t a_{a};
  sc_t b_{b};

  ::KokkosBlas::axpby(b_, mv1, a_, mv);
}

}}//end namespace pressio::ops
#endif  // OPS_KOKKOS_OPS_RANK2_UPDATE_HPP_
