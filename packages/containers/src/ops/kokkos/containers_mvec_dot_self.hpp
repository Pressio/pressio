/*
//@HEADER
// ************************************************************************
//
// containers_mvec_dot_self.hpp
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
#ifndef CONTAINERS_SRC_OPS_KOKKOS_MULTI_VECTOR_DOT_SELF_HPP_
#define CONTAINERS_SRC_OPS_KOKKOS_MULTI_VECTOR_DOT_SELF_HPP_

#include "../../multi_vector/containers_multi_vector_meta.hpp"
#include "KokkosBlas3_gemm.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*
 * A dot A
 * multi_vector dot self: equivalent to doing A^T A
 * returns a dense container kokkos wrapper matrix
 * in the form of a containers::Matrix< Kokkos::View<> >
 */


template <
  typename mvec_t,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_kokkos<mvec_t>::value
    > * = nullptr
  >
void dot_self(const mvec_t & A,
	      containers::Matrix<
		Kokkos::View<
		typename containers::details::traits<mvec_t>::scalar_t**,
		typename containers::details::traits<mvec_t>::layout,
		typename containers::details::traits<mvec_t>::execution_space
	       >
	      > & C)
{

  using sc_t = typename containers::details::traits<mvec_t>::scalar_t;
  constexpr auto zero = ::pressio::utils::constants::zero<sc_t>();
  constexpr auto one = ::pressio::utils::constants::one<sc_t>();
  const char ctA = 'T';
  const char ctB = 'N';
  KokkosBlas::gemm(&ctA, &ctB, one, *A.data(), *A.data(), zero, *C.data());
}


template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_kokkos<mvec_t>::value and
    containers::meta::is_matrix_wrapper_kokkos<result_t>::value and
    std::is_same<
      typename containers::details::traits<result_t>::execution_space,
      typename containers::details::traits<mvec_t>::execution_space
    >::value and
    std::is_same<
      typename containers::details::traits<result_t>::layout,
      typename containers::details::traits<mvec_t>::layout
    >::value and
    std::is_same<
      typename containers::details::traits<result_t>::scalar_t,
      typename containers::details::traits<mvec_t>::scalar_t
    >::value
    > * = nullptr
  >
result_t dot_self(const mvec_t & mvA)
{
  // using dm_t = Kokkos::View<sc_t**, layout, exe_space>;
  // using res_t = containers::Matrix<dm_t>;

  const auto numVecsA = mvA.numVectors();
  result_t C("dot_self_mat_res", numVecsA, numVecsA);
  dot_self(mvA, C);
  return C;
}


}}}//end namespace pressio::containers::ops
#endif
#endif
