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

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#ifndef CONTAINERS_SRC_OPS_TPETRA_MULTI_VECTOR_DOT_SELF_HPP_
#define CONTAINERS_SRC_OPS_TPETRA_MULTI_VECTOR_DOT_SELF_HPP_

#include "../../multi_vector/containers_multi_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*
 * A dot A
 * multi_vector dot self: equivalent to doing A^T A
 */

// --------------------------------------------------------
// specialize when result is an eigen dense matrix wrapper
// --------------------------------------------------------
template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_tpetra<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<mvec_t, result_t>::value
    > * = nullptr
  >
void dot_self(const mvec_t & A, result_t & C)
{
  // how many vectors are in A and mvB
  const auto numVecsA = A.globalNumVectors();
  assert(C.rows() == numVecsA);
  assert(C.cols() == numVecsA);

  // A dot A = A^T*A, which yields a symmetric matrix
  // only need to compute half and fill remaining entries accordingly
  for (std::size_t i=0; i<(std::size_t)numVecsA; i++){
    // colI is a Teuchos::RCP<Vector<...>>
    auto colI = A.data()->getVector(i);
    for (std::size_t j=i; j<(std::size_t)numVecsA; j++){
      auto colJ = A.data()->getVector(j);
      C(i,j) = colI->dot(*colJ);
      C(j,i) = C(i,j);
    }
  }
}

template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_tpetra<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value and
    ::pressio::containers::details::traits<result_t>::is_dynamic and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<mvec_t, result_t>::value
    > * = nullptr
  >
result_t dot_self(const mvec_t & mvA)
{
  const auto numVecsA = mvA.globalNumVectors();
  result_t C(numVecsA, numVecsA);
  dot_self(mvA, C);
  return C;
}


// --------------------------------------------------------
// specialize when result is a kokkos view wrapper
// --------------------------------------------------------
template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_tpetra<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_kokkos<result_t>::value and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<mvec_t, result_t>::value and
    std::is_same<
      typename containers::details::traits<mvec_t>::device_t,
      typename containers::details::traits<result_t>::device_t
      >::value
    > * = nullptr
  >
void dot_self(const mvec_t & A, result_t & C)
{

  using scalar_t = typename ::pressio::containers::details::traits<mvec_t>::scalar_t;
  using map_t    = typename ::pressio::containers::details::traits<mvec_t>::data_map_t;
  using tpetra_mv_t = typename ::pressio::containers::details::traits<mvec_t>::wrapped_t;
  const auto indexBase = A.data()->getMap()->getIndexBase();
  const auto comm = A.data()->getMap()->getComm();
  // C should be symmetric
  assert( C.rows() == C.cols() );
  const auto n = C.rows();
  Teuchos::RCP<const map_t> replMap(new map_t(n, indexBase, comm, Tpetra::LocallyReplicated));
  // create multivector that views the Kokkos matrix
  tpetra_mv_t Cmv(replMap, *C.data());

  constexpr auto beta = ::pressio::utils::constants::zero<scalar_t>();
  constexpr auto alpha = ::pressio::utils::constants::one<scalar_t>();
  // do the operation C = A^T A
  Cmv.multiply(Teuchos::ETransp::TRANS, Teuchos::ETransp::NO_TRANS, alpha, *A.data(), *A.data(), beta);
}

template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_tpetra<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_kokkos<result_t>::value and
    ::pressio::containers::details::traits<result_t>::is_dynamic and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<mvec_t, result_t>::value
    > * = nullptr
  >
result_t dot_self(const mvec_t & mvA)
{
  const auto numVecsA = mvA.globalNumVectors();
  result_t C("dummyLabel", numVecsA, numVecsA);
  dot_self(mvA, C);
  return C;
}

}}}//end namespace pressio::containers::ops
#endif
#endif
