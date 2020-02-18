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
#ifndef CONTAINERS_SRC_OPS_EPETRA_MULTI_VECTOR_DOT_SELF_HPP_
#define CONTAINERS_SRC_OPS_EPETRA_MULTI_VECTOR_DOT_SELF_HPP_

namespace pressio{ namespace containers{ namespace ops{

/*
 * A dot A
 * multi_vector dot self: equivalent to doing A^T A
 */

// store result into eigen dense matrix wrapper
template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_epetra<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value 
    > * = nullptr
  >
void dot_self(const mvec_t & A, result_t & C)
{
  static_assert(containers::meta::are_scalar_compatible<mvec_t, result_t>::value,
    "Types are not scalar compatible");

  // how many vectors are in A
  const auto numVecsA = A.numVectors();
  const auto & Adata = *A.data();
  assert(C.extent(0) == numVecsA);
  assert(C.extent(1) == numVecsA);

  // A dot A = A^T*A, which yields a symmetric matrix
  // only need to compute half and fill remaining entries accordingly
  for (auto i=0; i<numVecsA; i++){
    for (auto j=i; j<numVecsA; j++){
      Adata(i)->Dot( *(Adata(j)), &C(i,j) );
      // fill the lower triangular part
      C(j,i) = C(i,j);
    }
  }
}

// return result in the form of eigen dense matrix wrapper
template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_epetra<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value and
    ::pressio::containers::details::traits<result_t>::is_dynamic
   > * = nullptr
  >
result_t dot_self(const mvec_t & mvA)
{
  static_assert(containers::meta::are_scalar_compatible<mvec_t, result_t>::value,
    "Types are not scalar compatible");
  
  const auto numVecsA = mvA.numVectors();
  result_t C(numVecsA, numVecsA);
  dot_self(mvA, C);
  return C;
}


}}}//end namespace pressio::containers::ops
#endif
#endif
