/*
//@HEADER
// ************************************************************************
//
// ops_mvec_dot_mvec.hpp
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
#ifndef OPS_SRC_OPS_TPETRA_BLOCK_MULTI_VECTOR_DOT_MULTI_VECTOR_HPP_
#define OPS_SRC_OPS_TPETRA_BLOCK_MULTI_VECTOR_DOT_MULTI_VECTOR_HPP_

namespace pressio{ namespace ops{

/*
 * multi_vector dot multi vector
 */

// result is stored into an Eigen dense matrix wrapper
template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_tpetra_block<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value 
    > * = nullptr
  >
void dot(const mvec_t & mvA, const mvec_t & mvB, result_t & C)
{
  static_assert(containers::meta::are_scalar_compatible<mvec_t, result_t>::value,
    "Types are not scalar compatible");

  // how many vectors are in mvA and mvB
  const auto numVecsA = mvA.globalNumVectors();
  const auto numVecsB = mvB.globalNumVectors();
  auto mvA_v = mvA.data()->getMultiVectorView();
  auto mvB_v = mvB.data()->getMultiVectorView();

  // compute dot between every column of A with every col of B
  for (std::size_t i=0; i<(std::size_t)numVecsA; i++)
  {
    // colI is a Teuchos::RCP<Vector<...>>
    const auto colI = mvA_v.getVector(i);
    for (std::size_t j=0; j<(std::size_t)numVecsB; j++)
    {
      const auto colJ = mvB_v.getVector(j);
      C(i,j) = colI->dot(*colJ);
    }
  }
}


// result is stored into an Eigen dense DYNAMIC matrix wrapper
template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra_block<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value and
    ::pressio::containers::details::traits<result_t>::is_dynamic 
    > * = nullptr
  >
result_t dot(const mvec_t & mvA, const mvec_t & mvB)
{
  const auto numVecsA = mvA.globalNumVectors();
  const auto numVecsB = mvB.globalNumVectors();
  result_t C(numVecsA, numVecsB);
  dot(mvA, mvB, C);
  return C;
}



// C += mvA^T mvB
template <
  typename mvec_t,
  typename expr_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_tpetra_block<mvec_t>::value and
    ::pressio::containers::meta::is_expression<expr_t>::value and
    ::pressio::containers::meta::are_scalar_compatible<mvec_t, expr_t>::value and
    ::pressio::containers::meta::is_matrix_wrapper_eigen<
      typename ::pressio::containers::details::traits<expr_t>::data_t
      >::value
    > * = nullptr
  >
void updateWithDot(const mvec_t & mvA, const mvec_t & mvB, expr_t & C)
{
  throw std::runtime_error("Error: updateWithDot for tpetrablock not yet supported");
  // // how many vectors are in mvA and mvB
  // const auto numVecsA = mvA.globalNumVectors();
  // const auto numVecsB = mvB.globalNumVectors();
  // auto mvA_v = mvA.data()->getMultiVectorView();
  // auto mvB_v = mvB.data()->getMultiVectorView();

  // // compute dot between every column of A with every col of B
  // for (std::size_t i=0; i<(std::size_t)numVecsA; i++){
  //   // colI is a Teuchos::RCP<Vector<...>>
  //   const auto colI = mvA_v.getVector(i);
  //   for (std::size_t j=0; j<(std::size_t)numVecsB; j++){
  //     const auto colJ = mvB_v.getVector(j);
  //     C(i,j) += colI->dot(*colJ);
  //   }
  // }
}


/* ----------------------------------------------
 * result_t = an expression
 ---------------------------------------------- */
template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_tpetra_block<mvec_t>::value and
    ::pressio::containers::meta::is_expression<result_t>::value and
    ::pressio::containers::meta::are_scalar_compatible<mvec_t, result_t>::value
    > * = nullptr
  >
void dot(const mvec_t & mvA, const mvec_t & mvB, result_t & C)
{
  throw std::runtime_error("Error container::ops::dot operation between tpetra_block, tpetra_block, and putting result into expression not yet supported");
}

}}//end namespace pressio::ops
#endif
#endif
