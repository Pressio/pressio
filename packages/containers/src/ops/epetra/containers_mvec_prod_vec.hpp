/*
//@HEADER
// ************************************************************************
//
// containers_mvec_prod_vec.hpp
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
#ifndef CONTAINERS_SRC_OPS_EPETRA_MULTI_VECTOR_PROD_VECTOR_HPP_
#define CONTAINERS_SRC_OPS_EPETRA_MULTI_VECTOR_PROD_VECTOR_HPP_

namespace pressio{ namespace containers{ namespace ops{

/* epetra multi_vector prod vector */


// begin namespace pressio::containers::ops::impl
namespace impl{
template <typename mvec_type, typename operand_type>
void _product_epetra_mv_sharedmem_vec(const mvec_type & mvA,
				      const operand_type & b,
				      containers::Vector<Epetra_Vector> & C){
  //zero out result
  ::pressio::containers::ops::set_zero(C);
  // how many vectors are in mvA
  const auto numVecs = mvA.numVectorsGlobal();
  // size of b
  assert(size_t(numVecs) == size_t(b.extent(0)));
  // the data map of the multivector
  const auto mvMap = mvA.data()->Map();
  // my number of rows
  const auto myNrows = mvMap.NumMyElements();

  // loop
  for (size_t i=0; i<(size_t)myNrows; i++){
    for (size_t j=0; j<(size_t)numVecs; j++){
      C[i] += mvA(i,j) * b[j];
    }
  }
}
}//end namespace pressio::containers::ops::impl



/* -------------------------------------------------------------------
 * specialize for epetra mv wrapper operating on a sharedmem vector wrapper
 *-------------------------------------------------------------------*/
template <
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
    (containers::meta::is_vector_wrapper_eigen<vec_type>::value or
     containers::meta::is_dense_vector_wrapper_teuchos<vec_type>::value)
    > * = nullptr
  >
void product(const mvec_type & mvA,
	     const vec_type & vecB,
	     containers::Vector<Epetra_Vector> & C)
{
  static_assert(containers::meta::wrappers_have_same_scalar<mvec_type, vec_type>::value,
    "Types are not scalar compatible");

  ::pressio::containers::ops::impl::_product_epetra_mv_sharedmem_vec(mvA, vecB, C);
}

// return result Epetra wrapper
template <
  typename mvec_type, typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
    (containers::meta::is_vector_wrapper_eigen<vec_type>::value or
     containers::meta::is_dense_vector_wrapper_teuchos<vec_type>::value)
    > * = nullptr
  >
containers::Vector<Epetra_Vector>
product(const mvec_type & mvA, const vec_type & vecB) 
{
  static_assert(containers::meta::wrappers_have_same_scalar<mvec_type, vec_type>::value,
    "Types are not scalar compatible");

  // here, mvA is distrubted, but vecB is NOT.
  // we interpret this as a linear combination of vectors

  // the data map of the multivector
  const auto mvMap = mvA.data()->Map();
  // result is an Epetra Vector with same distribution of mvA
  using res_t = containers::Vector<Epetra_Vector>;
  res_t c(mvMap);
  product(mvA, vecB, c);
  return c;
}



/* -------------------------------------------------------------------
 * specialize for epetra mv wrapper operating on an expression
 *-------------------------------------------------------------------*/
template <
  typename mvec_type,
  typename expr_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
    ::pressio::containers::meta::is_expression<expr_type>::value
    > * = nullptr
  >
void product(const mvec_type & mvA,
	     const expr_type & exprObj,
	     containers::Vector<Epetra_Vector> & C)
{
  static_assert(containers::meta::wrappers_have_same_scalar<mvec_type, expr_type>::value,
    "Types are not scalar compatible");
  ::pressio::containers::ops::impl::_product_epetra_mv_sharedmem_vec(mvA, exprObj, C);
}

template <
  typename mvec_type,
  typename expr_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
    ::pressio::containers::meta::is_expression<expr_type>::value
    > * = nullptr
  >
containers::Vector<Epetra_Vector>
product(const mvec_type & mvA, const expr_type & exprObj) 
{
  static_assert(containers::meta::wrappers_have_same_scalar<mvec_type, expr_type>::value,
    "Types are not scalar compatible");

  const auto mvMap = mvA.data()->Map();
  using res_t = containers::Vector<Epetra_Vector>;
  res_t c(mvMap);
  product(mvA, exprObj, c);
  return c;
}


}}}//end namespace pressio::containers::ops
#endif
#endif
