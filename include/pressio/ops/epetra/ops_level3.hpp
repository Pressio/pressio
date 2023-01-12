/*
//@HEADER
// ************************************************************************
//
// ops_level3.hpp
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

#ifndef OPS_EPETRA_OPS_LEVEL3_HPP_
#define OPS_EPETRA_OPS_LEVEL3_HPP_

namespace pressio{ namespace ops{

/*
 * C = beta * C + alpha*op(A)*op(B)
*/


/* -------------------------------------------------------------------
C = beta * C + alpha * A^T * B

A = epetra multivector
B = epetra multivector
C is an Eigen dense matrix
*-------------------------------------------------------------------*/
template <
  typename A_type, typename B_type, typename alpha_type, typename beta_type, typename C_type
  >
::pressio::mpl::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<B_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_multi_vector_epetra<A_type>::value
  && ::pressio::is_multi_vector_epetra<B_type>::value
  && ::pressio::is_dense_matrix_eigen<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, B_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  && std::is_convertible<alpha_type, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_type,  typename ::pressio::Traits<A_type>::scalar_type>::value
  >
product(::pressio::transpose /*unused*/,
	::pressio::nontranspose /*unused*/,
	const alpha_type alpha,
	const A_type & A,
	const B_type & B,
	const beta_type beta,
	C_type & C)
{

  const auto numVecsA = A.NumVectors();
  const auto numVecsB = B.NumVectors();
  assert( (std::size_t)::pressio::ops::extent(A,0) == (std::size_t)::pressio::ops::extent(B,0));
  assert( (std::size_t)::pressio::ops::extent(C,0) == (std::size_t)numVecsA );
  assert( (std::size_t)::pressio::ops::extent(C,1) == (std::size_t)numVecsB );

  using sc_t = typename ::pressio::Traits<A_type>::scalar_type;
  const auto zero = ::pressio::utils::Constants<sc_t>::zero();
  const sc_t alpha_(alpha);
  const sc_t beta_(beta);
  auto tmp = zero;
  // compute dot between every column of A with every col of B
  for (std::size_t i=0; i<(std::size_t)numVecsA; i++)
  {
    for (std::size_t j=0; j<(std::size_t)numVecsB; j++)
    {
      A(i)->Dot( *(B(j)), &tmp );
      tmp *= alpha_;
      C(i,j) = beta_ == zero ? zero : beta_ * C(i,j);
      C(i,j) += tmp;
    }
  }
}

/* -------------------------------------------------------------------
C = beta * C + alpha * A^T * B

A = multivector
B = multivector
C is an Eigen dense matrix returned by the function
*-------------------------------------------------------------------*/
template <
  typename C_type, typename A_type, typename B_type, typename alpha_type
  >
::pressio::mpl::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<B_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_multi_vector_epetra<A_type>::value
  && ::pressio::is_multi_vector_epetra<B_type>::value
  && ::pressio::is_dynamic_dense_matrix_eigen<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, B_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  && std::is_convertible<alpha_type, typename ::pressio::Traits<A_type>::scalar_type>::value,
  C_type
  >
product(::pressio::transpose modeA,
 ::pressio::nontranspose modeB,
 const alpha_type alpha,
 const A_type & A,
 const B_type & B)
{
  using sc_t = typename ::pressio::Traits<A_type>::scalar_type;
  constexpr auto zero = ::pressio::utils::Constants<sc_t>::zero();

  const auto numVecsA = A.NumVectors();
  const auto numVecsB = B.NumVectors();
  C_type C(numVecsA, numVecsB);
  product(modeA, modeB, alpha, A, B, zero, C);
  return C;
}


/* -------------------------------------------------------------------
special case A==B

C = beta * C + alpha * A^T * A

A = epetra multivector
C is an Eigen dense matrix
*-------------------------------------------------------------------*/
template <
  typename A_type, typename alpha_type, typename beta_type, typename C_type
  >
::pressio::mpl::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_multi_vector_epetra<A_type>::value
  && ::pressio::is_dense_matrix_eigen<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  && std::is_convertible<alpha_type, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_type,  typename ::pressio::Traits<A_type>::scalar_type>::value
  >
product(::pressio::transpose /*unused*/,
	::pressio::nontranspose /*unused*/,
	const alpha_type alpha,
	const A_type & A,
	const beta_type beta,
	C_type & C)
{

  assert( ::pressio::ops::extent(C, 0) == ::pressio::ops::extent(A, 1) );
  assert( ::pressio::ops::extent(C, 1) == ::pressio::ops::extent(A, 1) );

  // how many vectors are in A and B
  const int numVecsA = A.NumVectors();
  assert(C.rows() == numVecsA);
  assert(C.cols() == numVecsA);

  using sc_t = typename ::pressio::Traits<A_type>::scalar_type;
  constexpr auto zero = ::pressio::utils::Constants<sc_t>::zero();
  sc_t tmp = zero;
  const sc_t alpha_(alpha);
  const sc_t beta_(beta);
  const auto apply_beta = [beta_](sc_t c) -> sc_t {
    return beta_ == zero ? zero : beta_ * c;
  };

  // A dot A = A^T*A, which yields a symmetric matrix
  // only need to compute half and fill remaining entries accordingly
  for (int i=0; i< numVecsA; i++)
  {
    C(i,i) = apply_beta(C(i,i));
    A(i)->Dot( *(A(i)), &tmp );
    C(i,i) += alpha_*tmp;

    for (int j=i+1; j<numVecsA; j++)
    {
      C(i,j) = apply_beta(C(i,j));
      C(j,i) = apply_beta(C(j,i));

      A(i)->Dot( *(A(j)), &tmp );
      C(i,j) += alpha_*tmp;
      C(j,i) += alpha_*tmp;
    }
  }
}

template <
  typename C_type, typename A_type, typename alpha_type
  >
::pressio::mpl::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_multi_vector_epetra<A_type>::value
  && ::pressio::is_dynamic_dense_matrix_eigen<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  && std::is_convertible<alpha_type, typename ::pressio::Traits<A_type>::scalar_type>::value,
  C_type
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const alpha_type alpha,
	const A_type & A)
{
  using sc_t = typename ::pressio::Traits<A_type>::scalar_type;
  constexpr auto zero = ::pressio::utils::Constants<sc_t>::zero();
  C_type C(A.NumVectors(), A.NumVectors());
  product(modeA, modeB, alpha, A, zero, C);
  return C;
}

}}//end namespace pressio::ops
#endif  // OPS_EPETRA_OPS_LEVEL3_HPP_
