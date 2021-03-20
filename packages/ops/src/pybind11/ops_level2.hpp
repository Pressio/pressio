/*
//@HEADER
// ************************************************************************
//
// ops_level2.hpp
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

#ifndef OPS_PYBIND11_OPS_LEVEL2_HPP_
#define OPS_PYBIND11_OPS_LEVEL2_HPP_

namespace pressio{ namespace ops{

/*
 * y = beta * y + alpha*A*x
 * where y,A,x are tensor wrappers
*/
template < typename A_type, typename x_type, typename scalar_type, typename y_type>
::pressio::mpl::enable_if_t<
  containers::predicates::is_rank1_tensor_wrapper_pybind<y_type>::value and
  containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<A_type>::value and
  containers::predicates::is_rank1_tensor_wrapper_pybind<x_type>::value
>
product(::pressio::nontranspose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	y_type & y)
{
  static_assert
    (containers::predicates::are_scalar_compatible<A_type, x_type, y_type>::value,
     "Types are not scalar compatible");

  // NOTE: need to check if doing this import is expensive,
  // and assess whether we can use blas directly when we know
  // that objects involved are dense with not strange layout.
  pybind11::object pyblas_ = pybind11::module::import("scipy.linalg.blas");

  assert( y.extent(0) == A.extent(0) );
  assert( x.extent(0) == A.extent(1) );
  const auto & A_native = *A.data();
  const auto & x_native = *x.data();
  auto & y_native = *y.data();

  constexpr auto izero	    = ::pressio::utils::constants<int>::zero();
  constexpr auto ione	    = ::pressio::utils::constants<int>::one();
  constexpr auto transA	    = izero;
  constexpr auto overWritey = ione;
  pyblas_.attr("dgemv")(alpha, A_native, x_native, beta, y_native,
			izero, ione, izero, ione, transA, overWritey);
}

/*
 * y = beta * y + alpha*A^T*x
 * where y,A,x are tensor wrappers
*/
template < typename A_type, typename x_type, typename scalar_type, typename y_type>
::pressio::mpl::enable_if_t<
  containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<A_type>::value and
  containers::predicates::is_rank1_tensor_wrapper_pybind<x_type>::value and
  containers::predicates::is_rank1_tensor_wrapper_pybind<y_type>::value
  >
product(::pressio::transpose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	y_type & y)
{
  static_assert
    (containers::predicates::are_scalar_compatible<A_type, x_type, y_type>::value,
     "Types are not scalar compatible");

  // NOTE: need to check if doing this import is expensive,
  // and assess whether we can use blas directly when we know
  // that objects involved are dense with not strange layout.
  pybind11::object pyblas_ = pybind11::module::import("scipy.linalg.blas");

  assert( y.extent(0) == A.extent(1) );
  assert( x.extent(0) == A.extent(0) );
  const auto & AE = *A.data();
  const auto & xE = *x.data();
  auto & yE = *y.data();

  constexpr auto izero	    = ::pressio::utils::constants<int>::zero();
  constexpr auto ione	    = ::pressio::utils::constants<int>::one();
  constexpr auto transA	    = ione;
  constexpr auto overWritey = ione;
  pyblas_.attr("dgemv")(alpha, AE, xE, beta, yE,
			izero, ione, izero, ione, transA, overWritey);
}


/*
 * y = beta * y + alpha*A*x
 * where y,A are tensor wrappers
 * x is span expression
*/
template < typename A_type, typename x_type, typename scalar_type, typename y_type>
::pressio::mpl::enable_if_t<
  containers::predicates::is_rank1_tensor_wrapper_pybind<y_type>::value and
  containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<A_type>::value
>
product(::pressio::nontranspose mode,
	const scalar_type alpha,
	const A_type & A,
	const pressio::containers::expressions::SpanExpr<x_type> & x,
	const scalar_type beta,
	y_type & y)
{
  static_assert
    (containers::predicates::are_scalar_compatible<A_type, x_type, y_type>::value,
     "Types are not scalar compatible");

  assert( y.extent(0) == A.extent(0) );
  assert( x.extent(0) == A.extent(1) );

  const auto nArows = A.extent(0);
  const auto nAcols = A.extent(1);

  for (std::size_t i=0; i<nArows; ++i){
    y(i) = beta*y(i);
  }
  for (std::size_t j=0; j<nAcols; ++j){
    for (std::size_t i=0; i<nArows; ++i){
      y(i) += alpha * A(i,j) * x(j);
    }
  }
}

/*
 * y = beta * y + alpha*A^T*x
 * where A,x are tensor wrappers
 * y is span expression
*/
template < typename A_type, typename x_type, typename scalar_type, typename y_type>
::pressio::mpl::enable_if_t<
  containers::predicates::is_rank1_tensor_wrapper_pybind<x_type>::value and
  containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<A_type>::value
>
product(::pressio::transpose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	pressio::containers::expressions::SpanExpr<y_type> & y)
{
  static_assert
    (containers::predicates::are_scalar_compatible<A_type, x_type, y_type>::value,
     "Types are not scalar compatible");

  assert( y.extent(0) == A.extent(1) );
  assert( x.extent(0) == A.extent(0) );
  for (std::size_t i=0; i<A.extent(1); i++){
    y(i) = beta * y(i);
    for (std::size_t j=0; j<x.extent(0); j++){
      y(i) += alpha * A(j,i) * x(j);
    }
  }
}

}}//end namespace pressio::ops
#endif  // OPS_PYBIND11_OPS_LEVEL2_HPP_
