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

#ifndef OPS_EIGEN_OPS_LEVEL2_HPP_
#define OPS_EIGEN_OPS_LEVEL2_HPP_

namespace pressio{ namespace ops{

/*
 * y = beta * y + alpha*op(A)*x
*/

// Implementation notes.
// Eigen requires that:
// - scalar type of A, x and y is the same
// - alpha and beta types are same as scalar type
// - scalar can be constructed from double, eg. Scalar(1)
// Note: only applies to the overloads that inside use native Eigen operations

//-------------------------------
// op(A) = A
//-------------------------------
template <
  class A_type, class x_type, class y_type,
  class alpha_t, class beta_t
  >
::pressio::mpl::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && (::pressio::is_native_container_eigen<A_type>::value
   || ::pressio::is_expression_acting_on_eigen<A_type>::value)
  && (::pressio::is_vector_eigen<x_type>::value
   || ::pressio::is_expression_acting_on_eigen<x_type>::value)
  && (::pressio::is_vector_eigen<y_type>::value
   || ::pressio::is_expression_acting_on_eigen<y_type>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::nontranspose /*unused*/,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x,
	const beta_t & beta,
	y_type & y)
{
  assert( ::pressio::ops::extent(y, 0) == ::pressio::ops::extent(A, 0) );
  assert( ::pressio::ops::extent(x, 0) == ::pressio::ops::extent(A, 1) );

  using sc_t = typename ::pressio::Traits<y_type>::scalar_type;
  constexpr sc_t zero{0};
  const sc_t alpha_(alpha);
  const sc_t beta_(beta);
  const bool has_beta = beta_ != zero;
  auto & y_n = impl::get_native(y);
  const auto & A_n = impl::get_native(A);
  const auto & x_n = impl::get_native(x);
  if (alpha_ == zero) {
    ::pressio::ops::scale(y_n, beta_);
  } else {
    if (has_beta) { y_n = beta_ * y_n + alpha_ * A_n * x_n; }
    else { y_n = alpha_ * A_n * x_n; }
  }

}

//-------------------------------
// op(A) = A , construct result
//-------------------------------
template <
  class y_type, class A_type, class x_type, class alpha_t
  >
::pressio::mpl::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && (::pressio::is_native_container_eigen<A_type>::value
   || ::pressio::is_expression_acting_on_eigen<A_type>::value)
  && (::pressio::is_vector_eigen<x_type>::value
   || ::pressio::is_expression_acting_on_eigen<x_type>::value)
  && (::pressio::is_vector_eigen<y_type>::value
   || ::pressio::is_expression_acting_on_eigen<y_type>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value),
  y_type
  >
product(::pressio::nontranspose mode,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x)
{

  assert( ::pressio::ops::extent(x, 0) == ::pressio::ops::extent(A, 1) );

  y_type y( ::pressio::ops::extent(A, 0) );
  using sc_t = typename ::pressio::Traits<y_type>::scalar_type;
  product(mode, alpha, A, x, sc_t(0), y);
  return y;
}

//-------------------------------
// op(A) = A^T
//-------------------------------
template <
  class A_type, class x_type, class y_type,
  class alpha_t, class beta_t
  >
::pressio::mpl::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && (::pressio::is_native_container_eigen<A_type>::value
   || ::pressio::is_expression_acting_on_eigen<A_type>::value)
  && (::pressio::is_vector_eigen<x_type>::value
   || ::pressio::is_expression_acting_on_eigen<x_type>::value)
  && (::pressio::is_vector_eigen<y_type>::value
   || ::pressio::is_expression_acting_on_eigen<y_type>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::transpose /*unused*/,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x,
	const beta_t & beta,
	y_type & y)
{

  assert( ::pressio::ops::extent(y, 0) == ::pressio::ops::extent(A, 1) );
  assert( ::pressio::ops::extent(x, 0) == ::pressio::ops::extent(A, 0) );

  using sc_t = typename ::pressio::Traits<y_type>::scalar_type;
  constexpr sc_t zero{0};
  const sc_t alpha_(alpha);
  const sc_t beta_(beta);
  const bool has_beta = beta_ != zero;
  auto & y_n = impl::get_native(y);
  const auto & A_n = impl::get_native(A);
  const auto & x_n = impl::get_native(x);
  if (alpha_ == zero) {
    ::pressio::ops::scale(y_n, beta_);
  } else {
    if (has_beta) { y_n = beta_ * y_n + alpha_ * A_n.transpose() * x_n; }
    else { y_n = alpha_ * A_n.transpose() * x_n; }
  }
}

//-------------------------------
// op(A) = A^T, construct result
//-------------------------------
template <class y_type, class A_type, class x_type, class alpha_t>
::pressio::mpl::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && (::pressio::is_native_container_eigen<A_type>::value
   || ::pressio::is_expression_acting_on_eigen<A_type>::value)
  && (::pressio::is_vector_eigen<x_type>::value
   || ::pressio::is_expression_acting_on_eigen<x_type>::value)
  && (::pressio::is_vector_eigen<y_type>::value
   || ::pressio::is_expression_acting_on_eigen<y_type>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value),
  y_type
  >
product(::pressio::transpose mode,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x)
{

  assert( ::pressio::ops::extent(x, 0) == ::pressio::ops::extent(A, 0) );
  y_type y(::pressio::ops::extent(A, 1));

  using sc_t = typename ::pressio::Traits<y_type>::scalar_type;
  product(mode, alpha, A, x, sc_t{0}, y);
  return y;
}


}}//end namespace pressio::ops
#endif  // OPS_EIGEN_OPS_LEVEL2_HPP_
