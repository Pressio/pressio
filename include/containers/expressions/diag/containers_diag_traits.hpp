/*
//@HEADER
// ************************************************************************
//
// containers_diag_traits.hpp
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

#ifndef CONTAINERS_EXPRESSIONS_DIAG_CONTAINERS_DIAG_TRAITS_HPP_
#define CONTAINERS_EXPRESSIONS_DIAG_CONTAINERS_DIAG_TRAITS_HPP_

namespace pressio{ namespace containers{ namespace details{

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <typename matrix_type>
struct traits<
  ::pressio::containers::expressions::DiagExpr<matrix_type>,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::is_dense_matrix_wrapper_eigen<matrix_type>::value
    >
  >
  : public containers_shared_traits<
  typename details::traits<matrix_type>::wrapped_t,
  WrappedPackageIdentifier::Eigen, true, 1>
{

  static constexpr bool is_static = true;
  static constexpr bool is_dynamic  = !is_static;

  using wrapped_t = typename traits<matrix_type>::wrapped_t;
  using scalar_t  = typename traits<matrix_type>::scalar_t;
  using ordinal_t = typename traits<matrix_type>::ordinal_t;
  using size_t    = ordinal_t;

  // the reference type is conditional because the native expression
  // returns by value when object is const
  using reference_t =  typename std::conditional<
      std::is_const<matrix_type>::value, scalar_t , scalar_t &
      >::type;

  using const_reference_t = typename std::conditional<
    std::is_const<matrix_type>::value, scalar_t, scalar_t const &
  >::type;

  // type of the native expression
  using _native_expr_t = decltype(std::declval<wrapped_t>().diagonal( ) );
  using _const_native_expr_t=decltype(std::declval<const wrapped_t>().diagonal());

  using native_expr_t = typename std::conditional<
    std::is_const<matrix_type>::value,
    _const_native_expr_t,
    _native_expr_t
  >::type;

  using const_data_return_t = native_expr_t const *;
  using data_return_t = native_expr_t *;
};
#endif

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename matrix_type>
struct traits<
  ::pressio::containers::expressions::DiagExpr<matrix_type>,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::is_dense_matrix_wrapper_kokkos<matrix_type>::value
    >
  >
  : public containers_shared_traits<
  typename details::traits<matrix_type>::wrapped_t,WrappedPackageIdentifier::Kokkos,
  true, //true because kokkos is shared mem
  1
  >
{

  static constexpr bool is_static = true;
  static constexpr bool is_dynamic  = !is_static;
  using wrapped_t = typename traits<matrix_type>::wrapped_t;
  using scalar_t	= typename traits<matrix_type>::scalar_t;
  using execution_space = typename traits<matrix_type>::execution_space;
  using memory_space	= typename traits<matrix_type>::memory_space;
  using device_t	= typename traits<matrix_type>::device_t;
  using device_type	= typename traits<matrix_type>::device_t;
  using ordinal_t	= typename traits<matrix_type>::ordinal_t;
  using size_t		= ordinal_t;
  using reference_t	  = scalar_t &;
  using const_reference_t = scalar_t const &;

  using _native_expr_t	     = Kokkos::View<scalar_t*, Kokkos::LayoutStride>;
  using _const_native_expr_t = Kokkos::View<const scalar_t*, Kokkos::LayoutStride>;
  using native_expr_t = typename std::conditional<
    std::is_const<matrix_type>::value,
    _const_native_expr_t,
    _native_expr_t
  >::type;

  using const_data_return_t = native_expr_t const *;
  using data_return_t	    = native_expr_t *;
};
#endif

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <typename matrix_type>
struct traits<
  ::pressio::containers::expressions::DiagExpr<matrix_type>,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::is_rank2_tensor_wrapper_pybind<matrix_type>::value
    >
  >
  : public containers_shared_traits<
  typename details::traits<matrix_type>::wrapped_t,
  WrappedPackageIdentifier::Pybind,
  true, //true because kokkos is shared mem
  1
  >
{
  static constexpr bool is_static = true;
  static constexpr bool is_dynamic  = !is_static;
  using wrapped_t = typename traits<matrix_type>::wrapped_t;
  using scalar_t  = typename traits<matrix_type>::scalar_t;
  using ordinal_t = typename traits<matrix_type>::ordinal_t;
  using size_t    = ordinal_t;
  using reference_t =  scalar_t &;
  using const_reference_t = scalar_t const &;

  // // type of the native expression
  // // this should be a pybind11::array that views the diagonal of the matrix
  // // for time being we make this of the same type as the wrapped_type
  // // even if this is not fully right I think.
  // // it seems to be doing the right thing.
  // using native_expr_t = wrapped_t;
  // using const_data_return_t = native_expr_t const *;
  // using data_return_t = native_expr_t *;
};
#endif

}}}//end namespace pressio::containers::details
#endif  // CONTAINERS_EXPRESSIONS_DIAG_CONTAINERS_DIAG_TRAITS_HPP_
