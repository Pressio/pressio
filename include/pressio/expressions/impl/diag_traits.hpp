/*
//@HEADER
// ************************************************************************
//
// diag_traits.hpp
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

#ifndef EXPRESSIONS_IMPL_DIAG_TRAITS_HPP_
#define EXPRESSIONS_IMPL_DIAG_TRAITS_HPP_

namespace pressio{ namespace expressions{ namespace impl{

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <typename MatrixType>
struct DiagTraits<
  DiagExpr<MatrixType>,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_dense_matrix_eigen<MatrixType>::value
    >
  >
  : public ::pressio::impl::EigenTraits<
      typename ::pressio::mpl::remove_cvref_t<MatrixType>,
      1
    >
{
  // type of the native expression
  using _native_expr_type = decltype(std::declval<MatrixType>().diagonal());
  using _const_native_expr_type=decltype(std::declval<const MatrixType>().diagonal());

  using native_expr_type = typename std::conditional<
    std::is_const<MatrixType>::value,
    _const_native_expr_type,
    _native_expr_type
  >::type;

  using reference_type = decltype( std::declval<_native_expr_type>()(0) );
  using const_reference_type = decltype( std::declval<_const_native_expr_type>()(0) );
};
#endif

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename MatrixType>
struct DiagTraits<
  DiagExpr<MatrixType>,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_dense_matrix_kokkos<MatrixType>::value
    >
  >
  : public ::pressio::impl::KokkosTraits<
      typename ::pressio::mpl::remove_cvref_t<MatrixType>,
      1,
      true
    >
{
  using native_expr_type = Kokkos::View<
    typename ::pressio::mpl::remove_cvref_t<MatrixType>::traits::value_type*,
    Kokkos::LayoutStride
  >;

  // using _native_expr_t	     = Kokkos::View<scalar_t*, Kokkos::LayoutStride>;
  // using _const_native_expr_t = Kokkos::View<const scalar_t*, Kokkos::LayoutStride>;
  // using native_expr_t = typename std::conditional<
  //   std::is_const<MatrixType>::value,
  //   _const_native_expr_t,
  //   _native_expr_t
  // >::type;
};
#endif

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <typename MatrixType>
struct DiagTraits<
  DiagExpr<MatrixType>,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_array_pybind<MatrixType>::value
    >
  >
  : public ::pressio::impl::ContainersSharedTraits<PackageIdentifier::Pybind, true, 1>
{
  using mat_remove_cv_t = typename std::remove_cv<MatrixType>::type;
  using scalar_type  = typename ::pressio::Traits<mat_remove_cv_t>::scalar_type;
  using size_type    = typename ::pressio::Traits<mat_remove_cv_t>::ordinal_type;
  using reference_type =  scalar_type &;
  using const_reference_type = scalar_type const &;
};
#endif

}}}
#endif  // EXPRESSIONS_IMPL_DIAG_TRAITS_HPP_
