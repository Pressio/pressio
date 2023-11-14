/*
//@HEADER
// ************************************************************************
//
// is_expression.hpp
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

#ifndef EXPRESSIONS_IS_EXPRESSION_HPP_
#define EXPRESSIONS_IS_EXPRESSION_HPP_

namespace pressio{

/* is_expression_diag */
template <typename T>
struct is_expression_diag : std::false_type{};

template <typename T>
struct is_expression_diag<
  ::pressio::expressions::impl::DiagExpr<T>
  > : std::true_type{};

template <typename T>
struct is_expression_diag<
  const ::pressio::expressions::impl::DiagExpr<T>
  > : is_expression_diag<T>{};

template <typename T>
struct is_expression_diag<
  ::pressio::expressions::impl::DiagExpr<const T>
  > : std::true_type{};

template <typename T>
struct is_expression_diag<
  const ::pressio::expressions::impl::DiagExpr<const T>
  > : std::true_type{};


/* span */
template <typename T>
struct is_expression_span : std::false_type{};

template <typename T>
struct is_expression_span<
  ::pressio::expressions::impl::SpanExpr<T>
  > : std::true_type{};

template <typename T>
struct is_expression_span<
  const ::pressio::expressions::impl::SpanExpr<T>
  > : is_expression_span<T>{};

template <typename T>
struct is_expression_span<
  ::pressio::expressions::impl::SpanExpr<const T>
  > : std::true_type{};

template <typename T>
struct is_expression_span<
  const ::pressio::expressions::impl::SpanExpr<const T>
  > : std::true_type{};


/* subspan */
template <typename T>
struct is_expression_subspan : std::false_type{};

template <typename T>
struct is_expression_subspan<
  ::pressio::expressions::impl::SubspanExpr<T>
  > : std::true_type{};

template <typename T>
struct is_expression_subspan<
  const ::pressio::expressions::impl::SubspanExpr<T>
  > : is_expression_subspan<T>{};

template <typename T>
struct is_expression_subspan<
  ::pressio::expressions::impl::SubspanExpr<const T>
  > : std::true_type{};

template <typename T>
struct is_expression_subspan<
  const ::pressio::expressions::impl::SubspanExpr<const T>
  > : std::true_type{};


/* detect any expression */
template <typename T, typename enable = void>
struct is_expression : std::false_type{};

template <typename T>
struct is_expression<
  T,
  mpl::enable_if_t<
    is_expression_span<T>::value or
    is_expression_diag<T>::value or
    is_expression_subspan<T>::value
    >
  > : std::true_type{};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <typename T>
struct is_expression_acting_on_eigen: public std::false_type {};

template <typename T>
struct is_expression_acting_on_eigen<
  ::pressio::expressions::impl::DiagExpr<T>
  >
{
  static constexpr auto value = ::pressio::is_native_container_eigen<T>::value;
};

template <typename T>
struct is_expression_acting_on_eigen<
  ::pressio::expressions::impl::SpanExpr<T>
  >
{
  static constexpr auto value = ::pressio::is_native_container_eigen<T>::value;
};

template <typename T>
struct is_expression_acting_on_eigen<
  ::pressio::expressions::impl::SubspanExpr<T>
  >
{
  static constexpr auto value = ::pressio::is_native_container_eigen<T>::value;
};
#endif // PRESSIO_ENABLE_TPL_EIGEN


#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename T>
struct is_expression_acting_on_kokkos: public std::false_type {};

template <typename T>
struct is_expression_acting_on_kokkos<
  ::pressio::expressions::impl::DiagExpr<T>
  >
{
  static constexpr auto value = ::pressio::is_native_container_kokkos<T>::value;
};

template <typename T>
struct is_expression_acting_on_kokkos<
  ::pressio::expressions::impl::SpanExpr<T>
  >
{
  static constexpr auto value = ::pressio::is_native_container_kokkos<T>::value;
};

template <typename T>
struct is_expression_acting_on_kokkos<
  ::pressio::expressions::impl::SubspanExpr<T>
  >
{
  static constexpr auto value = ::pressio::is_native_container_kokkos<T>::value;
};
#endif // PRESSIO_ENABLE_TPL_KOKKOS

}
#endif  // EXPRESSIONS_IS_EXPRESSION_HPP_
