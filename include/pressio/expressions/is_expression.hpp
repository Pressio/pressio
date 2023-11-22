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

#define PRESSIO_IMPL_IS_EXPRESSION(NAMEA, NAMEB) \
  template <typename T> struct is_expression_##NAMEA : std::false_type{}; \
  template <typename T> struct is_expression_##NAMEA<			\
    ::pressio::expressions::impl::NAMEB##Expr<T> > : std::true_type{};	\
  template <typename T> struct is_expression_##NAMEA<			\
    const ::pressio::expressions::impl::NAMEB##Expr<T> > : std::true_type{};	\
  template <typename T> struct is_expression_##NAMEA<			\
    ::pressio::expressions::impl::NAMEB##Expr<const T> > : std::true_type{}; \
  template <typename T> struct is_expression_##NAMEA<			\
   const ::pressio::expressions::impl::NAMEB##Expr<const T> > : std::true_type{}; \

PRESSIO_IMPL_IS_EXPRESSION(diagonal, Diagonal)
PRESSIO_IMPL_IS_EXPRESSION(span, Span)
PRESSIO_IMPL_IS_EXPRESSION(subspan, Subspan)
PRESSIO_IMPL_IS_EXPRESSION(column, Column)

/* detect any expression */
template <typename T, typename enable = void>
struct is_expression : std::false_type{};

template <typename T>
struct is_expression<
  T,
  mpl::enable_if_t<
    is_expression_span<T>::value
    || is_expression_diagonal<T>::value
    || is_expression_subspan<T>::value
    || is_expression_column<T>::value
    >
  > : std::true_type{};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <typename T>
struct is_expression_acting_on_eigen: public std::false_type {};

template <typename T>
struct is_expression_acting_on_eigen<
  ::pressio::expressions::impl::DiagonalExpr<T>
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

template <typename T>
struct is_expression_acting_on_eigen<
  ::pressio::expressions::impl::ColumnExpr<T>
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
  ::pressio::expressions::impl::DiagonalExpr<T>
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

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template <typename T>
struct is_expression_acting_on_tpetra: public std::false_type {};

template <typename T>
struct is_expression_acting_on_tpetra<
  ::pressio::expressions::impl::ColumnExpr<T>
  >
{
  static constexpr auto value = ::pressio::is_multi_vector_tpetra<T>::value;
};

template <typename T>
struct is_expression_column_acting_on_tpetra{
  static constexpr auto value = ::pressio::is_expression_column<T>::value
    && is_expression_acting_on_tpetra<T>::value;
};

template <typename T>
struct is_expression_acting_on_tpetra_block: public std::false_type {};

template <typename T>
struct is_expression_acting_on_tpetra_block<
  ::pressio::expressions::impl::ColumnExpr<T>
  >
{
  static constexpr auto value = ::pressio::is_multi_vector_tpetra_block<T>::value;
};

template <typename T>
struct is_expression_column_acting_on_tpetra_block{
  static constexpr auto value = ::pressio::is_expression_column<T>::value
    && is_expression_acting_on_tpetra_block<T>::value;
};
#endif

}
#endif  // EXPRESSIONS_IS_EXPRESSION_HPP_
