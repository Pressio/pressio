/*
//@HEADER
// ************************************************************************
//
// containers_span_function.hpp
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

#ifndef CONTAINERS_EXPRESSIONS_PUBLIC_FUNCTIONS_HPP_
#define CONTAINERS_EXPRESSIONS_PUBLIC_FUNCTIONS_HPP_

#include "impl/span_traits.hpp"
#include "impl/span_classes.hpp"
#include "impl/subspan_traits.hpp"
#include "impl/subspan_classes.hpp"
#include "impl/diag_traits.hpp"
#include "impl/diag_classes.hpp"
#include "impl/asdiagonalmatrix_traits.hpp"
#include "impl/asdiagonalmatrix_classes.hpp"

namespace pressio{ 

//----------------
// span
//----------------
template <typename T, typename ... Args>
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
  // for kokkos, we don't need this overload 
 mpl::enable_if_t<!::pressio::is_vector_kokkos<T>::value, expressions::impl::SpanExpr<T>>
#else
 expressions::impl::SpanExpr<T> 
#endif
span(T & vecObj, Args&& ... args)
{
  static_assert(0 < sizeof...(Args), 
    "span must be called with arguments specifying the span bounds.");
  static_assert(::pressio::Traits<T>::rank==1, 
    "span can only be applied to a rank-1 object.");

  using return_t = expressions::impl::SpanExpr<T>;
  return return_t(vecObj, std::forward<Args>(args)... );
}

template <typename T, typename ... Args>
expressions::impl::SpanExpr<const T> span(const T & vecObj, Args&& ... args)
{  
  static_assert(0 < sizeof...(Args), 
    "span must be called with arguments specifying the span bounds.");
  static_assert(::pressio::Traits<T>::rank==1, 
    "span can only be applied to a rank-1 object.");

  using return_t = expressions::impl::SpanExpr<const T>;
  return return_t(vecObj, std::forward<Args>(args)... );
}

//----------------
// subspan
//----------------
template <typename T, typename ... Args>
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
  // for kokkos, we don't need this overload 
 mpl::enable_if_t<!::pressio::is_dense_matrix_kokkos<T>::value, expressions::impl::SubspanExpr<T>>
#else
 expressions::impl::SubspanExpr<T> 
#endif
subspan(T & obj, Args&& ... args)
{
  static_assert(0 < sizeof...(Args), 
    "subspan must be called with arguments specifying the bounds.");
  static_assert(::pressio::Traits<T>::rank==2, 
    "subspan can only be applied to a rank-2 object.");

  using return_t = expressions::impl::SubspanExpr<T>;
  return return_t(obj, std::forward<Args>(args)... );
}

template <typename T, typename ... Args>
expressions::impl::SubspanExpr<const T> subspan(const T & obj, Args&& ... args)
{
  static_assert(0 < sizeof...(Args), 
    "subspan must be called with arguments specifying the bounds.");
  static_assert(::pressio::Traits<T>::rank==2, 
    "subspan can only be applied to a rank-2 object.");

  using return_t = expressions::impl::SubspanExpr<const T>;
  return return_t(obj, std::forward<Args>(args)... );
}

//----------------
// diag
//----------------
template <typename T, typename ... Args>
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
  // for kokkos, we don't need this overload 
 mpl::enable_if_t<!::pressio::is_dense_matrix_kokkos<T>::value, expressions::impl::DiagExpr<T>>
#else
 expressions::impl::DiagExpr<T> 
#endif
diag(T & obj)
{
  static_assert(::pressio::Traits<T>::rank==2, 
    "diag can only be applied to a rank-2 object.");

  using return_t = expressions::impl::DiagExpr<T>;
  return return_t(obj);
}

template <typename T, typename ... Args>
expressions::impl::DiagExpr<const T> diag(const T & obj)
{
  static_assert(::pressio::Traits<T>::rank==2, 
    "diag can only be applied to a rank-2 object.");

  using return_t = expressions::impl::DiagExpr<const T>;
  return return_t(obj);
}

//------------------
// asDiagonalMatrix
//------------------
template <typename T>
expressions::impl::AsDiagonalMatrixExpr<T> asDiagonalMatrix(T & vecObj)
{
  static_assert(::pressio::Traits<T>::rank==1, 
    "AsDiagonalMatrix can only be applied to a rank-1 object.");
  using return_t = expressions::impl::AsDiagonalMatrixExpr<T>;
  return return_t(vecObj);
}

template <typename T>
expressions::impl::AsDiagonalMatrixExpr<const T> asDiagonalMatrix(const T & vecObj)
{
  static_assert(::pressio::Traits<T>::rank==1, 
    "AsDiagonalMatrix can only be applied to a rank-1 object.");
  using return_t = expressions::impl::AsDiagonalMatrixExpr<const T>;
  return return_t(vecObj);
}

}
#endif  // CONTAINERS_EXPRESSIONS_SPAN_CONTAINERS_SPAN_FUNCTION_HPP_
