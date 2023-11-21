/*
//@HEADER
// ************************************************************************
//
// public_functions.hpp
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

#ifndef EXPRESSIONS_COLUMN_HPP_
#define EXPRESSIONS_COLUMN_HPP_

#include "impl/column_traits.hpp"
#include "impl/column_classes.hpp"

namespace pressio{

template <class T, class IndexType>
auto column(T & operand, IndexType colIndex)
{
  // note that this works also when T is const-qualified
  // because that qualification carries over to the impl

  constexpr bool constraint = false
#ifdef PRESSIO_ENABLE_TPL_EIGEN
    || is_dense_matrix_eigen<T>::value
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
    || is_multi_vector_tpetra<T>::value
#endif
    ;
  static_assert(constraint, "pressio::column() currently supported only for "
		"an Eigen dynamic matrix or a tpetra multivector");

  static_assert(Traits< std::remove_const_t<T> >::rank==2,
		"column can only be applied to a rank-2 object.");
  static_assert(std::is_integral_v<IndexType>);

  return expressions::impl::ColumnExpr<T>(operand, colIndex);
}

}
#endif  // EXPRESSIONS_SPAN_HPP_
