/*
//@HEADER
// ************************************************************************
//
// solvers_has_all_needed_methods_for_hessian.hpp
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

#ifndef SOLVERS_SRC_META_SOLVERS_HAS_ALL_NEEDED_METHODS_FOR_HESSIAN_HPP_
#define SOLVERS_SRC_META_SOLVERS_HAS_ALL_NEEDED_METHODS_FOR_HESSIAN_HPP_

namespace pressio{ namespace solvers{ namespace meta {

template <
  typename T, typename jac_t, typename result_t, typename scalar_t,
  typename enable = void
  >
struct has_all_needed_methods_for_hessian
  : std::false_type{};

template <typename T, typename jac_t, typename result_t, typename scalar_t>
struct has_all_needed_methods_for_hessian<
  T, jac_t, result_t, scalar_t,
  mpl::enable_if_t<
    /*::pressio::ops::meta::has_nonvoid_method_product_mat_mat<
      T, ::pressio::transpose, ::pressio::nontranspose, scalar_t, jac_t, jac_t, result_t
      >::value and*/
    ::pressio::ops::meta::has_void_method_product_mat_mat<
      T, ::pressio::transpose, ::pressio::nontranspose, scalar_t, jac_t, jac_t, result_t
      >::value
    >
  > : std::true_type{};

}}} //pressio::solvers::meta
#endif
