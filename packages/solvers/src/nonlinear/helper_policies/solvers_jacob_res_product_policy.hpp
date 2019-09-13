/*
//@HEADER
// ************************************************************************
//
// solvers_jacob_res_product_policy.hpp
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

#ifndef SOLVERS_IMPL_JACOBIAN_RESIDUAL_PRODUCT_POLICY_HPP
#define SOLVERS_IMPL_JACOBIAN_RESIDUAL_PRODUCT_POLICY_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../../../CONTAINERS_OPS"

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template<typename J_t, typename enable = void>
struct JacobianTranspResProdHelper;


// when J is a matrix wrapper, then J^T*R
// is computed doing regular mat-vec product
template<typename J_t>
struct JacobianTranspResProdHelper<
  J_t,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_matrix_wrapper<J_t>::value
    >>{

  template <typename resid_t, typename result_t>
  static void evaluate(J_t & J, resid_t & R, result_t & result)
  {
    constexpr bool transposeJ = true;
    ::pressio::containers::ops::product<J_t, resid_t, result_t,
				transposeJ>(J, R, result);
  }
};


// when J is multivector wrapper, then J^T*R
// can be computed doing the DOT of J*R
template<typename J_t>
struct JacobianTranspResProdHelper<
  J_t,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper<J_t>::value
    >>{

  template <typename resid_t, typename result_t>
  static void evaluate(J_t & J, resid_t & R, result_t & result) {
    ::pressio::containers::ops::dot(J, R, result);
  }
};


}}}} //end namespace pressio::solvers::iterative::impl
#endif
