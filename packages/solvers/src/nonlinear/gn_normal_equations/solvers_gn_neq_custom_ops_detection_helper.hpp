/*
//@HEADER
// ************************************************************************
//
// solvers_gn_custom_ops_detection_helper.hpp
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

#ifndef solvers_src_solvers_gn_neq_custom_ops_detection_helper_hpp_
#define solvers_src_solvers_gn_neq_custom_ops_detection_helper_hpp_

#include "../../meta/custom_ops_detection/solvers_has_all_needed_static_dot_self_overloads.hpp"
#include "../../meta/custom_ops_detection/solvers_has_all_needed_static_norm_methods.hpp"
#include "../../meta/custom_ops_detection/solvers_has_static_method_dot_three_args_return_void.hpp"

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template <
  typename gradient_t,
  typename residual_t,
  typename jacobian_t,
  typename hessian_t,
  typename ...Args
  >
struct GaussNewtonNEqCustomOpsDetectionHelper{

  using gradient_native_t = typename ::pressio::containers::details::traits<gradient_t>::wrapped_t;
  using residual_native_t = typename ::pressio::containers::details::traits<residual_t>::wrapped_t;
  using jacobian_native_t = typename ::pressio::containers::details::traits<jacobian_t>::wrapped_t;

  /*
   * detect if Args contain valid ops type with static methods to compute hessian
   */
  using icHessMV = ::pressio::mpl::variadic::find_if_ternary_pred_t<
    jacobian_native_t, hessian_t, ::pressio::solvers::meta::has_all_needed_static_dot_self_overloads, Args...>;
  using hess_op_t = ::pressio::mpl::variadic::at_or_t<void, icHessMV::value, Args...>;


  /*
   * detect if Args contain valid ops type with static methods to compute norm of residual
   */
  using resid_scalar_t	  = typename ::pressio::containers::details::traits<residual_t>::scalar_t;
  using icNorm		  = ::pressio::mpl::variadic::find_if_ternary_pred_t<
    residual_native_t, resid_scalar_t, ::pressio::solvers::meta::has_all_needed_static_norm_methods, Args...>;
  using norm_op_t	  = ::pressio::mpl::variadic::at_or_t<void, icNorm::value, Args...>;


  /*
   * detect if Args contain valid ops type with static methods to compute gradient J^T R
   */
  using icGrad = ::pressio::mpl::variadic::find_if_quaternary_pred_t<
    jacobian_native_t, residual_native_t, gradient_native_t,
    ::pressio::solvers::meta::has_static_method_dot_three_args_return_void, Args...>;
  using grad_op_t = ::pressio::mpl::variadic::at_or_t<void, icGrad::value, Args...>;


  //find if there is a single type that contains all methods for all ops
  static constexpr bool foundIt = std::is_same<hess_op_t, norm_op_t>::value and
				  std::is_same<hess_op_t, grad_op_t>::value and
				  !std::is_void<hess_op_t>::value;

  // here, I know that there is a single type for all ops, so it does not matter which one to use
  using ops_t = hess_op_t;

  using type = typename std::conditional< foundIt, ops_t, void >::type;

  static_assert( std::is_void<type>::value or 
    (!std::is_void<type>::value and ::pressio::containers::meta::is_multi_vector_wrapper<jacobian_t>::value),
     "For GaussNewton normal-eq solver, custom ops are currently supported when the \
jacobian is a multi-vector wrapper. If you are using this for doing ROM, this most likely\
means you wrapped the Jacobian type of your basis with a matrix not a multi-vector.");

};

}}}}//end namespace pressio::solvers::iterative::impl
#endif
