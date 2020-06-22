/*
//@HEADER
// ************************************************************************
//
// solvers_compose.hpp
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

#ifndef PRESSIO_SOLVERS_COMPOSE_HPP_
#define PRESSIO_SOLVERS_COMPOSE_HPP_

#include "./operators/solvers_hessian_gradient_operators.hpp"
#include "./correction_mixin/solvers_hessian_gradient_corrector.hpp"
#include "./correction_mixin/solvers_qr_corrector.hpp"

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

// ----------------------------------------------------------------------------
// *** COMPOSE CORRECTOR ***
// ----------------------------------------------------------------------------
template<typename tag, typename ... Args>
struct composeCorrector;

template<typename system_t, typename state_t, typename h_t, typename g_t, typename lin_solver_t>
struct composeCorrector<
  GaussNewton,
  mpl::enable_if_t<
    pressio::solvers::meta::system_meets_residual_jacobian_api<system_t>::value or
    pressio::solvers::meta::system_meets_fused_residual_jacobian_api<system_t>::value>,
  system_t, state_t, h_t, g_t, lin_solver_t
  >
{
  using r_t = typename system_t::residual_type;
  using j_t = typename system_t::jacobian_type;
  static constexpr auto norm = pressio::solvers::Norm::L2;
  using operators_t = HessianGradientOperatorsRJApi<h_t, g_t, r_t, j_t>;
  using type	    = HessianGradientCorrector<operators_t, state_t, lin_solver_t, norm>;
};


template<typename system_t, typename state_t, typename h_t, typename g_t, typename lin_solver_t>
struct composeCorrector<
  GaussNewton,
  mpl::enable_if_t<
    pressio::solvers::meta::system_meets_hessian_gradient_api<system_t>::value or
    pressio::solvers::meta::system_meets_fused_hessian_gradient_api<system_t>::value>,
  system_t, state_t, h_t, g_t, lin_solver_t
  >
{
  static constexpr auto norm = pressio::solvers::Norm::L2;
  using operators_t = HessianGradientOperatorsHGApi<h_t, g_t>;
  using type	    = HessianGradientCorrector<operators_t, state_t, lin_solver_t, norm>;
};


template<typename system_t, typename state_t, typename h_t, typename g_t, typename lin_solver_t>
struct composeCorrector<
  LM,
  mpl::enable_if_t<
    pressio::solvers::meta::system_meets_residual_jacobian_api<system_t>::value or
    pressio::solvers::meta::system_meets_fused_residual_jacobian_api<system_t>::value>,
  system_t, state_t, h_t, g_t, lin_solver_t
  >
{
  using r_t = typename system_t::residual_type;
  using j_t = typename system_t::jacobian_type;
  static constexpr auto norm = pressio::solvers::Norm::L2;
  using operators_t = LMHessianGradientOperatorsRJApi<h_t, g_t, r_t, j_t>;
  using type	    = HessianGradientCorrector<operators_t, state_t, lin_solver_t, norm>;
};

// ----------------------------------------------------------------------------
// *** COMPOSE SOLVER TYPE ***
// ----------------------------------------------------------------------------
template<
  typename system_t,
  typename tag,
  template< typename...> class update,
  template< typename...> class looper,
  typename ... Args>
struct compose{
  using type = void;
};

template<
  typename system_t,
  typename tag,
  template<typename...> class update_t,
  template<typename...> class looper_t,
  typename ... Args
  >
struct compose<
  system_t, tag, update_t, looper_t,
  mpl::enable_if_t<
    std::is_same<tag, GaussNewton>::value or std::is_same<tag, LM>::value
    >, Args...
  >
{
  // GN or LM with neq need a valid linear solver for hess/grad system
  using ic2 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::solvers::meta::is_legitimate_linear_solver_for_least_squares_solver, Args...>;
  using linear_solver_t = ::pressio::mpl::variadic::at_or_t<void, ic2::value, Args...>;
  static_assert(!std::is_void<linear_solver_t>::value and ic2::value < sizeof... (Args),
  		"A valid linear solver type must be passed to GN with normal equations");

  using scalar_t = typename system_t::scalar_type;
  using state_t = typename system_t::state_type;
  // gradient is same as state_t
  using grad_t = state_t; // the gradient is same type as state
  // hessian_t is extracted from linear solver
  using hess_t = typename linear_solver_t::matrix_type;

  using corr_mixin = typename composeCorrector<
    tag, void, system_t, state_t, hess_t, grad_t, linear_solver_t>::type;

  // TODO: assert that the update is admissible for the tag
  using update_mixin  = update_t<scalar_t, state_t, corr_mixin>;

  using type = looper_t<scalar_t, update_mixin>;
};

}}}}
#endif



/* template< */
/*   typename system_t, */
/*   typename state_t, */
/*   template<typename...> class update, */
/*   template<typename...> class looper, */
/*   typename solver_t, */
/*   typename ... Args */
/*   > */
/* typename compose<system_t, LM, update, looper, void, solver_t, Args...>::type */
/* composer(const system_t & sys, const state_t & state, solver_t & solver, Args && ...args) */
/* { */
/*   using res_t = typename compose<system_t, LM, update, looper, void, solver_t, Args...>::type; */
/*   return res_t(sys, state, solver, std::forward<Args>(args)...); */
/* } */
