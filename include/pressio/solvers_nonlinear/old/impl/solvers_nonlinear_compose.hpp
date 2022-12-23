/*
//@HEADER
// ************************************************************************
//
// solvers_nonlinear_compose.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_SOLVERS_NONLINEAR_COMPOSE_HPP_
#define SOLVERS_NONLINEAR_IMPL_SOLVERS_NONLINEAR_COMPOSE_HPP_

#include "./operators/solvers_weighting_irwls.hpp"
#include "./operators/solvers_gn_hessian_gradient_operators_hg_api.hpp"
#include "./operators/solvers_gn_hessian_gradient_operators_rj_api.hpp"
#include "./operators/solvers_gn_hessian_gradient_operators_weighted_rj_api.hpp"
#include "./operators/solvers_lm_hessian_gradient_operators_hg_api.hpp"
#include "./operators/solvers_lm_hessian_gradient_operators_rj_api.hpp"
#include "./operators/solvers_residual_jacobian_operators.hpp"

#include "./correctors/solvers_hessian_gradient_corrector.hpp"
#include "./correctors/solvers_qr_corrector.hpp"
#include "./correctors/solvers_rj_corrector.hpp"
#include "solver.hpp"

namespace pressio{ namespace nonlinearsolvers{ namespace impl{

////////////////////////////////////////////////////////////
//
// *** GN or LM ***
//
////////////////////////////////////////////////////////////

template< class SystemType, class tag, class ... Args>
struct Compose{ using type = void; };

// ***
// GN or LM with Neq and HG api
// ***
template<class SystemType, class tag, class LinearSolverType>
#ifdef PRESSIO_ENABLE_CXX20
requires (SystemWithHessianAndGradient<mpl::remove_cvref_t<SystemType>>
	  || SystemWithFusedHessianAndGradient<mpl::remove_cvref_t<SystemType>>)
#endif
struct Compose<
  SystemType, tag,
  mpl::enable_if_t<
#if not defined PRESSIO_ENABLE_CXX20
    (pressio::nonlinearsolvers::SystemWithHessianAndGradient<mpl::remove_cvref_t<SystemType>>::value or
     pressio::nonlinearsolvers::SystemWithFusedHessianAndGradient<mpl::remove_cvref_t<SystemType>>::value)
    and
#endif
    (std::is_same<tag, GaussNewton>::value or std::is_same<tag, LM>::value)
    >,
  LinearSolverType
  >
{
  using system_type = mpl::remove_cvref_t<SystemType>;
  using linear_solver_type = mpl::remove_cvref_t<LinearSolverType>;

  using residual_norm_type = typename system_type::residual_norm_type;
  using state_t  = typename system_type::state_type;
  using grad_t   = typename system_type::gradient_type;

  using operators_t =
    typename std::conditional<
    std::is_same<tag, GaussNewton>::value,
    HessianGradientOperatorsHGApi<typename system_type::hessian_type, grad_t, residual_norm_type>,
    LMHessianGradientOperatorsHGApi<typename system_type::hessian_type, grad_t, residual_norm_type>
    >::type;

  using corrector_t = HessianGradientCorrector<operators_t, state_t, LinearSolverType>;
  using type = Solver<tag, corrector_t>;
};

// ------------------------------------------------------------
// GN or LM with Neq and r/j api
// ------------------------------------------------------------
template<class SystemType, class tag, class LinearSolverType>
#ifdef PRESSIO_ENABLE_CXX20
requires (SystemWithResidualAndJacobian<mpl::remove_cvref_t<SystemType>>
       || SystemWithFusedResidualAndJacobian<mpl::remove_cvref_t<SystemType>>)
#endif
struct Compose<
  SystemType, tag,
  mpl::enable_if_t<
#if not defined PRESSIO_ENABLE_CXX20
    (::pressio::nonlinearsolvers::SystemWithResidualAndJacobian<mpl::remove_cvref_t<SystemType>>::value or
     ::pressio::nonlinearsolvers::SystemWithFusedResidualAndJacobian<mpl::remove_cvref_t<SystemType>>::value)
    and
#endif
    (std::is_same<tag, GaussNewton>::value or std::is_same<tag, LM>::value)
    >, LinearSolverType
  >
{
  using system_type = mpl::remove_cvref_t<SystemType>;
  using state_t  = typename system_type::state_type;
  using r_t	 = typename system_type::residual_type;
  using j_t	 = typename system_type::jacobian_type;

  // hessian_t is extracted from linear solver
  using linear_solver_type = mpl::remove_cvref_t<LinearSolverType>;
  using hess_t = typename mpl::remove_cvref_t<linear_solver_type>::matrix_type;
  // gradient is same type as state
  using grad_t = state_t;

  using operators_t =
    typename std::conditional<
    std::is_same<tag, GaussNewton>::value,
    HessianGradientOperatorsRJApiNoWeighting<hess_t, grad_t, r_t, j_t>,
    LMHessianGradientOperatorsRJApi<
      hess_t, grad_t, r_t, j_t, HessianGradientOperatorsRJApiNoWeighting
      >
    >::type;

  using corrector_t = HessianGradientCorrector<operators_t, state_t, LinearSolverType>;
  using type = Solver<tag, corrector_t>;
};

template<class weighting_functor_t, class tag>
struct ComposeOperators;

template<>
struct ComposeOperators<void, GaussNewton>{
  template<class ...Args> using type = HessianGradientOperatorsRJApiNoWeighting<Args...>;
};

template<class weighting_functor_t>
struct ComposeOperators<weighting_functor_t, GaussNewton>{
  template<class ...Args> using type = WeightedHessianGradientOperatorsRJApi<Args..., weighting_functor_t>;
};

template<>
struct ComposeOperators<void, LM>{
  template<class ...Args>
  using type = LMHessianGradientOperatorsRJApi<Args..., HessianGradientOperatorsRJApiNoWeighting>;
};

template<class weighting_functor_t>
struct ComposeOperators<weighting_functor_t, LM>{
  template<class ...Args>
  using type = LMHessianGradientOperatorsRJApi<Args...,
					       WeightedHessianGradientOperatorsRJApi,
					       weighting_functor_t>;
};

// ---------------------------------------------------------------------------
// *** weighted GN or LM with Neq and r/j api ***
// ---------------------------------------------------------------------------
template<class SystemType, class tag, class LinearSolverType, class WeightingOperator>
struct Compose<
  SystemType, tag,
  mpl::enable_if_t<std::is_same<tag, GaussNewton>::value or std::is_same<tag, LM>::value>,
  LinearSolverType, WeightingOperator
  >
{
  using system_type = mpl::remove_cvref_t<SystemType>;
  using linear_solver_type = mpl::remove_cvref_t<LinearSolverType>;

  using state_t = typename system_type::state_type;
  using r_t = typename system_type::residual_type;
  using j_t = typename system_type::jacobian_type;

  // hessian_t is extracted from linear solver
  using hess_t = typename linear_solver_type::matrix_type;
  // gradient is same type as state
  using grad_t = state_t;

  using operators_t = typename ComposeOperators<WeightingOperator, tag>::template type<
    hess_t, grad_t, r_t, j_t>;

  using corrector_t = HessianGradientCorrector<operators_t, state_t, LinearSolverType>;
  using type = Solver<tag, corrector_t>;
};

// ------------------------------------------------------------
// COMPOSE IRW GN with Neq and r/j api, M=I
// ------------------------------------------------------------
template<class SystemType, class LinearSolverType>
struct ComposeIrwGaussNewton
{
  using system_type = mpl::remove_cvref_t<SystemType>;
  using state_t = typename system_type::state_type;
  using r_t = typename system_type::residual_type;
  using j_t = typename system_type::jacobian_type;
  using grad_t = state_t;

  using weighting_t = ::pressio::nonlinearsolvers::impl::IrwWeightingOperator<r_t, j_t>;
  using type = typename Compose<SystemType, GaussNewton,
				void, LinearSolverType, weighting_t>::type;
};

// ------------------------------------------------------------

template<class SystemType, class ... Args>
using ComposeGaussNewton_t = typename Compose<SystemType, GaussNewton, void, Args...>::type;

template<class SystemType, class ... Args>
using ComposeLM_t = typename Compose<SystemType, LM, void, Args...>::type;

template<class SystemType, class ... Args>
using ComposeLevenbergMarquardt_t = ComposeLM_t<SystemType, Args...>;

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_SOLVERS_NONLINEAR_COMPOSE_HPP_
