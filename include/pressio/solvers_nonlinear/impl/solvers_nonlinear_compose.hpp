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
////////////////////////////////////////////////////////////
//
// *** NewtonRaphson ***
//
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

template<typename SystemType, typename LinearSolverType>
struct ComposeNewtonRaphson
{
  // Newton-Raphson requires r/j API
  static_assert
  (::pressio::nonlinearsolvers::compliant_with_residual_jacobian_api<SystemType>::value or
   ::pressio::nonlinearsolvers::compliant_with_fused_residual_jacobian_api<SystemType>::value,
   "Newton-Raphson: a system with residual/jacobian API is required.");

  using scalar_t = typename SystemType::scalar_type;
  using state_t  = typename SystemType::state_type;
  using r_t = typename SystemType::residual_type;
  using j_t = typename SystemType::jacobian_type;

  static_assert
  (::pressio::nonlinearsolvers::admissible_state<state_t>::value,
   "Newton-Raphson: invalid state type");

  // check the solver_t passed is valid
  static_assert
  (::pressio::nonlinearsolvers::admissible_linear_solver_for_newton_raphson<
   mpl::remove_cvref_t<LinearSolverType>, state_t>::value,
   "Newton-Raphson: invalid linear solver");

  using operators_t = ResidualJacobianOperators<r_t, j_t, scalar_t>;
  using corrector_t = RJCorrector<operators_t, state_t, LinearSolverType>;
  using type = Solver<NewtonRaphson, corrector_t>;
};

template<typename SystemType, typename ... Args>
using ComposeNewtonRaphson_t = typename ComposeNewtonRaphson<SystemType, Args...>::type;


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//
// *** COMPOSE GN QR ***
//
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

template<typename SystemType, typename SolverType, class ...Args>
struct ComposeGNQR
{
  static_assert
  (::pressio::nonlinearsolvers::compliant_with_residual_jacobian_api<SystemType>::value or
   ::pressio::nonlinearsolvers::compliant_with_fused_residual_jacobian_api<SystemType>::value,
   "QR-based GaussNewton requires the residual/jacobian API.");

  // if we get here, system meets the r/j api, need to check if the
  // SolverType passed is valid for QR
  using scalar_t = typename SystemType::scalar_type;
  using state_t = typename SystemType::state_type;
  using r_t = typename SystemType::residual_type;
  using j_t = typename SystemType::jacobian_type;

  static_assert
  (::pressio::nonlinearsolvers::admissible_state<state_t>::value,
   "invalid state type");

  static_assert
  (::pressio::nonlinearsolvers::admissible_qr_solver_for_gn_qr<
   mpl::remove_cvref_t<SolverType>, state_t, j_t, r_t>::value,
   "The solver type passed to compose a QR-based GN solver is not valid");

  using operators_t = ResidualJacobianOperators<r_t, j_t, scalar_t>;
  using corrector_t  = QRCorrector<operators_t, state_t, SolverType>;
  using type        = Solver<GaussNewtonQR, corrector_t>;
};

template<class ...Args>
using ComposeGaussNewtonQR_t = typename ComposeGNQR<Args...>::type;


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//
// *** COMPOSE for GN or LM ***
//
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

template<typename SystemType, typename tag, typename ... Args>
struct Compose{
  using type = void;
};


// ***
// COMPOSE GN or LM with Neq and HG api
// ***
template<
  typename SystemType,
  typename tag,
  typename LinearSolverType
  >
struct Compose<
  SystemType, tag,
  mpl::enable_if_t<
    (pressio::nonlinearsolvers::compliant_with_hessian_gradient_api<SystemType>::value or
     pressio::nonlinearsolvers::compliant_with_fused_hessian_gradient_api<SystemType>::value) and
    (std::is_same<tag, GaussNewton>::value or std::is_same<tag, LM>::value)
    >,
  LinearSolverType
  >
{
  using scalar_t = typename SystemType::scalar_type;
  using state_t = typename SystemType::state_type;
  using grad_t = state_t;

  static_assert
  (::pressio::nonlinearsolvers::admissible_state<state_t>::value,
   "invalid state type");

  static_assert
  (::pressio::nonlinearsolvers::admissible_linear_solver_for_nonlinear_least_squares<
   mpl::remove_cvref_t<LinearSolverType>, state_t>::value,
   "A valid linear solver type must be passed to GN or LM with normal equations");

  // hessian_t is extracted from linear solver
  using hess_t = typename mpl::remove_cvref_t<LinearSolverType>::matrix_type;

  using operators_t =
    typename std::conditional<
    std::is_same<tag, GaussNewton>::value,
    HessianGradientOperatorsHGApi<hess_t,   grad_t, scalar_t>,
    LMHessianGradientOperatorsHGApi<hess_t, grad_t, scalar_t>
    >::type;

  using corrector_t = HessianGradientCorrector<operators_t, state_t, LinearSolverType>;
  using type = Solver<tag, corrector_t>;
};


// ------------------------------------------------------------
// COMPOSE GN or LM with Neq and r/j api, M=I
// ------------------------------------------------------------
template<
  typename SystemType,
  typename tag,
  typename LinearSolverType
  >
struct Compose<
  SystemType, tag,
  mpl::enable_if_t<
    (::pressio::nonlinearsolvers::compliant_with_residual_jacobian_api<SystemType>::value or
     ::pressio::nonlinearsolvers::compliant_with_fused_residual_jacobian_api<SystemType>::value)
    and
    (std::is_same<tag, GaussNewton>::value or std::is_same<tag, LM>::value)
    >,
  LinearSolverType
  >
{
  //if we get here, SystemType meets r/j api
  using scalar_t = typename SystemType::scalar_type;
  using state_t = typename SystemType::state_type;
  using r_t = typename SystemType::residual_type;
  using j_t = typename SystemType::jacobian_type;
  using grad_t = state_t;

  static_assert
  (::pressio::nonlinearsolvers::admissible_state<state_t>::value,
   "invalid state type");

  static_assert
  (::pressio::nonlinearsolvers::admissible_linear_solver_for_nonlinear_least_squares<
   mpl::remove_cvref_t<LinearSolverType>, state_t>::value,
   "A valid linear solver type must be passed to GN or LM with normal equations");

  // hessian_t is extracted from linear solver
  using hess_t = typename mpl::remove_cvref_t<LinearSolverType>::matrix_type;

  using operators_t =
    typename std::conditional<
    std::is_same<tag, GaussNewton>::value,
    HessianGradientOperatorsRJApiNoWeighting<hess_t, grad_t, r_t, j_t, scalar_t>,
    LMHessianGradientOperatorsRJApi<
      hess_t, grad_t, r_t, j_t, scalar_t, HessianGradientOperatorsRJApiNoWeighting
      >
    >::type;

  using corrector_t = HessianGradientCorrector<operators_t, state_t, LinearSolverType>;
  using type = Solver<tag, corrector_t>;
};

// ---------------------------------------------------------------------------
// *** COMPOSE GN or LM with Neq and r/j api and valid M  ***
// ---------------------------------------------------------------------------
template< typename weighting_functor_t, typename tag>
struct ComposeOperators;

template<>
struct ComposeOperators<void, GaussNewton>
{
  template<typename ...Args>
  using type = HessianGradientOperatorsRJApiNoWeighting<Args...>;
};

template<typename weighting_functor_t>
struct ComposeOperators<weighting_functor_t, GaussNewton>
{
  template<typename ...Args>
  using type = WeightedHessianGradientOperatorsRJApi<Args..., weighting_functor_t>;
};

template<>
struct ComposeOperators<void, LM>
{
  template<typename ...Args>
  using type = LMHessianGradientOperatorsRJApi<
    Args..., HessianGradientOperatorsRJApiNoWeighting>;
};

template<typename weighting_functor_t>
struct ComposeOperators<weighting_functor_t, LM>
{
  template<typename ...Args>
  using type = LMHessianGradientOperatorsRJApi<
    Args..., WeightedHessianGradientOperatorsRJApi, weighting_functor_t>;
};


// ------------------------
template<
  typename SystemType,
  typename tag,
  typename LinearSolverType,
  typename extra_t
  >
struct Compose<
  SystemType, tag,
  mpl::enable_if_t<
    (::pressio::nonlinearsolvers::compliant_with_residual_jacobian_api<SystemType>::value or
     ::pressio::nonlinearsolvers::compliant_with_fused_residual_jacobian_api<SystemType>::value) and
    (std::is_same<tag, GaussNewton>::value or std::is_same<tag, LM>::value)
    >,
  LinearSolverType, extra_t
  >
{
  //if we get here, SystemType meets r/j api
  using scalar_t = typename SystemType::scalar_type;
  using state_t = typename SystemType::state_type;
  using r_t = typename SystemType::residual_type;
  using j_t = typename SystemType::jacobian_type;
  using grad_t = state_t;

  static_assert
  (::pressio::nonlinearsolvers::admissible_state<state_t>::value,
   "invalid state type");

  static_assert
  (::pressio::nonlinearsolvers::admissible_linear_solver_for_nonlinear_least_squares<
   mpl::remove_cvref_t<LinearSolverType>, state_t>::value,
   "A valid linear solver type must be passed to GN or LM with normal equations");

  // hessian_t is extracted from linear solver
  using hess_t = typename mpl::remove_cvref_t<LinearSolverType>::matrix_type;

  // if we get here, it means that:
  // - system meets r/j api
  // - we need to compose GN or LM with normal equations
  // - we have one optional template paramter, which can
  //   either be a udops or the weighting

  using extra_nocvref_t = mpl::remove_cvref_t<extra_t>;

  // using ud_ops_t  =
  //   typename std::conditional<
  //   ::pressio::nonlinearsolvers::constraints::ops_normal_equations_rj_api<
  //     extra_nocvref_t, scalar_t, hess_t, grad_t, j_t, r_t>::value,
  //   extra_t, void >::type;

  using weighting_functor_t =
    typename std::conditional<
    (::pressio::nonlinearsolvers::admissible_least_squares_weighting_operator<extra_nocvref_t, r_t, j_t>::value),
    extra_t, void >::type;

  // // the extra_t should be one of the two, cannot be both
  // static_assert
  // ( mpl::not_void<weighting_functor_t>::value,
  //   "Both weighting_functor_t are non-void, something is wrong.");

  using operators_t =
    typename ComposeOperators<weighting_functor_t, tag>::template type<
    hess_t, grad_t, r_t, j_t, scalar_t>;

  using corrector_t = HessianGradientCorrector<operators_t, state_t, LinearSolverType>;
  using type = Solver<tag, corrector_t>;
};

// ------------------------------------------------------------
// COMPOSE IRW GN with Neq and r/j api, M=I
// ------------------------------------------------------------
template<
  typename SystemType,
  typename LinearSolverType
  >
struct Compose<
  SystemType, IrwGaussNewton,
  mpl::enable_if_t<
    ::pressio::nonlinearsolvers::compliant_with_residual_jacobian_api<SystemType>::value or
    ::pressio::nonlinearsolvers::compliant_with_fused_residual_jacobian_api<SystemType>::value
    >,
  LinearSolverType
  >
{
  //if we get here, SystemType meets r/j api
  using scalar_t = typename SystemType::scalar_type;
  using state_t = typename SystemType::state_type;
  using r_t = typename SystemType::residual_type;
  using j_t = typename SystemType::jacobian_type;
  using grad_t = state_t;

  using weighting_t = ::pressio::nonlinearsolvers::impl::IrwWeightingOperator<r_t, j_t, scalar_t>;
  using type = typename Compose<SystemType, GaussNewton, void, LinearSolverType, weighting_t>::type;
};

// ------------------------------------------------------------
template<typename SystemType, typename ... Args>
using ComposeGaussNewton_t = typename Compose<SystemType, GaussNewton, void, Args...>::type;

template<typename SystemType, typename ... Args>
using ComposeLM_t = typename Compose<SystemType, LM, void, Args...>::type;

template<typename SystemType, typename ... Args>
using ComposeLevenbergMarquardt_t = ComposeLM_t<SystemType, Args...>;

template<typename SystemType, typename lin_s_t>
using ComposeIrwGaussNewton = Compose<SystemType, IrwGaussNewton, void, lin_s_t>;

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_SOLVERS_NONLINEAR_COMPOSE_HPP_
