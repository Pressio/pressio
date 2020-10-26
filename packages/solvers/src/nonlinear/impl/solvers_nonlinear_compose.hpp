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

#include "./operators/solvers_gn_hessian_gradient_operators_hg_api.hpp"
#include "./operators/solvers_gn_hessian_gradient_operators_no_weighting.hpp"
#include "./operators/solvers_gn_hessian_gradient_operators_with_weighting.hpp"
#include "./operators/solvers_lm_hessian_gradient_operators.hpp"
#include "./operators/solvers_residual_jacobian_operators.hpp"

#include "./correction_mixins/solvers_hessian_gradient_corrector.hpp"
#include "./correction_mixins/solvers_qr_corrector.hpp"
#include "./correction_mixins/solvers_rj_corrector.hpp"
#include "solver.hpp"

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//
// *** NewtonRaphson ***
//
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

template<typename system_t, typename linear_solver_t>
struct composeNewRaph
{

  // Newton-Raphson requires r/j API
  static_assert
  (::pressio::solvers::concepts::system_residual_jacobian<system_t>::value or
   ::pressio::solvers::concepts::system_fused_residual_jacobian<system_t>::value,
   "To use NewtonRaphson, your system must meet the residual/jacobian API.");

  using scalar_t = typename system_t::scalar_type;
  using state_t  = typename system_t::state_type;
  using r_t = typename system_t::residual_type;
  using j_t = typename system_t::jacobian_type;

  static_assert
  (::pressio::containers::predicates::is_vector_wrapper<state_t>::value,
   "Newton-Raphson solver: the state type must be a pressio vector wrapper.");

  // check the solver_t passed is valid
  static_assert
  (::pressio::solvers::concepts::linear_solver_for_newton_raphson<
   linear_solver_t, state_t>::value,
   "Invalid linear solver type passed to NewtonRaphson");

  using operators_t = ResidualJacobianOperators<r_t, j_t>;
  using corr_mixin = RJCorrector<operators_t, state_t, linear_solver_t>;
  using type = Solver<NewtonRaphson, corr_mixin>;
};


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//
// *** COMPOSE GN QR ***
//
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

template<typename enable, typename system_t, typename ...Args>
struct composeGNQR
{
  static_assert
  (::pressio::solvers::concepts::system_residual_jacobian<system_t>::value or
   ::pressio::solvers::concepts::system_fused_residual_jacobian<system_t>::value,
   "QR-based GaussNewton requires the system with the residual/jacobian API.");
};

template<typename system_t, typename solver_t>
struct composeGNQR<
  mpl::enable_if_t<
    pressio::solvers::concepts::system_residual_jacobian<system_t>::value or
    pressio::solvers::concepts::system_fused_residual_jacobian<system_t>::value
    >, system_t, solver_t
  >
{
  // if we get here, system meets the r/j api, need to check if the
  // solver_t passed is valid for QR
  using scalar_t = typename system_t::scalar_type;
  using state_t = typename system_t::state_type;
  using r_t = typename system_t::residual_type;
  using j_t = typename system_t::jacobian_type;
  static_assert
  (::pressio::solvers::concepts::qr_solver_for_gn_qr<solver_t, state_t, j_t, r_t>::value,
   "The solver type passed to compose a QR-based GN solver is not valid");

  static_assert
  (::pressio::containers::predicates::is_vector_wrapper<state_t>::value,
   "Nonlinear least-squares solver: the state type must be a pressio vector wrapper.");

  // if we get here, all is ok
  using qr_solver_matrix_t =
    typename ::pressio::qr::details::traits<solver_t>::matrix_t;

  using operators_t = ResidualJacobianOperators<r_t, j_t>;
  using corr_mixin  = QRCorrector<operators_t, state_t, solver_t>;
  using type        = Solver<GaussNewtonQR, corr_mixin>;
};


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//
// *** COMPOSE for GN or LM ***
//
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

template<typename system_t, typename tag, typename ... Args>
struct compose{
  using type = void;
};


// ***
// COMPOSE GN or LM with Neq and HG api
// ***
template<
  typename system_t,
  typename tag,
  typename linear_solver_t
  >
struct compose<
  system_t, tag,
  mpl::enable_if_t<
    (pressio::solvers::concepts::system_hessian_gradient<system_t>::value or
     pressio::solvers::concepts::system_fused_hessian_gradient<system_t>::value) and
    (std::is_same<tag, GaussNewton>::value or std::is_same<tag, LM>::value)
    >,
  linear_solver_t
  >
{
  using scalar_t = typename system_t::scalar_type;
  using state_t = typename system_t::state_type;
  using grad_t = state_t;

  static_assert
  (::pressio::containers::predicates::is_vector_wrapper<state_t>::value,
   "Nonlinear least-squares solver: the state type must be a pressio vector wrapper.");

  static_assert
  (::pressio::solvers::concepts::linear_solver_for_nonlinear_least_squares<
   linear_solver_t, state_t>::value,
   "A valid linear solver type must be passed to GN or LM with normal equations");

  // hessian_t is extracted from linear solver
  using hess_t = typename linear_solver_t::matrix_type;

  using operators_t =
    typename std::conditional<
    std::is_same<tag, GaussNewton>::value,
    HessianGradientOperatorsHGApi<hess_t,   grad_t>,
    LMHessianGradientOperatorsHGApi<hess_t, grad_t>
    >::type;

  using corr_mixin = HessianGradientCorrector<operators_t, state_t, linear_solver_t>;
  using type = Solver<tag, corr_mixin>;
};


// ------------------------------------------------------------
// COMPOSE GN or LM with Neq and r/j api, pressio ops and M=I
// ------------------------------------------------------------
template<
  typename system_t,
  typename tag,
  typename linear_solver_t
  >
struct compose<
  system_t, tag,
  mpl::enable_if_t<
    (::pressio::solvers::concepts::system_residual_jacobian<system_t>::value or
     ::pressio::solvers::concepts::system_fused_residual_jacobian<system_t>::value)
    and
    (std::is_same<tag, GaussNewton>::value or std::is_same<tag, LM>::value)
    >,
  linear_solver_t
  >
{
  //if we get here, system_t meets r/j api
  using scalar_t = typename system_t::scalar_type;
  using state_t = typename system_t::state_type;
  using r_t = typename system_t::residual_type;
  using j_t = typename system_t::jacobian_type;
  using grad_t = state_t;

  static_assert
  (::pressio::containers::predicates::is_vector_wrapper<state_t>::value,
   "Nonlinear least-squares solver: the state type must be a pressio vector wrapper.");

  static_assert
  (::pressio::solvers::concepts::linear_solver_for_nonlinear_least_squares<
   linear_solver_t, state_t>::value,
   "A valid linear solver type must be passed to GN or LM with normal equations");

  // hessian_t is extracted from linear solver
  using hess_t = typename linear_solver_t::matrix_type;

  using operators_t =
    typename std::conditional<
    std::is_same<tag, GaussNewton>::value,
    HessianGradientOperatorsRJApiNoWeighting<hess_t, grad_t, r_t, j_t>,
    LMHessianGradientOperatorsRJApi<hess_t, grad_t, r_t, j_t, void,
				    HessianGradientOperatorsRJApiNoWeighting>
    >::type;

  using corr_mixin = HessianGradientCorrector<operators_t, state_t, linear_solver_t>;
  using type = Solver<tag, corr_mixin>;
};

// ---------------------------------------------------------------------------
// *** COMPOSE GN or LM with Neq and r/j api and optional ops and valid M  ***
// ---------------------------------------------------------------------------
template< typename weighting_functor_t, typename tag>
struct composeOperators;

template<>
struct composeOperators<void, GaussNewton>{
  template<typename ...Args>
  using type = HessianGradientOperatorsRJApiNoWeighting<Args...>;
};

template<typename weighting_functor_t>
struct composeOperators<weighting_functor_t, GaussNewton>{
  template<typename ...Args>
  using type = WeightedHessianGradientOperatorsRJApi<Args..., weighting_functor_t>;
};

template<>
struct composeOperators<void, LM>{
  template<typename ...Args>
  using type = LMHessianGradientOperatorsRJApi<
    Args..., HessianGradientOperatorsRJApiNoWeighting>;
};

template<typename weighting_functor_t>
struct composeOperators<weighting_functor_t, LM>{
  template<typename ...Args>
  using type = LMHessianGradientOperatorsRJApi<
    Args..., WeightedHessianGradientOperatorsRJApi, weighting_functor_t>;
};


template<
  typename system_t,
  typename tag,
  typename linear_solver_t,
  typename extra_t
  >
struct compose<
  system_t, tag,
  mpl::enable_if_t<
    (::pressio::solvers::concepts::system_residual_jacobian<system_t>::value or
     ::pressio::solvers::concepts::system_fused_residual_jacobian<system_t>::value) and
    (std::is_same<tag, GaussNewton>::value or std::is_same<tag, LM>::value) and
    !std::is_void<extra_t>::value
    >,
  linear_solver_t, extra_t
  >
{
  //if we get here, system_t meets r/j api
  using scalar_t = typename system_t::scalar_type;
  using state_t = typename system_t::state_type;
  using r_t = typename system_t::residual_type;
  using j_t = typename system_t::jacobian_type;
  using grad_t = state_t;

  static_assert
  (::pressio::containers::predicates::is_vector_wrapper<state_t>::value,
   "Nonlinear least-squares solver: the state type must be a pressio vector wrapper.");

  static_assert
  (::pressio::solvers::concepts::linear_solver_for_nonlinear_least_squares<
   linear_solver_t, state_t>::value,
   "A valid linear solver type must be passed to GN or LM with normal equations");

  // hessian_t is extracted from linear solver
  using hess_t = typename linear_solver_t::matrix_type;

  // if we get here, it means that:
  // - system meets r/j api
  // - we need to compose GN or LM with normal equations
  // - we have one optional template paramter, which can
  //   either be a udops or the weighting

  using ud_ops_t  =
    typename std::conditional<
    ::pressio::solvers::concepts::ops_normal_equations_rj_api<
      extra_t, scalar_t, hess_t, grad_t, j_t, r_t>::value,
    extra_t, void >::type;

  using weighting_functor_t =
    typename std::conditional<
    (::pressio::solvers::concepts::least_squares_weighting_operator<extra_t, r_t, j_t>::value),
    extra_t, void >::type;

  // the extra_t should be one of the two, cannot be both
  static_assert
  ( std::is_void<ud_ops_t>::value or std::is_void<weighting_functor_t>::value,
    "Both ud_ops_t and weighting_functor_t are non-void, something is wrong.");

  using operators_t =
    typename composeOperators<weighting_functor_t, tag>::template type<
    hess_t, grad_t, r_t, j_t, ud_ops_t>;

  using corr_mixin = HessianGradientCorrector<operators_t, state_t, linear_solver_t>;
  using type = Solver<tag, corr_mixin>;
};


template<
  typename system_t,
  typename tag,
  typename linear_solver_t,
  typename ud_ops_t,
  typename weighting_functor_t
  >
struct compose<
  system_t, tag,
  mpl::enable_if_t<
    (::pressio::solvers::concepts::system_residual_jacobian<system_t>::value or
     ::pressio::solvers::concepts::system_fused_residual_jacobian<system_t>::value) and
    (std::is_same<tag, GaussNewton>::value or std::is_same<tag, LM>::value) and
    !std::is_void<ud_ops_t>::value and
    !std::is_void<weighting_functor_t>::value
    >,
  linear_solver_t, ud_ops_t, weighting_functor_t
  >
{
  //if we get here, system_t meets r/j api

  using scalar_t = typename system_t::scalar_type;
  using state_t = typename system_t::state_type;
  using r_t = typename system_t::residual_type;
  using j_t = typename system_t::jacobian_type;
  using grad_t = state_t;

  static_assert
  (::pressio::containers::predicates::is_vector_wrapper<state_t>::value,
   "Nonlinear least-squares solver: the state type must be a pressio vector wrapper.");

  static_assert
  (::pressio::solvers::concepts::linear_solver_for_nonlinear_least_squares<
   linear_solver_t, state_t>::value,
   "A valid linear solver type must be passed to GN or LM with normal equations");

  // hessian_t is extracted from linear solver
  using hess_t = typename linear_solver_t::matrix_type;

  // if we get here, it means that:
  // - system meets r/j api
  // - we need to compose GN or LM with normal equations
  static_assert
  (::pressio::solvers::concepts::ops_normal_equations_rj_api<
   ud_ops_t, scalar_t, hess_t, grad_t, j_t, r_t>::value,
   "The ops type is not admissible for normal equations.");

  static_assert
  (::pressio::solvers::concepts::least_squares_weighting_operator<
   weighting_functor_t, r_t, j_t>::value,
   "The weighting_functor is not admissible");

  using operators_t =
    typename composeOperators<
    weighting_functor_t, tag>::template type<hess_t, grad_t, r_t, j_t, ud_ops_t>;

  using corr_mixin = HessianGradientCorrector<operators_t, state_t, linear_solver_t>;
  using type = Solver<tag, corr_mixin>;
};

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_SOLVERS_NONLINEAR_COMPOSE_HPP_
