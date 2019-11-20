/*
//@HEADER
// ************************************************************************
//
// solvers_gauss_newton.hpp
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

#ifndef SOLVERS_GAUSS_NEWTON_HESSIAN_GRADIENT_API_HPP_
#define SOLVERS_GAUSS_NEWTON_HESSIAN_GRADIENT_API_HPP_

#include "../solvers_fwd.hpp"
#include "../base/solvers_nonlinear_base.hpp"
#include "../base/solvers_iterative_base.hpp"
#include "../nonlinear/helper_policies/solvers_converged_criterior_policy.hpp"
#include "../nonlinear/helper_policies/solvers_norm_helper_policy.hpp"
#include "../nonlinear/helper_policies/solvers_line_search_policy.hpp"
#include "../nonlinear/helper_policies/solvers_get_matrix_size_helper.hpp"

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{
namespace experimental{

template <
  typename system_type,
  typename linear_solver_type,
  typename scalar_type,
  typename line_search_type,
  typename convergence_when_t
  >
class GaussNewtonHessianGradientApi
  : public NonLinearSolverBase<
  GaussNewtonHessianGradientApi<system_type,
				linear_solver_type,
				scalar_type,
				line_search_type,
				convergence_when_t>
  >,
    public IterativeBase<
  GaussNewtonHessianGradientApi<system_type,
				linear_solver_type,
				scalar_type,
				line_search_type,
				convergence_when_t>, scalar_type
  >
{

  using this_t = GaussNewtonHessianGradientApi<
    system_type, linear_solver_type, scalar_type, line_search_type, convergence_when_t>;

  // need to be friend of base (crpt)
  using non_lin_sol_base_t = NonLinearSolverBase<this_t>;
  friend non_lin_sol_base_t;

  // iterative base
  using iterative_base_t = IterativeBase<this_t, scalar_type>;
  using typename iterative_base_t::iteration_t;

  // typedefs from the system
  using state_t    = typename system_type::state_type;
  using hessian_t  = typename system_type::hessian_type;
  using gradient_t = typename system_type::gradient_type;

  static_assert( std::is_same<state_t, gradient_t>::value,
		 "Currently, for the Gauss-Newton solver with hessian/gradient API, \
the type of the gradient approximation, i.e. jacobian^T * residual, must be the same as the state type.\
If you get this error, most likely it is because of  the typedefs inside your system class.");

  static_assert( std::is_same<line_search_type, ::pressio::solvers::iterative::gn::noLineSearch>::value,
		 "Currently, the Gauss-Newton solver with hessian/gradient API, does not support a line-search.");

  // using norm_t = typename NormSelectorHelper<converged_when_t>::norm_t;
  // using norm_evaluator_t = ComputeNormHelper<norm_t>;
  // using is_converged_t = IsConvergedHelper<converged_when_t>;

  // --- data members ---
  linear_solver_type & linSolver_ = {};
  hessian_t hess_		  = {};
  gradient_t grad_		  = {};
  // delta is the correction
  state_t delta_		  = {};
  // ytrail needed if/when line search is used
  state_t ytrial_		  = {};

  // old norm of state
  scalar_type normState0_	  = {0};

  // old and current norm of J^T R
  scalar_type normGrad0_	  = {0};
  scalar_type normGrad_	  = {0};

  // norm of the correction
  scalar_type norm_delta_	  = {0};

public:
  GaussNewtonHessianGradientApi() = delete;
  GaussNewtonHessianGradientApi(const GaussNewtonHessianGradientApi &) = delete;
  ~GaussNewtonHessianGradientApi() = default;

  GaussNewtonHessianGradientApi(const system_type  & system,
				const state_t	  & stateIn,
				linear_solver_type & linearSolverIn)
    : linSolver_(linearSolverIn),
      hess_( system.createHessianObject(stateIn) ),
      grad_( system.createGradientObject(stateIn) ),
      delta_(stateIn),
      ytrial_(stateIn)
  {}

private:

  void solveImpl(const system_type & sys, state_t & stateInOut)
  {
//     constexpr auto one = ::pressio::utils::constants::one<scalar_type>();
//     constexpr auto negOne = ::pressio::utils::constants::negOne<scalar_type>();
//     const auto convCondDescr = std::string(is_converged_t::description_);

//     const auto maxNonLIt  = iterative_base_t::maxIters_;
//     const auto tolerance  = iterative_base_t::tolerance_;

//     // compute J^T J and J^T R in one shot
//     sys.computeHessianAndGradient(stateInOut, hess_, grad_);

//     // compute the initial norm of y (the state)
//     norm_evaluator_t::evaluate(yStateInOut, normState0_);
//     norm_delta_ = {0};

//     iteration_t iStep = 0;
//     while (++iStep <= maxNonLIt)
//     {
//       // norm of J^T R
//       norm_evaluator_t::evaluate(grad_, normGrad_);
//       if (iStep==1) normGrad0_ = normGrad_;

//       linSolver.solveAllowMatOverwrite(hess_, grad_, delta_);
//       norm_evaluator_t::evaluate(delta_, norm_delta_);

// #ifdef PRESSIO_ENABLE_DEBUG_PRINT
//       ::pressio::utils::io::print_stdout(std::scientific,
// 					 "||J^T R|| =", normGrad_,
// 					 "||J^T R||(r) =", normGrad_/normGrad0_,
// 					 "||dy|| =", norm_delta_,
// 					 utils::io::reset(),
// 					 "\n");
// #endif

//       // exit with error if NaNs detected in solution update dy
//       if (std::isnan(norm_delta_)){
//       	throw std::runtime_error("Nonlinear solver: NEQ-based Gausss Newton: NaNs detected in solution update dy");
//       }

//       // // compute multiplicative factor if needed (not supported for now)
//       // lsearch_helper_t::evaluate(alpha, yStateInOut, ytrial, delta_, resid, jacob, sys);

//       // solution update: y = y + alpha*dy
//       ::pressio::containers::ops::do_update(yStateInOut, one, delta_, alpha);

//       // check convergence (whatever method user decided)
//       // TODO: we need to make sure the convergence criterion is compatible for this solver
//       // since here we only have the norm of J^T R and not the norm of R
//       const auto flag = is_converged_t::evaluate(yStateInOut, delta_, norm_delta_,
// 						 one, one, /* dummy entry for norm of R */
// 						 normGrad, normGrad0, iStep,
//       						 maxNonLIt, tolerance);

//       // if we have converged, query the observer
//       if (flag) {
//       	break;
//       }

//       // store new norm into old variable
//       normState0_ = norm_delta_;
//       // compute hessian and gradient
//       sys.computeHessianAndGradient(yStateInOut, hess_, HTResid_);

//     }//loop

  }
};

}}}}}//end namespace pressio::solvers::iterative::impl::experimental
#endif
