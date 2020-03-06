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

#include "../../helpers/solvers_converged_criterior_policy.hpp"
#include "../../helpers/solvers_norm_dispatcher.hpp"
#include "../../helpers/solvers_line_search_policy.hpp"
#include "../../helpers/solvers_get_matrix_size_helper.hpp"

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template <
  typename system_type,
  typename linear_solver_type,
  typename scalar_type,
  typename line_search_type,
  typename converged_when
  >
class GaussNewtonHessianGradientApi
  : public NonLinearSolverBase<
  GaussNewtonHessianGradientApi<system_type,
				linear_solver_type,
				scalar_type,
				line_search_type,
				converged_when>
  >,
    public IterativeBase<
  GaussNewtonHessianGradientApi<system_type,
				linear_solver_type,
				scalar_type,
				line_search_type,
				converged_when>, scalar_type
  >
{

  using this_t = GaussNewtonHessianGradientApi<
    system_type, linear_solver_type, scalar_type, line_search_type, converged_when>;

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

  using is_converged_t = IsConvergedHelper<converged_when>;

  // --- data members ---
  linear_solver_type & linSolver_ = {};
  ::pressio::solvers::Norm normType_ = {};

  hessian_t hessian_		  = {};
  gradient_t gradient_		  = {};
  // delta is the correction
  state_t correction_		  = {};

  // norm of correction
  scalar_type normCorrection_	  = {0};

  // old and current norm of J^T R
  scalar_type normGradient0_	  = {0};
  scalar_type normGradient_	  = {0};

  scalar_type normResidual0_	  = {0};
  scalar_type normResidual_	  = {0};

  NormDispatcher<void> normDispatcher_ = {};

public:
  GaussNewtonHessianGradientApi() = delete;
  GaussNewtonHessianGradientApi(const GaussNewtonHessianGradientApi &) = delete;
  ~GaussNewtonHessianGradientApi() = default;

  GaussNewtonHessianGradientApi(const system_type  & system,
				const state_t	  & stateIn,
				linear_solver_type & linearSolverIn,
				const ::pressio::solvers::Norm normType = ::pressio::solvers::defaultNormType)
    : linSolver_(linearSolverIn),
      normType_(normType),
      hessian_( system.createHessianObject(stateIn) ),
      gradient_( system.createGradientObject(stateIn) ),
      correction_(stateIn)
  {}

private:

  void solveImpl(const system_type & sys, state_t & stateInOut)
  {
    constexpr auto one = ::pressio::utils::constants::one<scalar_type>();
    constexpr auto negOne = ::pressio::utils::constants::negOne<scalar_type>();
    const auto convCondDescr = std::string(is_converged_t::description_);

    const auto maxNonLIt  = iterative_base_t::maxIters_;
    const auto tolerance  = iterative_base_t::tolerance_;

    // compute J^T J and J^T R in one shot
    sys.computeHessianAndGradient(stateInOut, hessian_, gradient_, normType_, normResidual_);

    ::pressio::ops::scale(gradient_, negOne);
    // gradient_.scale(negOne);

    normResidual0_ = normResidual_;

    // compute the initial norm of y (the state)
    normDispatcher_.evaluate(stateInOut, normCorrection_, normType_);
    normCorrection_ = {0};

    // the alpha is 1, but this can change from the line search
    scalar_type alpha = one;

    iteration_t iStep = 0;
    while (++iStep <= maxNonLIt)
    {
      // norm of J^T R
      normDispatcher_.evaluate(gradient_, normGradient_, normType_);
      if (iStep==1) normGradient0_ = normGradient_;

      linSolver_.solveAllowMatOverwrite(hessian_, gradient_, correction_);
      normDispatcher_.evaluate(correction_, normCorrection_, normType_);

#ifdef PRESSIO_ENABLE_DEBUG_PRINT
      ::pressio::utils::io::print_stdout(std::scientific,
					 "||R|| =", normResidual_,
					 "||R||(r) =", normResidual_/normResidual0_,
					 "||J^T R|| =", normGradient_,
					 "||J^T R||(r) =", normGradient_/normGradient0_,
					 "||dy|| =", normCorrection_,
					 utils::io::reset(),
					 "\n");
#endif

      // exit with error if NaNs detected in solution update dy
      if (std::isnan(normCorrection_)){
      	throw std::runtime_error("Nonlinear solver: NEQ-based Gausss Newton: NaNs detected in solution update dy");
      }

      // // compute multiplicative factor if needed (not supported for now)
      // lsearch_helper_t::evaluate(alpha, stateInOut, ytrial, correction_, resid, jacob, sys);

      // solution update: state = state + alpha*correction
      ::pressio::ops::do_update(stateInOut, one, correction_, alpha);

      // check convergence (whatever method user decided)
      // TODO: we need to make sure the convergence criterion is compatible for this solver
      // since here we only have the norm of J^T R and not the norm of R
      const auto flag = is_converged_t::evaluate(stateInOut, correction_, normCorrection_,
						 normResidual_, normResidual0_,
						 normGradient_, normGradient0_, iStep,
      						 maxNonLIt, tolerance);

      // if we have converged, query the observer
      if (flag) {
      	break;
      }

      // compute hessian and gradient
      sys.computeHessianAndGradient(stateInOut, hessian_, gradient_, normType_, normResidual_);
      ::pressio::ops::scale(gradient_, negOne);
      // gradient_.scale(negOne);
    }//loop
  }
};

}}}}//end namespace pressio::solvers::iterative::impl
#endif
