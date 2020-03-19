/*
//@HEADER
// ************************************************************************
//
// solvers_gauss_newton_normal_eq_impl.hpp
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

#ifndef SOLVERS_GAUSS_NEWTON_NORMAL_EQ_IMPL_HPP
#define SOLVERS_GAUSS_NEWTON_NORMAL_EQ_IMPL_HPP

#include "../../helpers/solvers_converged_criterior_policy.hpp"
#include "../../helpers/solvers_line_search_policy.hpp"
#include "../../helpers/solvers_residual_observer_when_solver_converged.hpp"
#include "../../helpers/solvers_residual_observer_each_solver_step.hpp"
#include "../../helpers/solvers_get_matrix_size_helper.hpp"

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template <typename hessian_t>
mpl::enable_if_t< !::pressio::containers::meta::is_dense_matrix_wrapper_eigen<hessian_t>::value >
printConditionNumber(const hessian_t & H)
{
  //no op
}

template <typename hessian_t>
mpl::enable_if_t< ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<hessian_t>::value >
printConditionNumber(const hessian_t & H)
{
  using eig_mat = typename ::pressio::containers::details::traits<hessian_t>::wrapped_t;
  Eigen::JacobiSVD<eig_mat> svd(*H.data());
  const auto cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);

  auto fmt2 = utils::io::magenta() + utils::io::bold();
  ::pressio::utils::io::print_stdout(fmt2, "GN_HessianCondition =", cond, utils::io::reset(), "\n");
}



template <
  typename line_search_t,
  typename converged_when_tag,
  typename system_t,
  typename hessian_t,
  typename lin_solver_t,
  typename iteration_t,
  typename scalar_t,
  typename hessian_dispatcher_t,
  typename gradient_dispatcher_t,
  typename norm_dispatcher_t,
  typename observer_t = utils::impl::empty
  >
void gauss_newton_neq_solve(const system_t & sys,
			    typename system_t::state_type & stateInOut,
			    typename system_t::state_type & ytrial,
			    typename system_t::residual_type & residual,
			    typename system_t::jacobian_type & jacobian,
			    typename system_t::state_type & correction,
			    typename system_t::state_type & gradient,
			    hessian_t & hessian,
			    lin_solver_t & linSolver,
			    iteration_t maxNonLIt,
			    scalar_t tolerance,
			    const observer_t * observer,
			    std::string & convCondDescr,
			    ::pressio::solvers::Norm normType,
			    const hessian_dispatcher_t & hessianDispatcher,
			    const gradient_dispatcher_t & gradientDispatcher,
			    const norm_dispatcher_t & normDispatcher)
{

  using residual_t	= typename system_t::residual_type;
  //using jacobian_t	= typename system_t::jacobian_type;

  // policy to checking convergence
  using is_converged_t = IsConvergedHelper<converged_when_tag>;

  // policy to observing residual at each GN step
  using residual_observer_each_step = ResidualObserverEachSolverStep<observer_t, residual_t>;

  // policy to observing residual when converged before exiting
  using residual_observer_when_conv = ResidualObserverWhenSolverConverged<observer_t, residual_t>;

  /* policy for computing line search factor (alpha) such that
   * the update is done with y = y + alpha correction
   * alpha = 1 default when user does not want line search*/
  using lsearch_helper_t = LineSearchHelper<line_search_t>;

  //-------------------------------------------------------

  constexpr auto one = ::pressio::utils::constants::one<scalar_t>();
  constexpr auto negOne = ::pressio::utils::constants::negOne<scalar_t>();
  convCondDescr = std::string(is_converged_t::description_);

#ifdef PRESSIO_ENABLE_DEBUG_PRINT
  // get precision before GN
  auto ss = std::cout.precision();
  // set to 14 for the GN prints
  std::cout.precision(14);
  auto reset = utils::io::reset();
  auto fmt1 = utils::io::cyan() + utils::io::underline();
  ::pressio::utils::io::print_stdout(fmt1, "GN normal eqns:", "criterion:",
				     convCondDescr, reset, "\n");
#endif

  // alpha for taking steps
  scalar_t alpha = {};
  // storing residual norm
  scalar_t normRes = {};
  scalar_t normRes0 = {};
  // storing gradient norm
  scalar_t normGrad = {};
  scalar_t normGrad0 = {};
  // norm of the correction
  scalar_t correctionNorm = {0};

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
  auto timer = Teuchos::TimeMonitor::getStackedTimer();
  timer->start("NEQ-based Gausss Newton");
#endif

  iteration_t iStep = 0;
  while (++iStep <= maxNonLIt)
  {

    // call residual observer at each gauss step (no op for dummy case)
    residual_observer_each_step::evaluate(observer, residual, iStep);

#ifdef PRESSIO_ENABLE_DEBUG_PRINT
    ::pressio::utils::io::print_stdout("\n");
    auto fmt = utils::io::underline();
    ::pressio::utils::io::print_stdout(fmt, "GN step", iStep, utils::io::reset(), "\n");
#endif

    // residual norm for current state
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("norm resid");
#endif
    normDispatcher.evaluate(residual, normRes, normType);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("norm resid");
#endif
    // store initial residual norm
    if (iStep==1) normRes0 = normRes;

    // compute LHS: J^T*J
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("hessian");
#endif
    hessianDispatcher.evaluate(jacobian, hessian);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("hessian");
#endif

#ifdef PRESSIO_ENABLE_DEBUG_PRINT
    auto fmt2 = utils::io::magenta() + utils::io::bold();
    ::pressio::utils::io::print_stdout(fmt2, "GN_HessianSize =",
				       hessian.extent(0), hessian.extent(1),
				       utils::io::reset(), "\n");
    printConditionNumber(hessian);
#endif

    // compute RHS: J^T*res
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("gradient");
#endif
    gradientDispatcher.evaluate(jacobian, residual, gradient);
    ::pressio::ops::scale(gradient, negOne);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("gradient");
#endif

    // gradient norm for current state
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("norm gradient");
#endif
    normDispatcher.evaluate(gradient, normGrad, normType);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("norm gradient");
#endif
    // store initial residual norm
    if (iStep==1) normGrad0 = normGrad;


    // solve normal equations
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("solve normeq");
#endif
    linSolver.solveAllowMatOverwrite(hessian, gradient, correction);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("solve normeq");
#endif
    // compute norm of the correction
    normDispatcher.evaluate(correction, correctionNorm, normType);

#ifdef PRESSIO_ENABLE_DEBUG_PRINT
    ::pressio::utils::io::print_stdout(std::scientific,
				    "||R|| =", normRes,
				    "||R||(r) =", normRes/normRes0,
				    "||J^T R|| =", normGrad,
				    "||J^T R||(r) =", normGrad/normGrad0,
				    "||dy|| =", correctionNorm,
				    utils::io::reset(),
				    "\n");
#endif

    // exit with error if NaNs detected in solution update dy
    if (std::isnan(correctionNorm)){
      throw std::runtime_error(
        "Nonlinear solver: NEQ-based Gausss Newton: NaNs detected in solution update dy");
    }

    // compute multiplicative factor if needed
    lsearch_helper_t::evaluate(alpha, stateInOut, ytrial, correction,
			       residual, jacobian, sys,
			       gradientDispatcher, normDispatcher);

    // solution update: y = y + alpha*correction
    ::pressio::ops::do_update(stateInOut, one, correction, alpha);

    // check convergence (whatever method user decided)
    const auto flag = is_converged_t::evaluate(stateInOut, correction, correctionNorm,
					       normRes, normRes0,
					       normGrad, normGrad0,
					       iStep, maxNonLIt, tolerance);

    // if we have converged, query the observer
    if (flag) {
      // observe residual (no op for dummy case)
      residual_observer_when_conv::evaluate(observer, residual);
      break;
    }

    // compute residual and jacobian
    sys.residual(stateInOut, residual);
    sys.jacobian(stateInOut, jacobian);

  }//loop

#if defined PRESSIO_ENABLE_DEBUG_PRINT
  std::cout.precision(ss);
  ::pressio::utils::io::print_stdout(std::fixed);
#endif

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
  timer->stop("NEQ-based Gausss Newton");
#endif

}// end

}}}} //end namespace pressio::solvers::iterative::impl
#endif
