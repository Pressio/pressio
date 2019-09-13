/*
//@HEADER
// ************************************************************************
//
// solvers_gauss_newton_qr_impl.hpp
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

#ifndef SOLVERS_GAUSS_NEWTON_QR_IMPL_HPP
#define SOLVERS_GAUSS_NEWTON_QR_IMPL_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../helper_policies/solvers_converged_criterior_policy.hpp"
#include "../helper_policies/solvers_norm_helper_policy.hpp"
#include "../helper_policies/solvers_line_search_policy.hpp"
#include "../helper_policies/solvers_get_matrix_size_helper.hpp"

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template <
  typename line_search_t,
  typename converged_when_tag,
  typename system_t,
  typename qr_solver_t,
  typename iteration_t,
  typename scalar_t
  >
void gauss_newton_qr_solve(const system_t & sys,
			   typename system_t::state_type & y,
			   typename system_t::state_type & ytrial,
			   typename system_t::residual_type & resid,
			   typename system_t::jacobian_type & jacob,
			   typename system_t::state_type & dy,
			   typename system_t::state_type & QTResid,
			   qr_solver_t & qrObj,
			   iteration_t maxNonLIt,
			   scalar_t tolerance,
			   scalar_t & normO,
			   scalar_t & norm_dy,
			   std::string & convCondDescr){

  using jacobian_t	= typename system_t::jacobian_type;

  // find out which norm to use
  using norm_t = typename NormSelectorHelper<converged_when_tag>::norm_t;

  // policy for  evaluating the norm of a vector
  using norm_evaluator_t = ComputeNormHelper<norm_t>;

  // policy for checking convergence
  using is_converged_t = IsConvergedHelper<converged_when_tag>;

  /* policy for computing line search factor (alpha) such that
   * the update is done with y = y + alpha dy
   * alpha = 1 default when user does not want line search
   */
  using lsearch_helper = LineSearchHelper<line_search_t>;
  //-------------------------------------------------------

  // alpha for taking steps
  scalar_t alpha = {};
  // residual norm
  scalar_t normRes = {};
  // initial residual norm
  scalar_t normRes0 = {};
  // residual norm
  scalar_t normQTRes = {};
  // initial residual norm
  scalar_t normQTRes0 = {};

  constexpr auto one = ::pressio::utils::constants::one<scalar_t>();
  convCondDescr = std::string(is_converged_t::description_);

#ifdef DEBUG_PRINT
  // get precision before GN
  auto ss = std::cout.precision();
  // set to 14 for the GN prints
  std::cout.precision(14);
  // print GN is starting
  auto fmt1 = utils::io::cyan() + utils::io::underline();
  ::pressio::utils::io::print_stdout(fmt1, "GN with QR:", "criterion:",
				     convCondDescr, utils::io::reset(), "\n");
#endif

  // compute (whatever type) norm of y
  norm_evaluator_t::evaluate(y, normO);
  norm_dy = {};

#ifdef HAVE_TEUCHOS_TIMERS
  auto timer = Teuchos::TimeMonitor::getStackedTimer();
  timer->start("QR-based Gausss Newton");
#endif

  iteration_t iStep = 0;
  while (++iStep <= maxNonLIt)
  {
#ifdef DEBUG_PRINT
    ::pressio::utils::io::print_stdout("\n");
    auto fmt = utils::io::underline();
    ::pressio::utils::io::print_stdout(fmt, "GN step", iStep,
				    utils::io::reset(), "\n");
#endif

    // residual norm for current state
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("norm resid");
    norm_evaluator_t::evaluate(resid, normRes);
    timer->stop("norm resid");
#else
    norm_evaluator_t::evaluate(resid, normRes);
#endif
    // store initial residual norm
    if (iStep==1) normRes0 = normRes;


#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("QR factorize");
    qrObj.computeThin(jacob);
    timer->stop("QR factorize");
#else
    qrObj.computeThin(jacob);
#endif

#ifdef DEBUG_PRINT
    auto fmt2 = utils::io::magenta() + utils::io::bold();
    ::pressio::utils::io::print_stdout(fmt2, "GN_JSize =",
    ::pressio::solvers::impl::MatrixGetSizeHelper<jacobian_t>::globalRows(jacob),
    ::pressio::solvers::impl::MatrixGetSizeHelper<jacobian_t>::globalCols(jacob),
				    utils::io::reset(),
				    "\n");
#endif

    // compute: Q^T Residual
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("QR project");
    qrObj.project(resid, QTResid);
    timer->stop("QR project");
#else
    qrObj.project(resid, QTResid);
#endif

    // projected residual norm for current state
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("norm QTResid");
#endif
    norm_evaluator_t::evaluate(QTResid, normQTRes);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("norm QTResid");
#endif

    // store initial residual norm
    if (iStep==1) normQTRes0 = normQTRes;

    // compute correction: dy
    // by solving R dy = - Q^T Residual
    QTResid.scale(static_cast<scalar_t>(-1));
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("QR R-solve");
    qrObj.solve(QTResid, dy);
    timer->stop("QR R-solve");
#else
    qrObj.solve(QTResid, dy);
#endif

    // norm of the correction
    norm_evaluator_t::evaluate(dy, norm_dy);

#ifdef DEBUG_PRINT
    ::pressio::utils::io::print_stdout(std::scientific,
				    "||R|| =", normRes,
				    "||R||(r) =", normRes/normRes0,
				    "||Q^T R|| =", normQTRes,
				    "||Q^T R||(r) =", normQTRes/normQTRes0,
				    "||dy|| =", norm_dy,
				    "\n");
#endif

    // compute multiplicative factor if needed
    lsearch_helper::evaluate(alpha, y, ytrial, dy, resid, jacob, sys);

    // solution update: y = y + alpha*dy
    ::pressio::containers::ops::do_update(y, one, dy, alpha);

    // check convergence (whatever method user decided)
    const auto flag = is_converged_t::evaluate(y, dy,
					       norm_dy, normRes, normRes0,
					       normQTRes, normQTRes0,
					       iStep, maxNonLIt, tolerance);
    if (flag) break;

    // store new norm into old variable
    normO = norm_dy;

    sys.residual(y, resid);
    sys.jacobian(y, jacob);

  }//loop

#if defined DEBUG_PRINT
  std::cout.precision(ss);
#endif

#ifdef HAVE_TEUCHOS_TIMERS
  timer->stop("QR-based Gausss Newton");
#endif

}


}}}} //end namespace pressio::solvers::iterative::implo
#endif
