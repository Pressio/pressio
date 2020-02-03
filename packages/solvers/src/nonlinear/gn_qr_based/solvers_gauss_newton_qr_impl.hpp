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
			   typename system_t::state_type & stateInOut,
			   typename system_t::state_type & ytrial,
			   typename system_t::residual_type & residual,
			   typename system_t::jacobian_type & jacobian,
			   typename system_t::state_type & correction,
			   typename system_t::state_type & QTResid,
			   qr_solver_t & qrObj,
			   iteration_t maxNonLIt,
			   scalar_t tolerance,
			   std::string & convCondDescr,
			   const ::pressio::solvers::Norm & normType){

  using jacobian_t	= typename system_t::jacobian_type;

  // policy for checking convergence
  using is_converged_t = IsConvergedHelper<converged_when_tag>;

  /* policy for computing line search factor (alpha) such that
   * the update is done with stateInOut = stateInOut + alpha correction
   * alpha = 1 default when user does not want line search
   */
  using lsearch_helper = LineSearchHelper<line_search_t>;
  //-------------------------------------------------------

  constexpr auto one = ::pressio::utils::constants::one<scalar_t>();
  convCondDescr = std::string(is_converged_t::description_);

#ifdef PRESSIO_ENABLE_DEBUG_PRINT
  // get precision before GN
  auto ss = std::cout.precision();
  // set to 14 for the GN prints
  std::cout.precision(14);
  // print GN is starting
  auto fmt1 = utils::io::cyan() + utils::io::underline();
  ::pressio::utils::io::print_stdout(fmt1, "GN with QR:", "criterion:",
				     convCondDescr, utils::io::reset(), "\n");
#endif

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
  // norm of correction
  scalar_t correctionNorm = {};

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
  auto timer = Teuchos::TimeMonitor::getStackedTimer();
  timer->start("QR-based Gauss Newton");
#endif

  iteration_t iStep = 0;
  while (++iStep <= maxNonLIt)
  {
#ifdef PRESSIO_ENABLE_DEBUG_PRINT
    ::pressio::utils::io::print_stdout("\n");
    auto fmt = utils::io::underline();
    ::pressio::utils::io::print_stdout(fmt, "GN step", iStep, utils::io::reset(), "\n");
#endif

    // residual norm for current state
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("norm resid");
#endif
    ComputeNormHelper::template evaluate<void>(residual, normRes, normType);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("norm resid");
#endif
    // store initial residual norm
    if (iStep==1) normRes0 = normRes;


#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("QR factorize");
#endif
    qrObj.computeThin(jacobian);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("QR factorize");
#endif

#ifdef PRESSIO_ENABLE_DEBUG_PRINT
    auto fmt2 = utils::io::magenta() + utils::io::bold();
    ::pressio::utils::io::print_stdout(fmt2, "GN_JSize =",
    ::pressio::solvers::impl::MatrixGetSizeHelper<jacobian_t>::globalRows(jacobian),
    ::pressio::solvers::impl::MatrixGetSizeHelper<jacobian_t>::globalCols(jacobian),
				       utils::io::reset(), "\n");
#endif

    // compute: Q^T Residual
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("QR project");
#endif
    qrObj.applyQTranspose(residual, QTResid);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("QR project");
#endif

    // projected residual norm for current state
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("norm QTResid");
#endif
    ComputeNormHelper::template evaluate<void>(QTResid, normQTRes, normType);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("norm QTResid");
#endif
    // store initial residual norm
    if (iStep==1) normQTRes0 = normQTRes;

    // compute correction: correction
    // by solving R correction = - Q^T Residual
    QTResid.scale(static_cast<scalar_t>(-1));
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("QR R-solve");
    qrObj.solve(QTResid, correction);
    timer->stop("QR R-solve");
#else
    qrObj.solve(QTResid, correction);
#endif

    // norm of the correction
    ComputeNormHelper::template evaluate<void>(correction, correctionNorm, normType);

#ifdef PRESSIO_ENABLE_DEBUG_PRINT
    ::pressio::utils::io::print_stdout(std::scientific,
				    "||R|| =", normRes,
				    "||R||(r) =", normRes/normRes0,
				    "||Q^T R|| =", normQTRes,
				    "||Q^T R||(r) =", normQTRes/normQTRes0,
				    "||dy|| =", correctionNorm,
				    "\n");
#endif

    // exit with error if NaNs detected in solution update correction
    if (std::isnan(correctionNorm))
    {
      throw std::runtime_error(
        "Nonlinear solver: QR-based Gauss Newton: NaNs detected in solution update correction");
    }

    // compute multiplicative factor if needed
    lsearch_helper::template evaluate<void>(alpha, stateInOut, ytrial, correction, residual, jacobian, sys);

    // solution update: stateInOut = stateInOut + alpha*correction
    ::pressio::containers::ops::do_update(stateInOut, one, correction, alpha);

    // check convergence (whatever method user decided)
    const auto flag = is_converged_t::evaluate(stateInOut, correction,
					       correctionNorm, normRes, normRes0,
					       normQTRes, normQTRes0,
					       iStep, maxNonLIt, tolerance);
    if (flag) break;

    sys.residual(stateInOut, residual);
    sys.jacobian(stateInOut, jacobian);

  }//loop

#if defined PRESSIO_ENABLE_DEBUG_PRINT
  std::cout.precision(ss);
#endif

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
  timer->stop("QR-based Gauss Newton");
#endif

}


}}}} //end namespace pressio::solvers::iterative::implo
#endif
