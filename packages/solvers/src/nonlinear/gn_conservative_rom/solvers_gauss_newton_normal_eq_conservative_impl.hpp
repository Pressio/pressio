/*
//@HEADER
// ************************************************************************
//
// solvers_gauss_newton_normal_eq_conservative_impl.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef SOLVERS_GAUSS_NEWTON_NORMAL_EQ_CONSERVATIVE_IMPL_HPP
#define SOLVERS_GAUSS_NEWTON_NORMAL_EQ_CONSERVATIVE_IMPL_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../solvers_meta_static_checks.hpp"
#include "../helper_policies/solvers_converged_criterior_policy.hpp"
#include "../helper_policies/solvers_hessian_helper_policy.hpp"
#include "../helper_policies/solvers_jacob_res_product_policy.hpp"
#include "../helper_policies/solvers_norm_helper_policy.hpp"
#include "../helper_policies/solvers_line_search_policy.hpp"
#include "../helper_policies/solvers_get_matrix_size_helper.hpp"

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template <
  typename system_t,
  typename iteration_t,
  typename scalar_t,
  typename lin_solver_t,
  typename line_search_t,
  typename converged_when_tag,
  typename cbar_t,
  typename mat_t,
  ::pressio::mpl::enable_if_t<
    //::pressio::solvers::details::system_traits<system_t>::is_system and
    containers::meta::is_vector_wrapper_eigen<
      typename system_t::state_type>::value and
    containers::meta::is_vector_wrapper<
      typename system_t::residual_type>::value
    > * =nullptr
  >
void gauss_newtom_neq_conserv_solve(const system_t & sys,
				    typename system_t::state_type & y,
				    typename system_t::state_type & ytrial,
				    typename system_t::residual_type & resid,
				    typename system_t::jacobian_type & jacob,
				    iteration_t maxNonLIt,
				    scalar_t tolerance,
				    typename system_t::state_type & dy,
				    lin_solver_t & linSolver,
				    scalar_t & normO,
				    scalar_t & norm_dy,
				    const cbar_t & cbarT,
				    mat_t & jTj,
				    mat_t & jTcbarT,
				    mat_t & cbarJ,
				    mat_t & zero,
				    typename system_t::residual_type & cbarTlambda,
				    typename system_t::state_type & jTr2,
				    typename system_t::state_type & cbarR,
				    mat_t & A,
				    typename system_t::state_type & b,
				    typename system_t::state_type & lambda,
				    typename system_t::state_type & y2,
				    std::string & convCondDescr)
{

  // find out which norm to use
  using norm_t = typename NormSelectorHelper<converged_when_tag>::norm_t;

  // functor for evaluating the norm of a vector
  using norm_evaluator_t = ComputeNormHelper<norm_t>;

  // functor for checking convergence
  using is_conv_helper_t = IsConvergedHelper<converged_when_tag>;

  // /* functor for computing line search factor (alpha) such that
  //  * the update is done with y = y + alpha dy
  //  * alpha = 1 default when user does not want line search
  //  */
  // using lsearch_helper = LineSearchHelper<line_search_t>;
  // lsearch_helper lineSearchHelper;
  //-------------------------------------------------------

  // alpha for taking steps
  scalar_t alpha = static_cast<scalar_t>(1);
  // storing residual norm
  scalar_t normRes = {};
  scalar_t normRes0 = {};
  // storing projected residual norm
  scalar_t normJTRes = {};
  scalar_t normJTRes0 = {};

  convCondDescr = std::string(is_conv_helper_t::description_);

#ifdef DEBUG_PRINT
  // get precision
  auto ss = std::cout.precision();
  // set to 14 for prints
  std::cout.precision(14);

  auto reset = utils::io::reset();
  auto fmt1 = utils::io::cyan() + utils::io::underline();
  ::pressio::utils::io::print_stdout(fmt1, "GN normal eqns conserv:",
				     "criterion:",
				     convCondDescr, reset, "\n");
#endif

  // compute (whatever type) norm of y
  norm_evaluator_t::evaluate(y, normO);
  norm_dy = static_cast<scalar_t>(0);
  auto normLambda = static_cast<scalar_t>(0);
  auto normCbarR = static_cast<scalar_t>(0);

  typename system_t::state_type dy_y(y.size());
  typename system_t::state_type dy_lambda(lambda.size());

#ifdef HAVE_TEUCHOS_TIMERS
  auto timer = Teuchos::TimeMonitor::getStackedTimer();
  timer->start("Gausss Newton Conserv");
#endif

  iteration_t iStep = 0;
  while (++iStep <= maxNonLIt)
  {

#ifdef DEBUG_PRINT
    ::pressio::utils::io::print_stdout("\n");
    auto fmt = utils::io::underline();
    ::pressio::utils::io::print_stdout(fmt, "step", iStep,
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

    // assemble LHS
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("lhs");
#endif
    ::pressio::containers::ops::dot_self(jacob, jTj);
    // ::pressio::utils::io::print_stdout(*jTj.data(), "\n");
    // ::pressio::utils::io::print_stdout("--------------\n");

    ::pressio::containers::ops::dot(jacob, cbarT, jTcbarT);
    // ::pressio::utils::io::print_stdout(*jTcbarT.data(), "\n");
    // ::pressio::utils::io::print_stdout("--------------\n");

    ::pressio::containers::ops::dot(cbarT, jacob, cbarJ);
    // ::pressio::utils::io::print_stdout(*cbarJ.data(), "\n");
    // ::pressio::utils::io::print_stdout("--------------\n");

    A.data()->block(0, 0, jTj.rows(), jTj.cols()) = *jTj.data();
    A.data()->block(0, jTj.cols(), jTcbarT.rows(), jTcbarT.cols()) = *jTcbarT.data();
    A.data()->block(jTj.rows(), 0, cbarJ.rows(), cbarJ.cols()) = *cbarJ.data();
    A.data()->block(jTj.rows(), jTj.cols(), zero.rows(), zero.cols()) = *zero.data();

    // ::pressio::utils::io::print_stdout(*A.data(), "\n");
    // ::pressio::utils::io::print_stdout("--------------\n");

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("lhs");
#endif

#ifdef DEBUG_PRINT
    // auto fmt1 = utils::io::magenta() + utils::io::bold();
    // ::pressio::utils::io::print_stdout(fmt1, "GN_JSize =",
    // ::pressio::solvers::impl::MatrixGetSizeHelper<jac_t>::globalRows(jacob),
    // ::pressio::solvers::impl::MatrixGetSizeHelper<jac_t>::globalCols(jacob),
    // 				    "\n");
    // // this print only works when hessian is a shared mem matrix
    // ::pressio::utils::io::print_stdout(fmt1, "GN_HessianSize =",
    // 				    H.rows(), H.cols(),
    // 				    utils::io::reset(), "\n");
#endif

    // compute RHS
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("rhs");
#endif

    ::pressio::containers::ops::dot(cbarT, resid, cbarR);
    norm_evaluator_t::evaluate(cbarR, normCbarR);

    ::pressio::containers::ops::product(cbarT, lambda, cbarTlambda);
    resid.data()->update(1.0, *cbarTlambda.data(), 1.0);
    ::pressio::containers::ops::dot(jacob, resid, jTr2);

    auto negOne = static_cast<scalar_t>(-1);
    jTr2.scale(negOne);
    cbarR.scale(negOne);
    b.data()->block(0, 0, jTr2.size(), 1) = *jTr2.data();
    b.data()->block(jTr2.size(), 0, cbarR.size(), 1) = *cbarR.data();

    // ::pressio::utils::io::print_stdout("-----cbarTlambda-----\n");
    // ::pressio::utils::io::print_stdout(*cbarTlambda.data(), "\n");
    // ::pressio::utils::io::print_stdout("--------------\n");

    // ::pressio::utils::io::print_stdout("-----jTr2-----\n");
    // ::pressio::utils::io::print_stdout(*jTr2.data(), "\n");
    // ::pressio::utils::io::print_stdout("--------------\n");

    // ::pressio::utils::io::print_stdout(*cbarR.data(), "\n");
    // ::pressio::utils::io::print_stdout("--------------\n");
    // ::pressio::utils::io::print_stdout(*b.data(), "\n");
    // ::pressio::utils::io::print_stdout("--------------\n");

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("rhs");
#endif

    // projected residual norm for current state
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("norm JTR");
#endif
    norm_evaluator_t::evaluate(jTr2, normJTRes);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("norm JTR");
#endif

    // store initial residual norm
    if (iStep==1) normJTRes0 = normJTRes;

    // solve normal equations
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("solve normeq");
#endif

    linSolver.solve(A, b, dy);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("solve normeq");
#endif

    *dy_y.data() = dy.data()->block(0, 0, y.size(), 1);
    *dy_lambda.data() = dy.data()->block(y.size(), 0, lambda.size(), 1);

    // norm of the correction
    norm_evaluator_t::evaluate(dy_y, norm_dy);
    norm_evaluator_t::evaluate(dy_lambda, normLambda);

#ifdef DEBUG_PRINT
    ::pressio::utils::io::print_stdout(std::scientific,
				    "||R|| =", normRes,
				    "||R||(r) =", normRes/normRes0,
				    "||J^T R|| =", normJTRes,
				    "||J^T R||(r) =", normJTRes/normJTRes0,
				    "||cbar_R|| = ", normCbarR,
				    "||dy|| =", norm_dy,
				    "||dlambda|| =", normLambda,
				    utils::io::reset(), "\n");
#endif

    // // compute multiplicative factor if needed
    // lineSearchHelper(alpha, y, ytrial, dy, resid, jacob, sys);

    // ::pressio::utils::io::print_stdout("-----dy-----\n");
    // ::pressio::utils::io::print_stdout(*dy.data(), "\n");
    // ::pressio::utils::io::print_stdout("--------------\n");

    // ::pressio::utils::io::print_stdout("-----y-----\n");
    // ::pressio::utils::io::print_stdout(*y.data(), "\n");
    // ::pressio::utils::io::print_stdout("--------------\n");

    // ::pressio::utils::io::print_stdout("-----y2-----\n");
    // ::pressio::utils::io::print_stdout(*y2.data(), "\n");
    // ::pressio::utils::io::print_stdout("--------------\n");

    y2 = y2 + alpha * dy;

    // solution update
    *y.data() = y2.data()->block(0, 0, y.size(), 1);
    *lambda.data() = y2.data()->block(y.size(), 0, lambda.size(), 1);

    // ::pressio::utils::io::print_stdout("-----ynew-----\n");
    // ::pressio::utils::io::print_stdout(*y.data(), "\n");
    // ::pressio::utils::io::print_stdout("--------------\n");

    // ::pressio::utils::io::print_stdout("-----y2new-----\n");
    // ::pressio::utils::io::print_stdout(*y2.data(), "\n");
    // ::pressio::utils::io::print_stdout("--------------\n");

    // check convergence (whatever method user decided)
    const auto flag = is_conv_helper_t::evaluate(y, dy,
						 norm_dy, normRes, normRes0,
						 normJTRes, normJTRes0,
						 iStep, maxNonLIt, tolerance);
    if (flag) break;

    // store new norm into old variable
    normO = norm_dy;

    sys.residual(y, resid);
    sys.jacobian(y, jacob);

  }//loop

#if defined DEBUG_PRINT
  std::cout.precision(ss);
  ::pressio::utils::io::print_stdout(std::fixed);
#endif

#ifdef HAVE_TEUCHOS_TIMERS
  timer->stop("Gausss Newton Conserv");
#endif

}

}}}} //end namespace pressio::solvers::iterative::impl
#endif
