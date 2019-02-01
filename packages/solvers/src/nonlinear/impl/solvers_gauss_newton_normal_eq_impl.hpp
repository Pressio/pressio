
#ifndef SOLVERS_GAUSS_NEWTON_NORMAL_EQ_IMPL_HPP
#define SOLVERS_GAUSS_NEWTON_NORMAL_EQ_IMPL_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../solvers_system_traits.hpp"
#include "../../solvers_meta_static_checks.hpp"
#include "../helper_policies/solvers_converged_criterior_policy.hpp"
#include "../helper_policies/solvers_hessian_helper_policy.hpp"
#include "../helper_policies/solvers_jacob_res_product_policy.hpp"
#include "../helper_policies/solvers_norm_helper_policy.hpp"
#include "../helper_policies/solvers_line_search_policy.hpp"

namespace rompp{ namespace solvers{ namespace iterative{ namespace impl{

template <
  typename system_t,
  typename hessian_t,
  typename iteration_t,
  typename scalar_t,
  typename lin_solver_t,
  typename line_search_t,
  typename converged_when_tag,
  core::meta::enable_if_t<
    ::rompp::solvers::details::system_traits<system_t>::is_system and
    core::meta::is_core_vector_wrapper<
      typename system_t::state_type>::value and
    core::meta::is_core_vector_wrapper<
      typename system_t::residual_type>::value
    > * =nullptr
  >
void gauss_newtom_neq_solve(const system_t & sys,
			    typename system_t::state_type & y,
			    typename system_t::state_type & ytrial,
			    typename system_t::residual_type & resid,
			    typename system_t::jacobian_type & jacob,
			    hessian_t & H,
			    typename system_t::state_type & JTR,
			    iteration_t maxNonLIt,
			    scalar_t tolerance,
			    typename system_t::state_type & dy,
			    lin_solver_t & linSolver,
			    scalar_t & normO,
			    scalar_t & normN){

  using jac_t	= typename system_t::jacobian_type;

  // find out which norm to use
  using norm_t = typename NormSelectorHelper<converged_when_tag>::norm_t;

  // functor for evaluating the norm of a vector
  using norm_evaluator_t = ComputeNormHelper<norm_t>;
  norm_evaluator_t normEvaluator;

  // functor for approximate hessian J^T*J
  using hessian_evaluator_t = HessianApproxHelper<jac_t>;
  hessian_evaluator_t hessEvaluator;

  // functor for J^T * residual
  using jtr_evaluator_t = JacobianTranspResProdHelper<jac_t>;
  jtr_evaluator_t jtrEvaluator;

  // functor for checking convergence
  using is_conv_helper_t = IsConvergedHelper<converged_when_tag>;
  is_conv_helper_t isConverged;

  /* functor for computing line search factor (alpha) such that
   * the update is done with y = y + alpha dy
   * alpha = 1 default when user does not want line search
   */
  using lsearch_helper = LineSearchHelper<line_search_t>;
  lsearch_helper lineSearchHelper;
  //-------------------------------------------------------
  // alpha for taking steps
  scalar_t alpha = {};
  // storing residaul norm
  scalar_t normRes = {};
  scalar_t normRes0 = {};

#ifdef DEBUG_PRINT
  // get precision before GN
  auto ss = std::cout.precision();
  // set to 14 for the GN prints
  std::cout.precision(14);

  auto reset = core::io::reset();
  auto fmt1 = core::io::cyan() + core::io::underline();
  const auto convString = std::string(is_conv_helper_t::description_);
  ::rompp::core::io::print_stdout(fmt1, "GN normal eqns:", "criterion:",
				  convString, reset, "\n");
#endif

  // compute (whatever type) norm of y
  normEvaluator(y, normO);
  normN = static_cast<scalar_t>(0);

#ifdef HAVE_TEUCHOS_TIMERS
  auto timer = Teuchos::TimeMonitor::getStackedTimer();
  timer->start("NEQ-based Gausss Newton");
#endif

  iteration_t iStep = 0;
  while (iStep++ <= maxNonLIt)
  {

#ifdef DEBUG_PRINT
    ::rompp::core::io::print_stdout("\n");
    auto fmt = core::io::underline();
    ::rompp::core::io::print_stdout(fmt, "GN step", iStep,
				    core::io::reset(), "\n");
#endif

    // residual norm for current state
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("norm resid");
    normEvaluator(resid, normRes);
    timer->stop("norm resid");
#else
    normEvaluator(resid, normRes);
#endif
    // store initial residual norm
    if (iStep==1) normRes0 = normRes;

    // compute LHS: J^T*J
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("hessian");
    hessEvaluator(jacob, H);
    timer->stop("hessian");
#else
    hessEvaluator(jacob, H);
#endif

    // compute RHS: J^T*res
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("JTR");
    jtrEvaluator(jacob, resid, JTR);
    JTR.scale(static_cast<scalar_t>(-1));
    timer->stop("JTR");
#else
    jtrEvaluator(jacob, resid, JTR);
    JTR.scale(static_cast<scalar_t>(-1));
#endif

    // solve normal equations
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("solve normeq");
    linSolver.solve(H, JTR, dy);
    timer->stop("solve normeq");
#else
    linSolver.solve(H, JTR, dy);
#endif

    // norm of the correction
    normEvaluator(dy, normN);

#ifdef DEBUG_PRINT
    ::rompp::core::io::print_stdout(std::scientific,
				    "||R|| =", normRes,
				    "||R||(r) =", normRes/normRes0,
				    "||dy|| =", normN,
				    "\n");
#endif

    // compute multiplicative factor if needed
    lineSearchHelper(alpha, y, ytrial, dy, resid, jacob, sys);

    // solution update
    y = y + alpha*dy;

    // check convergence (whatever method user decided)
    auto flag = isConverged(y, dy, normN, iStep,
			    maxNonLIt, tolerance);
    if (flag) break;

    // store new norm into old variable
    normO = normN;

    sys.residual(y, resid);
    sys.jacobian(y, jacob);

  }//loop

#if defined DEBUG_PRINT
  std::cout.precision(ss);
#endif

#ifdef HAVE_TEUCHOS_TIMERS
  timer->stop("NEQ-based Gausss Newton");
#endif

}

}}}} //end namespace rompp::solvers::iterative::implo
#endif
