
#ifndef SOLVERS_GAUSS_NEWTON_QR_IMPL_HPP
#define SOLVERS_GAUSS_NEWTON_QR_IMPL_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../solvers_system_traits.hpp"
#include "../../solvers_meta_static_checks.hpp"
#include "../mixins/solvers_converged_criterior_mixin.hpp"
#include "../mixins/solvers_norm_helper_mixin.hpp"
#include "../mixins/solvers_line_search_mixin.hpp"


namespace rompp{ namespace solvers{ namespace iterative{ namespace impl{

template <
  typename system_t,
  typename iteration_t,
  typename scalar_t,
  typename qr_obj_t,
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
void gauss_newtom_qr_solve(const system_t & sys,
			   typename system_t::state_type & y,
			   typename system_t::state_type & ytrial,
			   typename system_t::residual_type & resid,
			   typename system_t::jacobian_type & jacob,
			   iteration_t maxNonLIt,
			   scalar_t tolerance,
			   typename system_t::state_type & QTResid,
			   typename system_t::state_type & dy,
			   qr_obj_t & qrObj,
			   scalar_t & normO,
			   scalar_t & normN){

  // find out which norm to use
  using norm_t = typename NormSelectorHelper<converged_when_tag>::norm_t;

  // functor for evaluating the norm of a vector
  using norm_evaluator_t = ComputeNormHelper<norm_t>;
  norm_evaluator_t normEvaluator;

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

#ifdef DEBUG_PRINT
  // get precision before GN
  auto ss = std::cout.precision();
  // set to 14 for the GN prints
  std::cout.precision(14);

  auto reset = core::io::reset();
  auto fmt1 = core::io::cyan() + core::io::underline();
  const auto convString = std::string(is_conv_helper_t::description_);
  ::rompp::core::io::print_stdout(fmt1, "GN with QR:", "criterion:",
				  convString, reset, "\n");
#endif

  // compute (whatever type) norm of y
  normEvaluator(y, normO);
  normN = static_cast<scalar_t>(0);

#ifdef HAVE_TEUCHOS_TIMERS
  auto timer = Teuchos::TimeMonitor::getStackedTimer();
  timer->start("QR-based Gausss Newton");
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
    normEvaluator(resid, normRes);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("QR factorization");
    qrObj.computeThin(jacob);
    timer->stop("QR factorization");
#else
    qrObj.computeThin(jacob);
#endif

    // compute: Q^T Residual
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("QR projection");
    qrObj.project(resid, QTResid);
    timer->stop("QR projection");
#else
    qrObj.project(resid, QTResid);
#endif

    // compute correction: dy
    // by solving R dy = - Q^T Residual
    QTResid.scale(static_cast<scalar_t>(-1));
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("QR solve");
    qrObj.solve(QTResid, dy);
    timer->stop("QR solve");
#else
    qrObj.solve(QTResid, dy);
#endif

    // norm of the correction
    normEvaluator(dy, normN);

#ifdef DEBUG_PRINT
    ::rompp::core::io::print_stdout("norm(residual) =",
				    std::setprecision(14),
				    normRes, ",",
				    "norm(correction) =", normN,
				    "\n");
#endif

    // after correction is computed,
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

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("residual");
    sys.residual(y, resid);
    timer->stop("residual");
#else
    sys.residual(y, resid);
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("jacobian");
    sys.jacobian(y, jacob);
    timer->stop("jacobian");
#else
    sys.jacobian(y, jacob);
#endif
  }

#if defined DEBUG_PRINT
  std::cout.precision(ss);
#endif

#ifdef HAVE_TEUCHOS_TIMERS
  timer->stop("QR-based Gausss Newton");
#endif

}//


}}}} //end namespace rompp::solvers::iterative::implo
#endif
