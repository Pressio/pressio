
#ifndef SOLVERS_GAUSS_NEWTON_NORMAL_EQ_CONSERVATIVE_IMPL_HPP
#define SOLVERS_GAUSS_NEWTON_NORMAL_EQ_CONSERVATIVE_IMPL_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../solvers_system_traits.hpp"
#include "../../solvers_meta_static_checks.hpp"
#include "../helper_policies/solvers_converged_criterior_policy.hpp"
#include "../helper_policies/solvers_hessian_helper_policy.hpp"
#include "../helper_policies/solvers_jacob_res_product_policy.hpp"
#include "../helper_policies/solvers_norm_helper_policy.hpp"
#include "../helper_policies/solvers_line_search_policy.hpp"
#include "solvers_get_matrix_size_helper.hpp"

namespace rompp{ namespace solvers{ namespace iterative{ namespace impl{

template <
  typename system_t,
  typename iteration_t,
  typename scalar_t,
  typename lin_solver_t,
  typename line_search_t,
  typename converged_when_tag,
  typename cbar_t,
  typename mat_t,
  core::meta::enable_if_t<
    ::rompp::solvers::details::system_traits<system_t>::is_system and
    core::meta::is_vector_wrapper_eigen<
      typename system_t::state_type>::value and
    core::meta::is_core_vector_wrapper<
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
				    scalar_t & normN,
				    const cbar_t & cbarT,
				    mat_t & jTj_,
				    mat_t & jTcbarT_,
				    mat_t & cbarJ_,
				    mat_t & zero_,
				    typename system_t::state_type & cbarTlambda_,
				    typename system_t::state_type & jTr2_,
				    typename system_t::state_type & cbarR_,
				    mat_t & A_,
				    typename system_t::state_type & b_,
				    typename system_t::state_type & lambda_,
				    typename system_t::state_type & y2_)
{

  using jac_t	= typename system_t::jacobian_type;

  // find out which norm to use
  using norm_t = typename NormSelectorHelper<converged_when_tag>::norm_t;

  // functor for evaluating the norm of a vector
  using norm_evaluator_t = ComputeNormHelper<norm_t>;
  norm_evaluator_t normEvaluator;

  // // functor for approximate hessian J^T*J
  // using hessian_evaluator_t = HessianApproxHelper<jac_t>;
  // hessian_evaluator_t hessEvaluator;

  // // functor for J^T * residual
  // using jtr_evaluator_t = JacobianTranspResProdHelper<jac_t>;
  // jtr_evaluator_t jtrEvaluator;

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
  // get precision
  auto ss = std::cout.precision();
  // set to 14 for prints
  std::cout.precision(14);

  auto reset = core::io::reset();
  auto fmt1 = core::io::cyan() + core::io::underline();
  const auto convString = std::string(is_conv_helper_t::description_);
  ::rompp::core::io::print_stdout(fmt1, "GN normal eqns conserv:",
				  "criterion:",
				  convString, reset, "\n");
#endif

  // compute (whatever type) norm of y
  normEvaluator(y, normO);
  normN = static_cast<scalar_t>(0);

#ifdef HAVE_TEUCHOS_TIMERS
  auto timer = Teuchos::TimeMonitor::getStackedTimer();
  timer->start("Gausss Newton Conserv");
#endif

  iteration_t iStep = 0;
  while (iStep++ <= maxNonLIt)
  {

#ifdef DEBUG_PRINT
    ::rompp::core::io::print_stdout("\n");
    auto fmt = core::io::underline();
    ::rompp::core::io::print_stdout(fmt, "step", iStep,
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

    // assemble LHS
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("lhs");
#endif
    ::rompp::core::ops::dot_self(jacob, jTj_);
    ::rompp::core::ops::dot(jacob, cbarT, jTcbarT_);
    ::rompp::core::ops::dot(cbarT, jacob, cbarJ_);

    A_.block(0, 0, jTj_.rows(), jTj_.cols()) = jTj_;
    A_.block(0, jTj_.cols(), jTcbarT_.rows(), jTcbarT_.cols()) = jTcbarT_;

    A_.block(jTj_.rows(), 0, cbarJ_.rows(), cbarJ_.cols()) = cbarJ_;
    A_.block(jTj_.rows(), jTj_.cols(), zero_.rows(), zero_.cols()) = zero_;

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("lhs");
#endif

#ifdef DEBUG_PRINT
    // auto fmt1 = core::io::magenta() + core::io::bold();
    // ::rompp::core::io::print_stdout(fmt1, "GN_JSize =",
    // ::rompp::solvers::impl::MatrixGetSizeHelper<jac_t>::globalRows(jacob),
    // ::rompp::solvers::impl::MatrixGetSizeHelper<jac_t>::globalCols(jacob),
    // 				    "\n");
    // // this print only works when hessian is a shared mem matrix
    // ::rompp::core::io::print_stdout(fmt1, "GN_HessianSize =",
    // 				    H.rows(), H.cols(),
    // 				    core::io::reset(), "\n");
#endif

    // compute RHS
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("rhs");
#endif
    ::rompp::core::ops::dot(cbarT_, resid, jTr2_);
    ::rompp::core::ops::product(cbarT, lambda_, cbarTlambda_);
    resid += cbarTlambda_;
    ::rompp::core::ops::dot(jacob_, resid, jTr2_);

    b_.block(0, 0, jTr2_.rows(), 1) = -jTr2_;
    b_.block(jTr2_.rows(), 0, cbarR_.rows(), 1) = -cbarR_;

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("rhs");
#endif


    // solve normal equations
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("solve normeq");
#endif
    y2_.data()->block(0, 0, y.size(), 1) = *y.data();
    y2_.data()->block(y.size(), 0, lambda.size(), 1) = *lambda.data();

    linSolver.solve(A_, b_, dy);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("solve normeq");
#endif

    // norm of the correction
    normEvaluator(dy, normN);

#ifdef DEBUG_PRINT
    ::rompp::core::io::print_stdout(std::scientific,
				    "||R|| =", normRes,
				    "||R||(r) =", normRes/normRes0,
				    "||dy|| =", normN,
				    core::io::reset(),
				    "\n");
#endif

    // // compute multiplicative factor if needed
    // lineSearchHelper(alpha, y, ytrial, dy, resid, jacob, sys);

    // solution update
    y = y + alpha*dy.data()->block(0, 0, y.size(), 1);
    lambda = lambda + alpha*dy.data()->block(y.size(), 0, lambda.size(), 1);

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
  ::rompp::core::io::print_stdout(std::fixed);
#endif

#ifdef HAVE_TEUCHOS_TIMERS
  timer->stop("Gausss Newton Conserv");
#endif

}

}}}} //end namespace rompp::solvers::iterative::impl
#endif
