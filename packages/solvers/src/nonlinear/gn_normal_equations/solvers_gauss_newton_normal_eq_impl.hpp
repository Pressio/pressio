
#ifndef SOLVERS_GAUSS_NEWTON_NORMAL_EQ_IMPL_HPP
#define SOLVERS_GAUSS_NEWTON_NORMAL_EQ_IMPL_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../helper_policies/solvers_converged_criterior_policy.hpp"
#include "../helper_policies/solvers_hessian_helper_policy.hpp"
#include "../helper_policies/solvers_jacob_res_product_policy.hpp"
#include "../helper_policies/solvers_norm_helper_policy.hpp"
#include "../helper_policies/solvers_line_search_policy.hpp"
#include "../helper_policies/solvers_residual_observer_when_solver_converged.hpp"
#include "../helper_policies/solvers_residual_observer_each_solver_step.hpp"
#include "../helper_policies/solvers_get_matrix_size_helper.hpp"

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{


template <
  typename line_search_t,
  typename converged_when_tag,
  typename system_t,
  typename hessian_t,
  typename lin_solver_t,
  typename iteration_t,
  typename scalar_t,
  typename observer_t = utils::impl::empty
  >
void gauss_newton_neq_solve(const system_t & sys,
			    typename system_t::state_type & y,
			    typename system_t::state_type & ytrial,
			    typename system_t::residual_type & resid,
			    typename system_t::jacobian_type & jacob,
			    typename system_t::state_type & dy,
			    typename system_t::state_type & JTR,
			    hessian_t & H,
			    lin_solver_t & linSolver,
			    iteration_t maxNonLIt,
			    scalar_t tolerance,
			    scalar_t & normO,
			    scalar_t & norm_dy,
			    const observer_t * observer,
			    std::string & convCondDescr)
{

  using residual_t	= typename system_t::residual_type;
  using jacobian_t	= typename system_t::jacobian_type;

  // find out which norm to use
  using norm_t = typename NormSelectorHelper<converged_when_tag>::norm_t;

  // policy for evaluating the norm of a vector
  using norm_evaluator_t = ComputeNormHelper<norm_t>;

  // policy to approximate hessian J^T*J
  using hessian_evaluator_t = HessianApproxHelper<jacobian_t>;

  // policy to J^T * residual
  using jtr_evaluator_t = JacobianTranspResProdHelper<jacobian_t>;

  // policy to checking convergence
  using is_converged_t = IsConvergedHelper<converged_when_tag>;

  // policy to observing residual at each GN step
  using residual_observer_each_step = ResidualObserverEachSolverStep<
    observer_t, residual_t>;

  // policy to observing residual when converged before exiting
  using residual_observer_when_conv = ResidualObserverWhenSolverConverged<
    observer_t, residual_t>;

  /* policy for computing line search factor (alpha) such that
   * the update is done with y = y + alpha dy
   * alpha = 1 default when user does not want line search
   */
  using lsearch_helper_t = LineSearchHelper<line_search_t>;

  //-------------------------------------------------------

  // alpha for taking steps
  scalar_t alpha = {};
  // storing residual norm
  scalar_t normRes = {};
  scalar_t normRes0 = {};
  // storing projected residual norm
  scalar_t normJTRes = {};
  scalar_t normJTRes0 = {};

  constexpr auto one = ::pressio::utils::constants::one<scalar_t>();
  convCondDescr = std::string(is_converged_t::description_);

#ifdef DEBUG_PRINT
  // get precision before GN
  auto ss = std::cout.precision();
  // set to 14 for the GN prints
  std::cout.precision(14);
  auto reset = utils::io::reset();
  auto fmt1 = utils::io::cyan() + utils::io::underline();
  ::pressio::utils::io::print_stdout(fmt1, "GN normal eqns:", "criterion:",
				     convCondDescr, reset, "\n");
#endif

  // compute the initial norm of y (the state)
  norm_evaluator_t::evaluate(y, normO);
  norm_dy = {0};

#ifdef HAVE_TEUCHOS_TIMERS
  auto timer = Teuchos::TimeMonitor::getStackedTimer();
  timer->start("NEQ-based Gausss Newton");
#endif

  iteration_t iStep = 0;
  while (++iStep <= maxNonLIt)
  {

    // call residual observer at each gauss step (no op for dummy case)
    residual_observer_each_step::evaluate(observer, resid, iStep);

#ifdef DEBUG_PRINT
    ::pressio::utils::io::print_stdout("\n");
    auto fmt = utils::io::underline();
    ::pressio::utils::io::print_stdout(fmt, "GN step", iStep,
				    utils::io::reset(), "\n");
#endif

    // residual norm for current state
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("norm resid");
#endif
    norm_evaluator_t::evaluate(resid, normRes);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("norm resid");
#endif

    // store initial residual norm
    if (iStep==1) normRes0 = normRes;

    // compute LHS: J^T*J
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("hessian");
#endif
    hessian_evaluator_t::evaluate(jacob, H);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("hessian");
#endif

    // ::pressio::utils::io::print_stdout("HESSIAN" , "\n");
    // ::pressio::utils::io::print_stdout(std::fixed,
    // 				    *H.data() , "\n");

#ifdef DEBUG_PRINT
    auto fmt1 = utils::io::magenta() + utils::io::bold();
    ::pressio::utils::io::print_stdout(fmt1, "GN_JSize =",
    ::pressio::solvers::impl::MatrixGetSizeHelper<jacobian_t>::globalRows(jacob),
    ::pressio::solvers::impl::MatrixGetSizeHelper<jacobian_t>::globalCols(jacob),
				    "\n");
    // this print only works when hessian is a shared mem matrix
    ::pressio::utils::io::print_stdout(fmt1, "GN_HessianSize =",
				    H.rows(), H.cols(),
				    utils::io::reset(), "\n");
#endif

    // compute RHS: J^T*res
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("JTR");
#endif
    jtr_evaluator_t::evaluate(jacob, resid, JTR);
    JTR.scale(static_cast<scalar_t>(-1));
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("JTR");
#endif

    // projected residual norm for current state
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("norm JTR");
#endif
    norm_evaluator_t::evaluate(JTR, normJTRes);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("norm JTR");
#endif

    // store initial residual norm
    if (iStep==1) normJTRes0 = normJTRes;

    // ::pressio::utils::io::print_stdout("J^T R \n");
    // ::pressio::utils::io::print_stdout( std::fixed,
    // 				     *JTR.data() , "\n");

    // solve normal equations
#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("solve normeq");
#endif
    linSolver.solveAllowMatOverwrite(H, JTR, dy);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("solve normeq");
#endif

    // // // print the correction
    // ::pressio::utils::io::print_stdout("Correction dy \n");
    // ::pressio::utils::io::print_stdout(std::fixed,
    // 				    *dy.data());

    // compute norm of the correction
    norm_evaluator_t::evaluate(dy, norm_dy);

#ifdef DEBUG_PRINT
    ::pressio::utils::io::print_stdout(std::scientific,
				    "||R|| =", normRes,
				    "||R||(r) =", normRes/normRes0,
				    "||J^T R|| =", normJTRes,
				    "||J^T R||(r) =", normJTRes/normJTRes0,
				    "||dy|| =", norm_dy,
				    utils::io::reset(),
				    "\n");
#endif

    // compute multiplicative factor if needed
    lsearch_helper_t::evaluate(alpha, y, ytrial, dy, resid, jacob, sys);

    // solution update: y = y + alpha*dy
    ::pressio::containers::ops::do_update(y, one, dy, alpha);

    // check convergence (whatever method user decided)
    const auto flag = is_converged_t::evaluate(y, dy,
					       norm_dy, normRes, normRes0,
					       normJTRes, normJTRes0,
					       iStep, maxNonLIt, tolerance);

    // if we have converged, query the observer
    if (flag) {
      // observe residual (no op for dummy case)
      residual_observer_when_conv::evaluate(observer, resid);
      break;
    }

    // store new norm into old variable
    normO = norm_dy;

    // compute residual and jacobian
    sys.residual(y, resid);
    sys.jacobian(y, jacob);

  }//loop

#if defined DEBUG_PRINT
  std::cout.precision(ss);
  ::pressio::utils::io::print_stdout(std::fixed);
#endif

#ifdef HAVE_TEUCHOS_TIMERS
  timer->stop("NEQ-based Gausss Newton");
#endif

}// end

}}}} //end namespace pressio::solvers::iterative::impl
#endif
