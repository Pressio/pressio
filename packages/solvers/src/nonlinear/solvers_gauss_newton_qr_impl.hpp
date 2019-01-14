
#ifndef SOLVERS_GAUSS_NEWTON_QR_IMPL_HPP
#define SOLVERS_GAUSS_NEWTON_QR_IMPL_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../solvers_system_traits.hpp"
#include "../solvers_meta_static_checks.hpp"
#include "../../../CORE_OPS"
#include "../../../QR_BASIC"

namespace rompp{ namespace solvers{ namespace iterative{ namespace impl{

template <typename convergence_tag>
struct norm_type_to_use{
  using norm_t = L2Norm;
};

template <typename norm_type>
struct norm_type_to_use<
  converged_when::absoluteNormCorrectionBelowTol<norm_type>
  >{
  using norm_t = norm_type;
};
//---------------------------------------------------------

template <typename norm_t>
struct computeNormHelper;

template <>
struct computeNormHelper<::rompp::solvers::L2Norm>{
  template <typename vec_t, typename scalar_t>
  void operator()(const vec_t & vecIn, scalar_t & result){
    result = ::rompp::core::ops::norm2(vecIn);
  }
};

template <>
struct computeNormHelper<::rompp::solvers::L1Norm>{
  template <typename vec_t, typename scalar_t>
  void operator()(const vec_t & vecIn, scalar_t & result){
    result = ::rompp::core::ops::norm1(vecIn);
  }
};
//---------------------------------------------------------

template <typename conv_tag>
struct isConvergedHelper;

template <>
struct isConvergedHelper<converged_when::completingNumMaxIters>{
  template <typename state_t, typename step_t, typename scalar_t>
  bool operator()(const state_t & x, const state_t & dx,
		  scalar_t norm_dx, step_t step,
		  step_t maxIters, scalar_t tol){
    return step==maxIters;
  }
};

template <typename norm_t>
struct isConvergedHelper<
  converged_when::absoluteNormCorrectionBelowTol<norm_t>
  >{
  template <typename state_t, typename step_t, typename scalar_t>
  bool operator()(const state_t & x, const state_t & dx,
		  scalar_t norm_dx, step_t step,
		  step_t maxIters, scalar_t tol){
    return (norm_dx<tol);
  }
};
//---------------------------------------------------------



template <
  typename system_t,
  typename iteration_t,
  typename scalar_t,
  typename qr_obj_t,
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
			   typename system_t::state_type & x,
			   typename system_t::residual_type & resid,
			   typename system_t::jacobian_type & jacob,
			   iteration_t maxNonLIt,
			   scalar_t tolerance,
			   typename system_t::state_type & QTResid,
			   typename system_t::state_type & dx,
			   qr_obj_t & qrObj,
			   scalar_t & normO,
			   scalar_t & normN)
{

  // find out which norm to use
  using norm_t = typename norm_type_to_use<converged_when_tag>::norm_t;

  // functor for evaluating the norm of a vector
  using norm_evaluator_t = computeNormHelper<norm_t>;
  norm_evaluator_t normEvaluator;

  // functor for checking convergence
  using is_conv_helper_t = isConvergedHelper<converged_when_tag>;
  is_conv_helper_t isConverged;

  // storing residaul norm
  scalar_t normRes = {};

#ifdef DEBUG_PRINT
  ::rompp::core::debug::print("starting GN solve with",
			      "tol=",tolerance,",",
			      "maxIter=", maxNonLIt,"\n");
#endif

  // compute (whatever type) norm of x
  normEvaluator(x, normO);
  normN = static_cast<scalar_t>(0);

#ifdef HAVE_TEUCHOS_TIMERS
  auto timer = Teuchos::TimeMonitor::getStackedTimer();
  timer->start("QR-based Gausss Newton");
#endif

  iteration_t iStep = 0;
  while (iStep++ <= maxNonLIt)
    {

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

      // solve R dx = Q^T Residual
#ifdef HAVE_TEUCHOS_TIMERS
      timer->start("QR solve");
      qrObj.solve(QTResid, dx);
      timer->stop("QR solve");
#else
      qrObj.solve(QTResid, dx);
#endif

      normEvaluator(resid, normRes);
      normEvaluator(dx, normN);
#ifdef DEBUG_PRINT
      ::rompp::core::debug::print("GN step", iStep, ",",
				  "norm(dx)=", normN, ",",
				  "norm(residual)=", normRes, "\n");
#endif

      // update solution
      x -= dx;

      // check convergence (based on whatever method user decided)
      auto flag = isConverged(x, dx, normN, iStep, maxNonLIt, tolerance);
      if (flag) break;

      normO = normN;

#ifdef HAVE_TEUCHOS_TIMERS
      timer->start("residual");
      sys.residual(x, resid);
      timer->stop("residual");
#else
      sys.residual(x, resid);
#endif

#ifdef HAVE_TEUCHOS_TIMERS
      timer->start("jacobian");
      sys.jacobian(x, jacob);
      timer->stop("jacobian");
#else
      sys.jacobian(x, jacob);
#endif
    }

#ifdef HAVE_TEUCHOS_TIMERS
  timer->stop("QR-based Gausss Newton");
#endif

}//


}}}} //end namespace rompp::solvers::iterative::implo
#endif
