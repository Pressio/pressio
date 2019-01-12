
#ifndef SOLVERS_GAUSS_NEWTON_QR_IMPL_HPP
#define SOLVERS_GAUSS_NEWTON_QR_IMPL_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../solvers_system_traits.hpp"
#include "../solvers_meta_static_checks.hpp"
#include "../../../CORE_OPS"
#include "../../../QR_BASIC"

namespace rompp{ namespace solvers{ namespace iterative{ namespace impl{

template <
  typename system_t,
  typename iteration_t,
  typename scalar_t,
  typename qr_obj_t,
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

#ifdef DEBUG_PRINT
  int myRank = 0;

#ifdef HAVE_MPI
  int flag = 0;
  MPI_Initialized( &flag );
  if (flag==1)
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#endif

  if (myRank==0)
    std::cout << " starting Gauss-Newton solve "
	      << " tol = " << tolerance
	      << " maxIter = " << maxNonLIt
	      << std::endl;
#endif
  // here we should set dx=x, and do normO = norm2(dx)
  // but we can just put norm2(x) so we avoid an assignment operation
  normO = ::rompp::core::ops::norm2(x);
  normN = static_cast<scalar_t>(0);

#ifdef HAVE_TEUCHOS_TIMERS
  auto timer = Teuchos::TimeMonitor::getStackedTimer();
  timer->start("QR-based Gausss Newton");
#endif

  iteration_t iStep = 1;
  while (iStep++ < maxNonLIt)
    {

#ifdef HAVE_TEUCHOS_TIMERS
      timer->start("QR factorization");
      // QR decomposition of Jacobian
      qrObj.computeThin(jacob);
      timer->stop("QR factorization");
#else
      // QR decomposition of Jacobian
      qrObj.computeThin(jacob);
#endif


#ifdef HAVE_TEUCHOS_TIMERS
      timer->start("QR projection");
      // compute: Q^T Residual
      qrObj.project(resid, QTResid);
      timer->stop("QR projection");
#else
      // compute: Q^T Residual
      qrObj.project(resid, QTResid);
#endif


#ifdef HAVE_TEUCHOS_TIMERS
      timer->start("QR solve");
      // solve R dx = Q^T Residual
      qrObj.solve(QTResid, dx);
      timer->stop("QR solve");
#else
      // solve R dx = Q^T Residual
      qrObj.solve(QTResid, dx);
#endif

      // update solution
      x -= dx;

      normN = ::rompp::core::ops::norm2(dx);
#ifdef DEBUG_PRINT
      if (myRank==0)
	std::cout << " GN step=" << iStep
		  << " norm(dx)= " << normN
		  << std::endl;
#endif

      if (std::abs(normO - normN) < tolerance){
#ifdef DEBUG_PRINT
	if (myRank==0)
	  std::cout << " GN converged! "
		    << " final norm(dx)= " << normN
		    << std::endl;
#endif
      	break;
      }

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


}}}} //end namespace rompp::solvers::iterative::impl
#endif
