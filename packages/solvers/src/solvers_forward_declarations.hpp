
#ifndef SOLVERS_FORWARD_DECLARATIONS_HPP_
#define SOLVERS_FORWARD_DECLARATIONS_HPP_

#include "solvers_ConfigDefs.hpp"

namespace rompp{ namespace solvers{ namespace iterative{

namespace impl{

/*
 * this being inside impl should not be called by the user
 * it is only used for implementation
 */
template <
  typename scalar_t,
  typename lin_solver_tag,
  template <typename, typename> class lin_solver_t,
  typename line_search_t,
  typename when_converged_t = default_convergence,
  typename system_t	= void,
  typename state_t	= void,
  typename residual_t	= void,
  typename jacobian_t	= void,
  typename hessian_t	= void,
  typename observer_t	= void,
  typename enable	= void
  >
class GaussNewton;

/*
 * this being inside impl should not be called by the user
 * it is only used for implementation
 */
template <
  typename scalar_t,
  typename qr_type,
  typename line_search_t,
  typename when_converged_t = default_convergence,
  typename system_t	= void,
  typename state_t	= void,
  typename residual_t	= void,
  typename jacobian_t	= void,
  typename enable	= void
  >
class GaussNewtonQR;

}//end namespace rompp::solvers::iterative::impl


/* alias: GN solvers without line search */
template <
  typename scalar_t,
  typename lin_solver_tag,
  template <typename, typename> class lin_solver_t,
  typename when_converged_t = default_convergence,
  typename system_t	= void,
  typename hessian_t	= void,
  typename state_t	= void,
  typename residual_t	= void,
  typename jacobian_t	= void,
  typename observer_t	= void,
  typename enable	= void
  >
using GaussNewton = impl::GaussNewton
  <scalar_t, lin_solver_tag, lin_solver_t,
   gn::noLineSearch,
   when_converged_t, system_t, state_t,
   residual_t, jacobian_t,
   hessian_t, observer_t>;


/* alias: GN solvers with line search */
template <
  typename scalar_t,
  typename lin_solver_tag,
  template <typename, typename> class lin_solver_t,
  typename line_search_t,
  typename when_converged_t = default_convergence,
  typename system_t	= void,
  typename hessian_t	= void,
  typename state_t	= void,
  typename residual_t	= void,
  typename jacobian_t	= void,
  typename observer_t	= void,
  typename enable	= void
  >
using GaussNewtonLineSearch = impl::GaussNewton
  <scalar_t, lin_solver_tag, lin_solver_t,
   line_search_t, when_converged_t, system_t,
   state_t, residual_t, jacobian_t,
   hessian_t, observer_t>;


/* alias: QR-based GN solvers without line search */
template <
  typename scalar_t,
  typename qr_type,
  typename when_converged_t = default_convergence,
  typename system_t	= void,
  typename state_t	= void,
  typename residual_t	= void,
  typename jacobian_t	= void,
  typename enable	= void
  >
using GaussNewtonQR = impl::GaussNewtonQR
  <scalar_t, qr_type, gn::noLineSearch,
   when_converged_t, system_t, state_t,
   residual_t, jacobian_t>;


/* alias: QR-based GN solvers with line search */
template <
  typename scalar_t,
  typename qr_type,
  typename line_search_t,
  typename when_converged_t = default_convergence,
  typename system_t	= void,
  typename state_t	= void,
  typename residual_t	= void,
  typename jacobian_t	= void,
  typename enable	= void
  >
using GaussNewtonQRLineSearch = impl::GaussNewtonQR<
  scalar_t, qr_type, line_search_t,
  when_converged_t, system_t, state_t,
  residual_t, jacobian_t>;



namespace hacked{

/* solver for conservative ROM */
template <
  typename scalar_t,
  typename lin_solver_tag,
  template <typename, typename> class lin_solver_t,
  typename line_search_t,
  typename when_converged_t = default_convergence,
  typename system_t = void,
  typename cbar_t = void,
  typename enable = void
  >
class GaussNewtonConservative;

}//end namespace rompp::solvers::iterative::hacked

}}}//end namespace rompp::solvers::iterative

#endif
