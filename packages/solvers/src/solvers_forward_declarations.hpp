
#ifndef SOLVERS_FORWARD_DECLARATIONS_HPP_
#define SOLVERS_FORWARD_DECLARATIONS_HPP_

#include "solvers_ConfigDefs.hpp"
#include "solvers_convergence_tags.hpp"
#include "solvers_line_search_tags.hpp"

namespace rompp{ namespace solvers{ namespace iterative{

/* impl should not be called by the user */
namespace impl{

template <
  typename system_t,
  typename hessian_t,
  typename linear_solver_t,
  typename scalar_t,
  typename line_search_t,
  typename when_converged_t,
  typename resid_obs_t,
  typename enable = void
  >
class GaussNewton;

template <typename ... Args>
struct GNPicker;

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
template <typename ... Args>
using GaussNewton = typename impl::GNPicker<Args...>::type;


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
