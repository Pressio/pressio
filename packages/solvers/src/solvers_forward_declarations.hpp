
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
  typename qr_type,
  typename line_search_t,
  typename when_converged_t = default_convergence,
  typename system_t = void,
  typename state_t = void,
  typename residual_t = void,
  typename jacobian_t = void,
  typename enable = void
  >
class GaussNewtonQR;

}//end namespace rompp::solvers::iterative::impl


/*
 * alias for the user for a GN solvers without line search
 */
template <
  typename scalar_t,
  typename qr_type,
  typename when_converged_t = default_convergence,
  typename system_t = void,
  typename state_t = void,
  typename residual_t = void,
  typename jacobian_t = void,
  typename enable = void
  >
using GaussNewtonQR = impl::GaussNewtonQR<
  scalar_t, qr_type, gn::noLineSearch,
  when_converged_t, system_t, state_t,
  residual_t, jacobian_t
  >;

/*
 * alias for the user for a GN solvers with line search
 */
template <
  typename scalar_t,
  typename qr_type,
  typename line_search_t,
  typename when_converged_t = default_convergence,
  typename system_t = void,
  typename state_t = void,
  typename residual_t = void,
  typename jacobian_t = void,
  typename enable = void
  >
using GaussNewtonQRLineSearch = impl::GaussNewtonQR<
  scalar_t, qr_type, line_search_t,
  when_converged_t, system_t, state_t,
  residual_t, jacobian_t
  >;


}}}//end namespace rompp::solvers::iterative

#endif
