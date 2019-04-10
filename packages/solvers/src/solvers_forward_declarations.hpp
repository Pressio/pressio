
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

template <
  typename system_t,
  typename qr_solver_t,
  typename scalar_t,
  typename line_search_t,
  typename when_converged_t,
  typename enable = void
  >
class GaussNewtonQR;

template <typename ... Args>
struct GNNEQSpecializationPicker;

template <typename ... Args>
struct GNQRSpecializationPicker;

}//end namespace rompp::solvers::iterative::impl


/* alias: GN solvers normal equations */
template <typename ... Args>
using GaussNewton = typename impl::GNNEQSpecializationPicker<Args...>::type;

/* alias: QR-based GN solvers */
template <typename ... Args>
using GaussNewtonQR = typename impl::GNQRSpecializationPicker<Args...>::type;


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
