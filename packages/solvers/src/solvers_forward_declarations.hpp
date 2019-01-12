
#ifndef SOLVERS_FORWARD_DECLARATIONS_HPP_
#define SOLVERS_FORWARD_DECLARATIONS_HPP_

#include "solvers_ConfigDefs.hpp"

namespace rompp{ namespace solvers{ namespace iterative{

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
class GaussNewtonQR;


}}}//end namespace rompp::solvers::iterative

#endif
