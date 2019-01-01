
#ifndef SOLVERS_FORWARD_DECLARATIONS_HPP_
#define SOLVERS_FORWARD_DECLARATIONS_HPP_

#include "solvers_ConfigDefs.hpp"

namespace rompp{ namespace solvers{

template <
  typename scalar_t,
  typename qr_type,
  typename system_t = void,
  typename enable = void
  >
class GaussNewtonQR;


}}//end namespace rompp::solvers

#endif
