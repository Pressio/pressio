
#ifndef SOLVERS_IS_LEGITIMATE_RESIDUAL_OBSERVER_WHEN_CONVERGED_HPP
#define SOLVERS_IS_LEGITIMATE_RESIDUAL_OBSERVER_WHEN_CONVERGED_HPP

#include "../solvers_ConfigDefs.hpp"

namespace rompp{ namespace solvers{ namespace meta{

template <typename observer_t, typename residual_t, typename = void>
struct is_legitimate_residual_observer_when_solver_converged
  : std::false_type{};

template <typename observer_t, typename residual_t>
struct is_legitimate_residual_observer_when_solver_converged<
  observer_t, residual_t,
  ::rompp::mpl::void_t<
    decltype
    (
     std::declval<observer_t>().observeResidualWhenSolverConverged
     (
      std::declval<const residual_t &>()
      )
     )
    >
  > : std::true_type{};

}}} //end namespace rompp::solvers::meta
#endif
