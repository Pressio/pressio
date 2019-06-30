
#ifndef SOLVERS_IS_LEGITIMATE_RESIDUAL_OBSERVER_EACH_SOLVER_STEP_HPP
#define SOLVERS_IS_LEGITIMATE_RESIDUAL_OBSERVER_EACH_SOLVER_STEP_HPP

#include "../solvers_ConfigDefs.hpp"

namespace rompp{ namespace solvers{ namespace meta{

template <typename observer_t, typename residual_t, typename = void>
struct is_legitimate_residual_observer_each_solver_step
  : std::false_type{};

template <typename observer_t, typename residual_t>
struct is_legitimate_residual_observer_each_solver_step<
  observer_t, residual_t,
  ::rompp::mpl::void_t<
    decltype
    (
     std::declval<observer_t>().observeResidualEachGNStep
     (
      std::declval<int>(), // this int is for passing the solver step
      std::declval<const residual_t &>()
      )
     )
    >
  > : std::true_type{};

}}} //end namespace rompp::solvers::meta
#endif
