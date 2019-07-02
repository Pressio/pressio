
#ifndef SOLVERS_IS_LEGITIMATE_SYSTEM_NONLIN_SOLVER_HPP_
#define SOLVERS_IS_LEGITIMATE_SYSTEM_NONLIN_SOLVER_HPP_

#include "solvers_basic_meta.hpp"
#include "solvers_system_has_all_needed_jacobian_methods.hpp"
#include "solvers_system_has_all_needed_residual_methods.hpp"
// #include "../../../containers/src/vector/containers_vector_meta.hpp"

namespace pressio{ namespace solvers{ namespace meta {

template<typename system_type, typename enable = void>
struct is_legitimate_system_for_nonlinear_solver : std::false_type{};

template<typename system_type>
struct is_legitimate_system_for_nonlinear_solver
<
  system_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::mpl::is_detected<has_scalar_typedef, system_type>::value   and
    ::pressio::mpl::is_detected<has_state_typedef, system_type>::value    and
    ::pressio::mpl::is_detected<has_residual_typedef, system_type>::value and
    ::pressio::mpl::is_detected<has_jacobian_typedef, system_type>::value and
    // containers::meta::is_vector_wrapper<
    //   typename system_type::state_type
    //   >::value and
    // containers::meta::is_vector_wrapper<
    //   typename system_type::residual_type
    //   >::value and
    system_has_needed_residual_methods<
      system_type,
      typename system_type::state_type,
      typename system_type::residual_type
      >::value and
    system_has_needed_jacobian_methods<
      system_type,
      typename system_type::state_type,
      typename system_type::jacobian_type
    >::value
    >
  > : std::true_type{};


}}} // namespace pressio::solvers::meta
#endif
