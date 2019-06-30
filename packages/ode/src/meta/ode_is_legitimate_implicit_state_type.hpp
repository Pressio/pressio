
#ifndef ODE_IS_LEGITIMATE_IMPLICIT_STATE_TYPE_HPP_
#define ODE_IS_LEGITIMATE_IMPLICIT_STATE_TYPE_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../../../containers/src/vector/containers_vector_meta.hpp"

namespace rompp{ namespace ode{ namespace meta {

template<typename state_type, typename enable = void>
struct is_legitimate_implicit_state_type : std::false_type{};

template<typename state_type>
struct is_legitimate_implicit_state_type<
  state_type,
  typename std::enable_if<
    containers::meta::is_vector_wrapper<state_type>::value
#ifdef HAVE_PYBIND11
    or containers::meta::is_cstyle_array_pybind11<state_type>::value
#endif
    >::type > : std::true_type{};

}}} // namespace rompp::ode::meta
#endif
