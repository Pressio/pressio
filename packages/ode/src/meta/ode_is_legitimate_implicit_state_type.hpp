
#ifndef ODE_IS_LEGITIMATE_IMPLICIT_STATE_TYPE_HPP_
#define ODE_IS_LEGITIMATE_IMPLICIT_STATE_TYPE_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../../../containers/src/vector/containers_vector_meta.hpp"

namespace pressio{ namespace ode{ namespace meta {

template<typename T, typename enable = void>
struct is_legitimate_implicit_state_type : std::false_type{};

template<typename T>
struct is_legitimate_implicit_state_type<
  T,
  typename std::enable_if<
    containers::meta::is_vector_wrapper<T>::value
#ifdef HAVE_PYBIND11
    or containers::meta::is_cstyle_array_pybind11<T>::value
#endif
    >::type > : std::true_type{};

}}} // namespace pressio::ode::meta
#endif
