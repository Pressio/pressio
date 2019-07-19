
#ifndef ODE_IS_LEGITIMATE_EXPLICIT_VELOCITY_TYPE_HPP_
#define ODE_IS_LEGITIMATE_EXPLICIT_VELOCITY_TYPE_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../../../containers/src/vector/containers_vector_meta.hpp"
#include "../../../containers/src/meta/containers_native_pybind_array_meta.hpp"

namespace pressio{ namespace ode{ namespace meta {

/*
 * for now, we enable if a vector is wrapped into containers::Vector
 * or if it is pybind::array
 */
template<typename T, typename enable = void>
struct is_legitimate_explicit_velocity_type : std::false_type{};

template<typename T>
struct is_legitimate_explicit_velocity_type<T,
 typename std::enable_if<
   containers::meta::is_vector_wrapper<T>::value
#ifdef HAVE_PYBIND11
   or containers::meta::is_cstyle_array_pybind11<T>::value
#endif
   >::type
  > : std::true_type{};


}}} // namespace pressio::ode::meta
#endif
