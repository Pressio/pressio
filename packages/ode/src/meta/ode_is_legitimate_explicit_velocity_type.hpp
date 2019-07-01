
#ifndef ODE_IS_LEGITIMATE_EXPLICIT_VELOCITY_TYPE_HPP_
#define ODE_IS_LEGITIMATE_EXPLICIT_VELOCITY_TYPE_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../../../containers/src/vector/containers_vector_meta.hpp"
#include "../../../containers/src/meta/containers_native_pybind_array_meta.hpp"

namespace rompp{ namespace ode{ namespace meta {

/*
 * for now, we enable if a vector is wrapped into containers::Vector
 * or if it is pybind::array
 */
template<typename residual_type, typename enable = void>
struct is_legitimate_explicit_velocity_type : std::false_type{};

template<typename residual_type>
struct is_legitimate_explicit_velocity_type<residual_type,
 typename std::enable_if<
   containers::meta::is_vector_wrapper<residual_type>::value
#ifdef HAVE_PYBIND11
   or containers::meta::is_cstyle_array_pybind11<residual_type>::value
#endif
   >::type
  > : std::true_type{};


}}} // namespace rompp::ode::meta
#endif
