
#ifndef ODE_IS_LEGITIMATE_EXPLICIT_RESIDUAL_TYPE_HPP_
#define ODE_IS_LEGITIMATE_EXPLICIT_RESIDUAL_TYPE_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../../../algebra/src/vector/algebra_vector_meta.hpp"
#include "../../../algebra/src/meta/algebra_native_pybind_array_meta.hpp"

namespace rompp{ namespace ode{ namespace meta {

/*
 * for now, we enable if a vector is wrapped into algebra::Vector
 * or if it is pybind::array
 */

template<typename residual_type, typename enable = void>
struct is_legitimate_explicit_residual_type : std::false_type{};

template<typename residual_type>
struct is_legitimate_explicit_residual_type<residual_type,
 typename std::enable_if<
   algebra::meta::is_algebra_vector_wrapper<residual_type>::value
#ifdef HAVE_PYBIND11
   or algebra::meta::is_cstyle_array_pybind11<residual_type>::value
#endif
   >::type
  > : std::true_type{};


}}} // namespace rompp::ode::meta
#endif
