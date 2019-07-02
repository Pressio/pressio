
#ifndef ODE_IS_LEGITIMATE_IMPLICIT_RESIDUAL_TYPE_HPP_
#define ODE_IS_LEGITIMATE_IMPLICIT_RESIDUAL_TYPE_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../../../containers/src/vector/containers_vector_meta.hpp"

namespace pressio{ namespace ode{ namespace meta {

template<typename residual_type, typename enable = void>
struct is_legitimate_implicit_residual_type : std::false_type{};

template<typename residual_type>
struct is_legitimate_implicit_residual_type<residual_type,
 typename std::enable_if<
   containers::meta::is_vector_wrapper<residual_type>::value
#ifdef HAVE_PYBIND11
   or containers::meta::is_cstyle_array_pybind11<residual_type>::value
#endif
   >::type
  > : std::true_type{};

}}} // namespace pressio::ode::meta
#endif
