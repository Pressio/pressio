
#ifndef ODE_IS_LEGITIMATE_MODEL_FOR_IMPLICIT_ODE_HPP_
#define ODE_IS_LEGITIMATE_MODEL_FOR_IMPLICIT_ODE_HPP_

#include "../../../containers/src/meta/containers_meta_has_scalar_typedef.hpp"
#include "ode_has_state_typedef.hpp"
#include "ode_has_residual_typedef.hpp"
#include "ode_has_jacobian_typedef.hpp"
#include "ode_model_has_all_needed_residual_methods.hpp"
#include "ode_model_has_all_needed_jacobian_methods.hpp"

namespace rompp{ namespace ode{ namespace meta {

template<typename model_type,
	 typename enable = void>
struct is_legitimate_model_for_implicit_ode : std::false_type{};

template<typename model_type>
struct is_legitimate_model_for_implicit_ode<
  model_type,
  typename std::enable_if<
    ::rompp::containers::meta::has_scalar_typedef<model_type>::value and
    has_state_typedef<model_type>::value and
    has_residual_typedef<model_type>::value and
    has_jacobian_typedef<model_type>::value and
    // has residual methods
    model_has_needed_residual_methods<
      model_type,
      typename model_type::state_type,
      typename model_type::residual_type,
      typename model_type::scalar_type>::value and
    // has jacobian methods
    model_has_needed_jacobian_methods<
      model_type,
      typename model_type::state_type,
      typename model_type::jacobian_type,
      typename model_type::scalar_type>::value
    >::type
  > : std::true_type{};

}}} // namespace rompp::ode::meta
#endif
