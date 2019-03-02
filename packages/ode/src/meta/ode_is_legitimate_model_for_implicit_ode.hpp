
#ifndef ODE_IS_LEGITIMATE_MODEL_FOR_IMPLICIT_ODE_HPP_
#define ODE_IS_LEGITIMATE_MODEL_FOR_IMPLICIT_ODE_HPP_

#include "ode_basic_meta.hpp"
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
   // has to have scalar typedef
   core::meta::is_detected<has_scalar_typedef, model_type>::value and
   // has to have state typedef
   core::meta::is_detected<has_state_typedef, model_type>::value and
   // has to have residual typedef
   core::meta::is_detected<has_residual_typedef, model_type>::value and
   // has to have jacobian typedef
   core::meta::is_detected<has_jacobian_typedef, model_type>::value and
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


//------------------------------------------------------

// for now, just leave Backward Euler as the legitimate stepper
template<typename T,
	 typename enable = void>
struct is_legitimate_auxiliary_stepper : std::false_type{};

template<typename T>
struct is_legitimate_auxiliary_stepper<T,
  typename std::enable_if<
    ::rompp::ode::details::traits<typename T::base_t>::enum_id ==
    ::rompp::ode::ImplicitEnum::Euler
    >::type
  > : std::true_type{};



}}} // namespace rompp::ode::meta
#endif
