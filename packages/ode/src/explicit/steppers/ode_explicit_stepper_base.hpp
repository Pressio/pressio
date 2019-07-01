
#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_BASE_EXPLICIT_STEPPER_BASE_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_BASE_EXPLICIT_STEPPER_BASE_HPP_

#include "ode_explicit_stepper_traits.hpp"
#include "../policies/ode_is_legitimate_explicit_velocity_policy.hpp"
#include "../../ode_storage.hpp"
#include "../../ode_system_wrapper.hpp"

namespace rompp{ namespace ode{

/*
 * (1) constructors here should be private but we need
 * them public to enable interfacing with pybind11
 */

template<typename stepper_type>
class ExplicitStepperBase
{
private:
  using step_traits	= ode::details::traits<stepper_type>;
  using scalar_t	= typename step_traits::scalar_t;
  using state_t		= typename step_traits::state_t;
  using velocity_t	= typename step_traits::velocity_t;
  using model_t		= typename step_traits::model_t;
  using policy_t	= typename step_traits::velocity_policy_t;

  static_assert( meta::is_legitimate_explicit_state_type<state_t>::value,
  "OOPS: STATE_TYPE IN SELECTED EXPLICIT STEPPER IS NOT VALID");

  static_assert( meta::is_legitimate_explicit_velocity_type<velocity_t>::value,
  "OOPS: VELOCITY_TYPE IN SELECTED EXPLICIT STEPPER IS NOT VALID");

  static_assert( meta::is_legitimate_explicit_velocity_policy<
  		 policy_t>::value,
  "VELOCITY_POLICY NOT ADMISSIBLE: MAYBE NOT INHERITING FROM EXPLICIT POLICY BASE");

public:
  ExplicitStepperBase() = default;
  ~ExplicitStepperBase() = default;

public:
  typename step_traits::order_t order() const{
    return step_traits::order_value;
  }

  template<typename step_t>
  void operator()(state_t & yinout,
		  scalar_t t,
		  scalar_t dt,
		  step_t step){
    static_cast<stepper_type&>(*this).compute(yinout, t, dt, step);
  }

};//end class

}}//end namespace rompp::ode
#endif
