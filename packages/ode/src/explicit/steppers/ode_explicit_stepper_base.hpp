
#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_BASE_EXPLICIT_STEPPER_BASE_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_BASE_EXPLICIT_STEPPER_BASE_HPP_

#include "ode_explicit_stepper_traits.hpp"
#include "../policies/ode_is_legitimate_explicit_residual_policy.hpp"
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
  using step_traits	  = ode::details::traits<stepper_type>;
  using scalar_t	  = typename step_traits::scalar_t;
  using state_t		  = typename step_traits::state_t;
  using residual_t	  = typename step_traits::residual_t;
  using model_t		  = typename step_traits::model_t;
  using residual_policy_t = typename step_traits::residual_policy_t;

  static_assert( meta::is_legitimate_explicit_state_type<state_t>::value,
  "OOPS: STATE_TYPE IN SELECTED EXPLICIT STEPPER IS NOT VALID");

  static_assert( meta::is_legitimate_explicit_residual_type<residual_t>::value,
  "OOPS: RESIDUAL_TYPE IN SELECTED EXPLICIT STEPPER IS NOT VALID");

  static_assert( meta::is_legitimate_explicit_residual_policy<
  		 residual_policy_t>::value,
  "RESIDUAL_POLICY NOT ADMISSIBLE: MAYBE NOT INHERITING FROM EXPLICIT POLICY BASE");

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
