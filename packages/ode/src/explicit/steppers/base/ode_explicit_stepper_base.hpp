
#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_BASE_EXPLICIT_STEPPER_BASE_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_BASE_EXPLICIT_STEPPER_BASE_HPP_

#include "../ode_explicit_stepper_traits.hpp"
#include "../../policies/ode_explicit_policies_meta.hpp"
#include "../../../ode_storage.hpp"
#include "../../../ode_aux_data.hpp"
#include "../../../ode_basic_meta.hpp"

namespace rompp{ namespace ode{

template<typename stepper_type>
class ExplicitStepperBase
  : private core::details::CrtpBase<
  ExplicitStepperBase<stepper_type>>
{
private:
  using step_traits = ode::details::traits<stepper_type>;
  using scalar_t = typename step_traits::scalar_t;
  using state_t = typename step_traits::state_t; 
  using residual_t = typename step_traits::residual_t;
  using model_t = typename step_traits::model_t;
  using residual_policy_t = typename step_traits::residual_policy_t;

  using order_t = typename step_traits::order_t;
  static constexpr order_t order_value = step_traits::order_value;

  //do checking here that things are as supposed
  static_assert( meta::is_legitimate_explicit_state_type<state_t>::value,
  "OOPS: STATE_TYPE IN SELECTED EXPLICIT STEPPER IS NOT VALID");
  static_assert( meta::is_legitimate_explicit_residual_type<
  		 residual_t>::value,
  "OOPS: RESIDUAL_TYPE IN SELECTED EXPLICIT STEPPER IS NOT VALID");
  static_assert( meta::is_legitimate_explicit_residual_policy<
  		 residual_policy_t>::value,
  "RESIDUAL_POLICY NOT ADMISSIBLE, MAYBE \
NOT A CHILD OF EXPLICIT POLICY BASE");

public:
  order_t order() const{
    return order_value;
  }

  template<typename step_t>
  void operator()(state_t & yinout, scalar_t t,
		  scalar_t dt, step_t step){
    this->underlying()(yinout, t, dt, step);
  }

private:    
  ExplicitStepperBase() = default;
  ~ExplicitStepperBase() = default;

private:
  friend stepper_type;
  friend core::details::CrtpBase<ExplicitStepperBase<stepper_type>>;
  
};//end class

}}//end namespace rompp::ode
#endif
