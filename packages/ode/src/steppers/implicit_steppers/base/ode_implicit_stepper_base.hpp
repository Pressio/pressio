
#ifndef ODE_IMPLICIT_STEPPER_BASE_HPP_
#define ODE_IMPLICIT_STEPPER_BASE_HPP_

#include "ode_ConfigDefs.hpp"
#include "../ode_implicit_stepper_traits.hpp"
#include "../../../meta/ode_meta.hpp"
#include "../../../meta/ode_meta_implicit.hpp"
#include "../../../policies/implicit_policies/ode_implicit_policies_meta.hpp"

namespace ode{

template<typename stepper_type>
class implicitStepperBase
{
private:
  using traits = typename ode::details::traits<stepper_type>;

  using state_t = typename traits::state_t;
  using residual_t = typename traits::residual_t;
  using jacobian_t = typename traits::jacobian_t;
  using sc_t = typename traits::scalar_t;
  using model_t = typename traits::model_t;
  using time_t = typename traits::time_t;
  using residual_policy_t = typename traits::residual_policy_t;  
  using jacobian_policy_t = typename traits::jacobian_policy_t;  
  using order_t = typename traits::order_t; 
  static constexpr order_t order_value = ode::details::traits<stepper_type>::order_value;

  //do checking here that things are as supposed
  static_assert( meta::isLegitimateImplicitStateType<state_t>::value,
		 "OOPS: STATE_TYPE IN SELECTED IMPLICIT STEPPER IS NOT VALID");
  static_assert( meta::isLegitimateImplicitResidualType<residual_t>::value,
		 "OOPS: RESIDUAL_TYPE IN SELECTED IMPLICIT STEPPER IS NOT VALID");
  static_assert( meta::isLegitimateJacobianType<jacobian_t>::value,
		 "OOPS: JACOBIAN_TYPE IN SELECTED IMPLICIT STEPPER IS NOT VALID");
  static_assert( meta::isLegitimateTimeType<time_t>::value,
		 "OOPS: TIME_TYPE IN SELECTED IMPLICIT STEPPER IS NOT VALID");
  static_assert( meta::isLegitimateImplicitEulerResidualPolicy<residual_policy_t>::value,
		 "RESIDUAL_POLICY NOT ADMISSIBLE, MAYBE NOT A CHILD OF IMPLICIT POLICY BASE");
  static_assert( meta::isLegitimateImplicitEulerJacobianPolicy<jacobian_policy_t>::value,
		 "RESIDUAL_POLICY NOT ADMISSIBLE, MAYBE NOT A CHILD OF IMPLICIT POLICY BASE");

public:
  order_t order() const{
    return order_value;
  }
  void doStep(state_t & y, time_t t, time_t dt){
    this->stepper()->doStepImpl( y, t, dt );
  }
  void residual(const state_t & y, state_t & R){
    this->stepper()->residualImpl(y, R);
  }
  void jacobian(const state_t & y, jacobian_t & J){
    this->stepper()->jacobianImpl(y, J);
  }

private:
  implicitStepperBase(model_t & model,
		      residual_policy_t & res_policy_obj,
		      jacobian_policy_t & jac_policy_obj)
    : model_(&model),
      residual_policy_obj_(&res_policy_obj),
      jacobian_policy_obj_(&jac_policy_obj)
  {}
  
  ~implicitStepperBase(){}
  
private:
  friend stepper_type;

  stepper_type * stepper( ){
    return static_cast< stepper_type* >( this );
  }
  const stepper_type * stepper( void ) const{
    return static_cast< const stepper_type* >( this );
  }

protected:
  model_t * model_;
  residual_policy_t * residual_policy_obj_;
  jacobian_policy_t * jacobian_policy_obj_;

};//end class
}//end namespace  
#endif
