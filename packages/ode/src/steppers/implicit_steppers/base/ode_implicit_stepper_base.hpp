
#ifndef ODE_STEPPERS_IMPLICIT_STEPPERS_BASE_IMPLICIT_STEPPER_BASE_HPP_
#define ODE_STEPPERS_IMPLICIT_STEPPERS_BASE_IMPLICIT_STEPPER_BASE_HPP_

#include "ode_ConfigDefs.hpp"
#include "../ode_implicit_stepper_traits.hpp"
#include "../../../meta/ode_meta.hpp"
#include "../../../meta/ode_meta_implicit.hpp"
#include "../../../policies/meta/ode_implicit_policies_meta.hpp"
#include "../../../ode_storage.hpp"
#include "../../../ode_aux_data.hpp"

namespace ode{

template<typename stepper_type>
class ImplicitStepperBase
  : private core::details::CrtpBase<ImplicitStepperBase<stepper_type>>
{
private:
  using traits = typename ode::details::traits<stepper_type>;

  using state_t = typename traits::state_t;
  using residual_t = typename traits::residual_t;
  using jacobian_t = typename traits::jacobian_t;
  using sc_t = typename traits::scalar_t;
  using residual_policy_t = typename traits::residual_policy_t;  
  using jacobian_policy_t = typename traits::jacobian_policy_t;  

  using order_t = typename traits::order_t; 
  static constexpr order_t order_value =
    ode::details::traits<stepper_type>::order_value;

  //do checking here that things are as supposed
  static_assert( meta::isLegitimateImplicitStateType<state_t>::value,
       "OOPS: STATE_TYPE IN SELECTED IMPLICIT STEPPER IS NOT VALID");
  static_assert( meta::isLegitimateImplicitResidualType<residual_t>::value,
       "OOPS: RESIDUAL_TYPE IN SELECTED IMPLICIT STEPPER IS NOT VALID");
  static_assert( meta::isLegitimateJacobianType<jacobian_t>::value,
       "OOPS: JACOBIAN_TYPE IN SELECTED IMPLICIT STEPPER IS NOT VALID");

public:
  order_t order() const{
    return order_value;
  }

  template <typename solver_type,
	    typename step_t>
  void doStep(state_t & y, sc_t t,
	      sc_t dt, step_t step,
	      solver_type & solver){
    this->underlying().doStepImpl( y, t, dt, step, solver);
  }

  void residual(const state_t & y, state_t & R){
    this->underlying().residualImpl(y, R);
  }

  void jacobian(const state_t & y, jacobian_t & J){
    this->underlying().jacobianImpl(y, J);
  }

private:
  
  ImplicitStepperBase() = default;
  ~ImplicitStepperBase() = default;
  
private:
  friend stepper_type;
  friend core::details::CrtpBase<ImplicitStepperBase<stepper_type>>;

};//end class

}//end namespace  
#endif
