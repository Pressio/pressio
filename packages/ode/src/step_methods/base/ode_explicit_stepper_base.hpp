
#ifndef ODE_EXPLICIT_STEPPER_BASE_HPP_
#define ODE_EXPLICIT_STEPPER_BASE_HPP_

#include "ode_ConfigDefs.hpp"
#include "../ode_stepper_traits.hpp"
#include "../../ode_meta.hpp"
#include "../../policies/ode_explicit_policies_meta.hpp"


namespace ode{

template<typename stepper_type>
class explicitStepperBase
{
private:
  using step_traits = ode::details::traits<stepper_type>;
  using state_t = typename step_traits::state_t;
  using sc_t = typename step_traits::scalar_t;
  using order_t = typename step_traits::order_t;
  using time_t = typename step_traits::time_t;
  using residual_t = typename step_traits::residual_t;
  using residual_policy_t = typename step_traits::residual_policy_t;  
  using model_t = typename step_traits::model_t;
  static constexpr order_t order_value = step_traits::order_value;

  //do checking here that things are as supposed
  static_assert( meta::isLegitimateExplicitStateType<state_t>::value,
		 "OOPS: STATE_TYPE IN SELECTED EXPLICIT STEPPER IS NOT VALID");
  static_assert( meta::isLegitimateExplicitResidualType<residual_t>::value,
		 "OOPS: RESIDUAL_TYPE IN SELECTED EXPLICIT STEPPER IS NOT VALID");
  static_assert( meta::isLegitimateTimeType<time_t>::value,
		 "OOPS: TIME_TYPE IN SELECTED EXPLICIT STEPPER IS NOT VALID");
  static_assert( meta::derivesFromExplicitResidualPolicyBase<residual_policy_t>::value,
		 "RESIDUAL_POLICY_TYPE DOES NOT INHERIT FROM EXPLICIT POLICY BASE");
  
public:
  order_t order() const{
    return order_value;
  }

  void doStep(state_t &inout, time_t t, time_t dt ){
    this->stepper()->doStepImpl(inout, t, dt);
  }
  
private:
  friend stepper_type;
    
  // constructor when policy is NOT standard.
  // we then we store ptr to the incoming object
  explicitStepperBase(model_t & model,
		      residual_policy_t & res_policy_obj)
    : model_(&model),
      residual_policy_obj_(&res_policy_obj),
      standardPolicyActive_(false)
  {}

  // constructor when the policy is standard
  // when policy is standard, because standard policies
  // need to be default constructible, we simply create
  // a new object and store it, and delete it in destructor.
  template <typename U = residual_policy_t,
  	    typename std::enable_if<
  	      core::meta::is_default_constructible<U>::value 
  	      >::type * = nullptr
  	    >
  explicitStepperBase(model_t & model)
    : model_(&model),
      residual_policy_obj_(new residual_policy_t()),
      standardPolicyActive_(true)
  {}

  explicitStepperBase() = delete;

  ~explicitStepperBase(){
    if (standardPolicyActive_)
      delete residual_policy_obj_;
  }

  stepper_type * stepper( ){
    return static_cast< stepper_type* >( this );
  }
  const stepper_type * stepper( void ) const{
    return static_cast< const stepper_type* >( this );
  }

protected:
  model_t * model_;
  residual_policy_t * residual_policy_obj_;
private:
  bool standardPolicyActive_;
};


}//end namespace
  
#endif
