
#ifndef ODE_IMPLICIT_STEPPER_BASE_HPP_
#define ODE_IMPLICIT_STEPPER_BASE_HPP_

#include "ode_ConfigDefs.hpp"
#include "../ode_stepper_traits.hpp"
#include "../../ode_meta.hpp"
#include "../../policies/ode_implicit_policies_meta.hpp"


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
  using resizer_t = typename traits::resizer_t;
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
public:
  order_t order() const{
    return order_value;
  }
  void doStep(state_t &inout, time_t t, time_t dt){
    this->stepper()->doStepImpl( inout, t, dt );
  }
  void residual(const state_t & y, state_t & R){
    this->stepper()->residualImpl(y, R);
  }
  void jacobian(const state_t & y, jacobian_t & J){
    this->stepper()->jacobianImpl(y, J);
  }

private:
  //*********************************************************
  // residual policy = NOT standard (passed from derived)
  // jacobian policy = NOT standard (passed from derived)
  //*********************************************************
  implicitStepperBase(model_t & model,
		      residual_policy_t & res_policy_obj,
		      jacobian_policy_t & jac_policy_obj)
    : model_(&model),
      residual_policy_obj_(res_policy_obj),
      jacobian_policy_obj_(jac_policy_obj),
      standardResidualPolicyActive_(false), standardJacobianPolicyActive_(false)
  {}

  
  //*********************************************************
  // residual policy = standard (not passed from derived)
  // jacobian policy = NOT standard (passed from derived)
  //*********************************************************
  template <typename U = residual_policy_t,
  	    typename std::enable_if<
  	      core::meta::is_default_constructible<U>::value 
  	      >::type * = nullptr
  	    >
  implicitStepperBase(model_t & model,
		      /* residual policy obj not given because is standard*/
		      jacobian_policy_t & jac_policy_obj)
    : model_(&model),
      residual_policy_obj_(new U()),
      jacobian_policy_obj_(jac_policy_obj),
      standardResidualPolicyActive_(true), standardJacobianPolicyActive_(false)
  {}

  //*********************************************************
  // residual policy = NOT standard (passed from derived)
  // jacobian policy = standard (not passed from derived)
  //*********************************************************
  template <typename U = jacobian_policy_t,
  	    typename std::enable_if<
  	      core::meta::is_default_constructible<U>::value 
  	      >::type * = nullptr
  	    >
  implicitStepperBase(model_t & model,
		      residual_policy_t & res_policy_obj
		      /* jacobian policy obj not given because is standard*/)
    : model_(&model),
      residual_policy_obj_(res_policy_obj),
      jacobian_policy_obj_(new U()),
      standardResidualPolicyActive_(false), standardJacobianPolicyActive_(true)
  {}

  //*********************************************************
  // residual policy = standard (not passed from derived)
  // jacobian policy = standard (not passed from derived)
  //*********************************************************
  template <typename M = model_t,
	    typename U = residual_policy_t,
	    typename T = jacobian_policy_t,
  	    typename std::enable_if<
  	      core::meta::is_default_constructible<U>::value &&
  	      core::meta::is_default_constructible<T>::value 
  	      >::type * = nullptr
  	    >
  implicitStepperBase(M & model
		      /*residual policy obj not given because is standard
		        jacobian policy obj not given because is standard*/)
    : model_(&model),
      residual_policy_obj_(new U()),
      jacobian_policy_obj_(new T()),
      standardResidualPolicyActive_(true), standardJacobianPolicyActive_(true)
  {}

  
  ~implicitStepperBase(){
    if (standardResidualPolicyActive_)
      delete residual_policy_obj_;
    if (standardJacobianPolicyActive_)
      delete jacobian_policy_obj_;
  }
  
private:
  friend stepper_type;
  stepper_type * stepper( ){
    return static_cast< stepper_type* >( this );
  }
  const stepper_type * stepper( void ) const{
    return static_cast< const stepper_type* >( this );
  }

protected:
  resizer_t myResizer_;
  model_t * model_;
  residual_policy_t * residual_policy_obj_;
  jacobian_policy_t * jacobian_policy_obj_;
private:
  bool standardResidualPolicyActive_;
  bool standardJacobianPolicyActive_;

};


}//end namespace  
#endif
