
#ifndef ODE_IMPLICIT_EULER_STEPPER_IMPL_HPP_
#define ODE_IMPLICIT_EULER_STEPPER_IMPL_HPP_

#include "../base/ode_implicit_stepper_base.hpp"
#include "policies/implicit_policies/ode_implicit_policies_meta.hpp"

namespace ode{
namespace impl{

template<typename state_type,
	 typename residual_type,
	 typename jacobian_type,
	 typename scalar_type,
	 typename model_type,
	 typename time_type,
   typename sizer_type,
	 typename solver_policy_type,
	 typename residual_policy_type,
	 typename jacobian_policy_type>
class implicitEulerStepperImpl<state_type,
			       residual_type,
			       jacobian_type,
			       scalar_type,
			       model_type,
			       time_type, 
             sizer_type,
			       solver_policy_type,
			       residual_policy_type,
			       jacobian_policy_type>
  : public implicitStepperBase<
  implicitEulerStepperImpl<state_type, residual_type,
			   jacobian_type, scalar_type,
			   model_type, time_type, 
         sizer_type,
			   solver_policy_type,
			   residual_policy_type,
			   jacobian_policy_type> >
{
  static_assert( meta::isLegitimateImplicitEulerResidualPolicy<residual_policy_type>::value,
		 "IMPLICIT EULER RESIDUAL_POLICY NOT ADMISSIBLE, MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");
  static_assert( meta::isLegitimateImplicitEulerJacobianPolicy<jacobian_policy_type>::value,
		 "IMPLICIT EULER JACOBIAN_POLICY NOT ADMISSIBLE, MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");
  
private:
  using stepper_t = implicitEulerStepperImpl<state_type,
					     residual_type,
					     jacobian_type,
					     scalar_type,
					     model_type,
					     time_type, 
               sizer_type,
					     solver_policy_type,
					     residual_policy_type,
					     jacobian_policy_type>;
  using stepper_base_t = implicitStepperBase<stepper_t>;

protected:
  using stepper_base_t::model_;
  using stepper_base_t::solver_;
  using stepper_base_t::residual_obj_;
  using stepper_base_t::jacobian_obj_;
  using stepper_base_t::t_;
  using stepper_base_t::dt_;

protected:
  template < typename M = model_type,
	     typename S = solver_policy_type,
	     typename U = residual_policy_type,
	     typename T = jacobian_policy_type,
	     typename... Args>
  implicitEulerStepperImpl(M & model,
			   S & solver,
			   U & res_policy_obj,
			   T & jac_policy_obj,
			   Args&& ... rest)
    : stepper_base_t(model, solver, res_policy_obj, jac_policy_obj), 
      y_nm1_(std::forward<Args>(rest)...)
  {}
    
  implicitEulerStepperImpl() = delete;
  ~implicitEulerStepperImpl() = default;
  
protected:
  template<typename step_t>
  void doStepImpl(state_type & y, time_type t, time_type dt, step_t step )
  {
    dt_ = dt;
    t_ = t;
    y_nm1_ = y;
    solver_->solve(y, *this);
  }//end doStepImpl


  //------------------------------------------------
  // residual policy is STANDARD  
  //------------------------------------------------
  template <typename U = state_type,
	    typename T = residual_type,
	    typename
	    std::enable_if<
	      ode::meta::isImplicitEulerResidualStandardPolicy<T>::value
	      >::type * = nullptr
	    >
  void residualImpl(const U & y, T & R){
    // first compute space residual
    residual_obj_->compute(y, R, *model_, t_);

    // Here, R contains the application RHS or whatever, i.e. if
    //           dudt = f(x,u,...), R contains f(...)
    // or any modifications to it provided by the policy
    // now compute time discrete residual
    for (decltype(R.size()) i=0; i < R.size(); i++){
      R[i] = y[i] - y_nm1_[i] - dt_*R[i];
    }    
  }

  //------------------------------------------------
  // residual policy is ARBITRARY:
  // call additional to weight the time discrete residual
  // after it is computed
  //------------------------------------------------
  template <typename U = state_type,
  	    typename T = residual_type,
  	    typename
  	    std::enable_if<
  	      !ode::meta::isImplicitEulerResidualStandardPolicy<T>::value
  	      >::type * = nullptr
  	    >
  void residualImpl(const U & y, T & R){
    // first compute space residual
    residual_obj_->compute(y, R, *model_, t_);

    // Here, R contains the application RHS or whatever, i.e. if
    //           dudt = f(x,u,...), R contains f(...)
    // or any modifications to it provided by the policy
    // now compute time discrete residual
    for (decltype(R.size()) i=0; i < R.size(); i++){
      R[i] = y[i] - y_nm1_[i] - dt_*R[i];
    }
    residual_obj_->weightTimeDiscreteResidual(y, R, *model_, t_);    
  }

  
  //------------------------------------------------
  // jacobian policy is STANDARD
  //------------------------------------------------
  template <typename U = state_type,
	    typename T = jacobian_type,
	    typename
	    std::enable_if<
	      ode::meta::isImplicitEulerJacobianStandardPolicy<T>::value
	      >::type * = nullptr
	    >
  void jacobianImpl(const U & y, T & J)
  {
    // compute the jacobian of the target model, however this happens
    jacobian_obj_->compute(y, J, *model_, t_);

    // this is bad, but for now leave it
    jacobian_type A( J.rows(),J.cols() );
    A.setIdentity();

    J.scale(-dt_);
    J += A;
  }
  
  
  //------------------------------------------------
  // jacobian policy is ARBITRARY
  // call additional to weight the time discrete residual
  // after it is computed
  //------------------------------------------------
  template <typename U = state_type,
  	    typename T = residual_type,
  	    typename
  	    std::enable_if<
  	      !ode::meta::isImplicitEulerJacobianStandardPolicy<T>::value
  	      >::type * = nullptr
  	    >
  void jacobianImpl(const U & y, T & J)
  {
    // compute the jacobian of the target model, however this happens
    jacobian_obj_->compute(y, J, *model_, t_);

    // this is bad, but for now leave it
    jacobian_type A( J.rows(),J.cols() );
    A.setIdentity();
    J.scale(-dt_);
    J += A;

    // leftmultiply to scale
    jacobian_obj_->weightJacobian(y, J, *model_, t_);
  }
  

private:
  friend stepper_base_t;

  state_type y_nm1_;

}; //end class

}//end namespace impl
}//end namespace ode  
#endif 
