
#ifndef ODE_IMPLICIT_BDF3_STEPPER_IMPL_HPP_
#define ODE_IMPLICIT_BDF3_STEPPER_IMPL_HPP_

#include "../base/ode_ImplicitStepperBase.hpp"

namespace rompp{
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
	 typename aux_stepper_type, 
	 typename residual_policy_type,
	 typename jacobian_policy_type>
class implicitBDF3StepperImpl<state_type,
			      residual_type,
			      jacobian_type,
			      scalar_type,
			      model_type,
			      time_type,
			      sizer_type,
			      solver_policy_type,
			      aux_stepper_type,
			      residual_policy_type,
			      jacobian_policy_type>
  : public ImplicitStepperBase<
  implicitBDF3StepperImpl<state_type, residual_type,
			  jacobian_type, scalar_type,
			  model_type, time_type,
			  sizer_type,
			  solver_policy_type,
			  aux_stepper_type,
			  residual_policy_type,
			  jacobian_policy_type> >
{ 

  static_assert( meta::isLegitimateImplicitBDF3ResidualPolicy<
		 residual_policy_type>::value,
  		 "IMPLICIT BDF3 RESIDUAL_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");
  static_assert( meta::isLegitimateImplicitBDF3JacobianPolicy<
		 jacobian_policy_type>::value,
  		 "IMPLICIT BDF3 JACOBIAN_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");
 
  using stepper_t = implicitBDF3StepperImpl<state_type,
					    residual_type,
					    jacobian_type,
					    scalar_type,
					    model_type,
					    time_type,
					    sizer_type,
					    solver_policy_type,
					    aux_stepper_type,
					    residual_policy_type,
					    jacobian_policy_type>;
  using stepper_base_t = ImplicitStepperBase<stepper_t>;
  
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
	     typename V = aux_stepper_type,
	     typename U = residual_policy_type,
	     typename T = jacobian_policy_type,
	     typename... Args>
  implicitBDF3StepperImpl(M & model,
			  S & solver,
			  V & auxStepper,
			  U & res_policy_obj,
			  T & jac_policy_obj,
			  Args&& ... rest)
    : stepper_base_t(model, solver, res_policy_obj, jac_policy_obj),
      y_nm1_(std::forward<Args>(rest)...),
      y_nm2_(std::forward<Args>(rest)...),
      y_nm3_(std::forward<Args>(rest)...),
      auxStp_(&auxStepper)
  {}
    
  implicitBDF3StepperImpl() = delete;
  ~implicitBDF3StepperImpl() = default;
  
protected:
  template<typename step_t>
  void doStepImpl(state_type & y, time_type t, time_type dt, step_t step )
  {
    // first few steps, use auxiliary stepper
    if (step == 1){
      y_nm3_ = y;
      auxStp_->doStep(y, t, dt, step);
    }
    if (step == 2){
      y_nm2_ = y;
      auxStp_->doStep(y, t, dt, step);
    }
    if (step == 3){
      y_nm1_ = y;
      auxStp_->doStep(y, t, dt, step);
    }
    if (step >= 4){
      y_nm3_ = y_nm2_;
      y_nm2_ = y_nm1_;
      y_nm1_ = y;
      dt_ = dt;
      t_ = t;
      solver_->solve(y, *this);
    }

  }//end doStepImpl
  
  //------------------------------------------------
  // residual policy is STANDARD  
  //------------------------------------------------
  template <typename U = state_type,
      typename T = residual_type,
      typename V = residual_policy_type,
  	    typename
  	    std::enable_if<
  	      ode::meta::isImplicitBDF3ResidualStandardPolicy<V>::value
  	      >::type * = nullptr
  	    >
  void residualImpl(const U & y, T & R){
    residualImplImpl(y, R);
  }

  //------------------------------------------------
  // residual policy is ARBITRARY:
  // weight the time discrete residual after it is computed
  //------------------------------------------------
  template <typename U = state_type,
      typename T = residual_type,
      typename V = residual_policy_type,
  	    typename
  	    std::enable_if<
  	      !ode::meta::isImplicitBDF3ResidualStandardPolicy<V>::value
  	      >::type * = nullptr
  	    >
  void residualImpl(const U & y, T & R){
    residualImplImpl(y, R);
    residual_obj_->weightTimeDiscreteResidual(y, R, *model_, t_);
  }

  //------------------------------------------------
  // jacobian policy is STANDARD
  //------------------------------------------------
  template <typename U = state_type,
      typename T = jacobian_type,
      typename V = jacobian_policy_type,
  	    typename
  	    std::enable_if<
  	      ode::meta::isImplicitBDF3JacobianStandardPolicy<V>::value
  	      >::type * = nullptr
  	    >
  void jacobianImpl(const U & y, T & J)
  {
    jacobianImplImpl(y,J);
  }
  
  //------------------------------------------------
  // jacobian policy is ARBITRARY
  // call additional to weight the time discrete residual
  // after it is computed
  //------------------------------------------------
  template <typename U = state_type,
      typename T = jacobian_type,
      typename V = jacobian_policy_type,
  	    typename
  	    std::enable_if<
  	      !ode::meta::isImplicitBDF3JacobianStandardPolicy<V>::value
  	      >::type * = nullptr
  	    >
  void jacobianImpl(const U & y, T & J)
  {
    jacobianImplImpl(y,J);
    // leftmultiply to scale
    jacobian_obj_->weightJacobian(y, J, *model_, t_);
  }
  //------------------------------------------------
  

private:
  void residualImplImpl(const state_type & y, residual_type & R){
    // first compute space residual
    residual_obj_->compute(y, R, *model_, t_);

    // Here, R contains the application RHS or whatever, i.e. if
    //           dudt = f(x,u,...), R contains f(...)
    // or any modifications to it provided by the policy
    // now compute time discrete residual
    for (decltype(R.size()) i=0; i < R.size(); i++){
      R[i] = y[i] - c1_*y_nm1_[i] + c2_*y_nm2_[i] - c3_*y_nm3_[i] - c4_*dt_*R[i];
    }
  }

  void jacobianImplImpl(const state_type & y, jacobian_type & J)
  {
    // compute the jacobian of the target model, however this happens
    jacobian_obj_->compute(y, J, *model_, t_);

    // this is bad, but for now leave it
    jacobian_type A( J.rows(),J.cols() );
    A.setIdentity();
    J.scale(-c3_*dt_);
    J += A;
  }

private:
  friend stepper_base_t;
  
  state_type y_nm1_;
  state_type y_nm2_;
  state_type y_nm3_;
  aux_stepper_type * auxStp_;

  const scalar_type c1_ = static_cast<scalar_type>(18)/11;
  const scalar_type c2_ = static_cast<scalar_type>(9)/11;
  const scalar_type c3_ = static_cast<scalar_type>(2)/11;
  const scalar_type c4_ = static_cast<scalar_type>(6)/11;
  
}; //end class

}//end namespace impl
}//end namespace ode  
}//end namespace rompp
#endif 
