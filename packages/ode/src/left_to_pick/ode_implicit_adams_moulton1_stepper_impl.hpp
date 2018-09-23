
#ifndef ODE_IMPLICIT_ADAMS_MOULTON1_STEPPER_IMPL_HPP_
#define ODE_IMPLICIT_ADAMS_MOULTON1_STEPPER_IMPL_HPP_

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
	 typename residual_policy_type,
	 typename jacobian_policy_type>
class implicitAdamsMoulton1StepperImpl<state_type,
				       residual_type,
				       jacobian_type,
				       scalar_type,
				       model_type,
				       time_type,
				       sizer_type,
				       solver_policy_type,
				       residual_policy_type,
				       jacobian_policy_type>
  : public ImplicitStepperBase<
  implicitAdamsMoulton1StepperImpl<state_type, residual_type,
				   jacobian_type, scalar_type,
				   model_type, time_type,
				   sizer_type,
				   solver_policy_type,
				   residual_policy_type,
				   jacobian_policy_type> >
{ 
  static_assert( meta::isLegitimateImplicitAdamsMoulton1ResidualPolicy<
		 residual_policy_type>::value,
  		 "IMPLICIT ADAMS-MOULTON1 Residual_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");
  static_assert( meta::isLegitimateImplicitAdamsMoulton1JacobianPolicy<
		 jacobian_policy_type>::value,
  		 "IMPLICIT ADAMS-MOULTON1 JACOBIAN_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");
 
private:
  using stepper_t = implicitAdamsMoulton1StepperImpl<state_type,
						     residual_type,
						     jacobian_type,
						     scalar_type,
						     model_type,
						     time_type,
						     sizer_type,
						     solver_policy_type,
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
	     typename U = residual_policy_type,
	     typename T = jacobian_policy_type,
	     typename... Args>
  implicitAdamsMoulton1StepperImpl(M & model,
				   S & solver,
				   U & res_policy_obj,
				   T & jac_policy_obj,
				   Args&& ... rest)
    : stepper_base_t(model, solver, res_policy_obj, jac_policy_obj),
      y_nm1_(std::forward<Args>(rest)...),
      Rnm1_(std::forward<Args>(rest)...)
  {}
    
  implicitAdamsMoulton1StepperImpl() = delete;
  ~implicitAdamsMoulton1StepperImpl() = default;
  
protected:

  template<typename step_t>
  void doStepImpl(state_type & y, time_type t, time_type dt, step_t step )
  {
    dt_ = dt;
    t_ = t;
    y_nm1_ = y;
    residual_obj_->compute(y, Rnm1_, *model_, t_);
    solver_->solve(y, *this);
  }//end doStepImpl

  
  //------------------------------------------------
  // residual policy is STANDARD  
  //------------------------------------------------
  template <typename U = state_type,
	    typename T = residual_type,
	    typename V = residual_policy_type,
  	    typename
  	    std::enable_if<
  	      ode::meta::isImplicitAdamsMoulton1ResidualStandardPolicy<V>::value
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
  	      !ode::meta::isImplicitAdamsMoulton1ResidualStandardPolicy<V>::value
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
  	      ode::meta::isImplicitAdamsMoulton1JacobianStandardPolicy<V>::value
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
  	      !ode::meta::isImplicitAdamsMoulton1JacobianStandardPolicy<V>::value
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
  void residualImplImpl(const state_type & y, residual_type & R)
  {
    // first compute space residual
    residual_obj_->compute(y, R, *model_, t_);

    // On input: R contains the application RHS evaluated for y
    //           dudt = f(x,u,...), R contains f(...) 
    // or any modifications to it provided by the policy
    for (decltype(R.size()) i=0; i < R.size(); i++){
      R[i] = y[i] - y_nm1_[i] - c1*dt_*R[i] - c1*dt_*Rnm1_[i];
    }
  }

  void jacobianImplImpl(const state_type & y, jacobian_type & J)
  {
    // compute the jacobian of the target model, however this happens
    jacobian_obj_->compute(y, J, *model_, t_);

    // this is bad, but for now leave it
    jacobian_type A( J.rows(),J.cols() );
    A.setIdentity();
    J.scale(-c1*dt_);
    J += A;
  }
  
private:
  friend stepper_base_t;

  state_type y_nm1_;
  residual_type Rnm1_;

  const time_type c1 = static_cast<time_type>(1)/2;
  
}; //end class

}//end namespace impl
}//end namespace ode  
}//end namespace rompp
#endif 
