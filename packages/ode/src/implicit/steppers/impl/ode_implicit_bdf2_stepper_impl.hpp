
#ifndef ODE_STEPPERS_IMPLICIT_STEPPERS_IMPL_IMPLICIT_BDF2_STEPPER_IMPL_HPP_
#define ODE_STEPPERS_IMPLICIT_STEPPERS_IMPL_IMPLICIT_BDF2_STEPPER_IMPL_HPP_

#include "../base/ode_implicit_stepper_base.hpp"
#include "../../policies/meta/ode_implicit_bdf2_policies_meta.hpp"

namespace rompp{ namespace ode{ namespace impl{

template<typename state_type,
	 typename residual_type,
	 typename jacobian_type,
	 typename model_type,
	 typename aux_stepper_type, 
	 typename residual_policy_type,
	 typename jacobian_policy_type>
class ImplicitBDF2StepperImpl<state_type,
			      residual_type,
			      jacobian_type,
			      model_type,
			      aux_stepper_type,
			      residual_policy_type,
			      jacobian_policy_type>
  : public ImplicitStepperBase<
             ImplicitBDF2StepperImpl<state_type,
				     residual_type,
				     jacobian_type,
				     model_type,
				     aux_stepper_type,
				     residual_policy_type,
				     jacobian_policy_type> >,
    private OdeStorage<state_type, residual_type, 2>,
    private ImpOdeAuxData<model_type,
			  typename core::details::traits<state_type>::scalar_t,
			  residual_policy_type,
			  jacobian_policy_type>{
  
  static_assert( meta::is_legitimate_implicit_bdf2_residual_policy<
		 residual_policy_type>::value ||
		 meta::is_implicit_bdf2_residual_standard_policy<
		 residual_policy_type>::value,
		 "IMPLICIT BDF2 RESIDUAL_POLICY NOT ADMISSIBLE,\
MAYBE NOT A CHILD OR DERIVING FROM WRONG BASE");

  static_assert( meta::is_legitimate_implicit_bdf2_jacobian_policy<
		 jacobian_policy_type>::value,
		 "IMPLICIT BDF2 JACOBIAN_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OR DERIVING FROM WRONG BASE");


  using stepper_t = ImplicitEulerStepperImpl<state_type,
					     residual_type,
					     jacobian_type,
					     model_type,
					     residual_policy_type,
					     jacobian_policy_type>;
  using scalar_type  = typename core::details::traits<state_type>::scalar_t;
  using stepper_base_t = ImplicitStepperBase<stepper_t>;
  using storage_base_t = OdeStorage<state_type, residual_type, 2>;
  using auxdata_base_t = ImpOdeAuxData<model_type, scalar_type,
				       residual_policy_type,
				       jacobian_policy_type>;

public:
  // these aliases are needed by the solver
  using vector_type = state_type;
  using matrix_type = jacobian_type;
  
protected:
  using storage_base_t::auxStates_;
  using auxdata_base_t::model_;
  using auxdata_base_t::residual_obj_;
  using auxdata_base_t::jacobian_obj_;
  using auxdata_base_t::t_;
  using auxdata_base_t::dt_;

protected:
  template < typename M = model_type,
	     typename U = residual_policy_type,
	     typename T = jacobian_policy_type,
	     typename V = aux_stepper_type,
	     typename T3 = state_type>
  ImplicitBDF2StepperImpl(const M & model,
			  const U & res_policy_obj,
			  const T & jac_policy_obj,
			  V & auxStepper,
			  const T3 & y0)
    : storage_base_t(y0),
      auxdata_base_t(model, res_policy_obj, jac_policy_obj),
      auxStp_(&auxStepper){}
    
  ImplicitBDF2StepperImpl() = delete;
  virtual ~ImplicitBDF2StepperImpl() = default;
  
public:
  template<typename solver_type, typename step_t>
  void operator()(state_type & y, scalar_type t,
		  scalar_type dt, step_t step,
		  solver_type & solver){
    dt_ = dt;
    t_ = t;

    // first step, use auxiliary stepper
    if (step == 1){
      auxStates_[0] = y;
      (*auxStp_)(y, t, dt, step, solver);
    }
    if (step == 2){
      auxStates_[1] = y;
      y = solver.solve(*this, y);
    }
    // if (step >= 3){
    //   auxStates_[0] = auxStates_[1];
    //   auxStates_[1] = y;
    //   y = solver.solve(*this, y);
    // }
  }//end doStepImpl

public:
  
  void residualImpl(const state_type & y, residual_type & R)const{
    (*residual_obj_)(y, R, auxStates_, *model_, t_, dt_);
  }
  //--------------------------------------------------------

  void jacobianImpl(const state_type & y, jacobian_type & J)const{
    (*jacobian_obj_)(y, J, *model_, t_, dt_);
  }
  //--------------------------------------------------------

  residual_type residualImpl(const state_type & y)const{
    return (*residual_obj_)(y, auxStates_, *model_, t_, dt_);
  }
  //--------------------------------------------------------

  jacobian_type jacobianImpl(const state_type & y)const{
    return (*jacobian_obj_)(y, *model_, t_, dt_);
  }
  //--------------------------------------------------------

  
private:
  friend stepper_base_t;  
  aux_stepper_type * auxStp_;
  
}; //end class

}//end namespace impl
}//end namespace ode  
}//end namespace rompp
#endif 





// private:
//   void residualImplImpl(const state_type & y, residual_type & R){
//     // first compute space residual
//     residual_obj_->compute(y, R, *model_, t_);

//     // Here, R contains the application RHS or whatever, i.e. if
//     //           dudt = f(x,u,...), R contains f(...)
//     // or any modifications to it provided by the policy
//     // now compute time discrete residual
//     for (decltype(R.size()) i=0; i < R.size(); i++){
//       R[i] = y[i] - c1_*y_nm1_[i] + c2_*y_nm2_[i] - c3_*dt_*R[i];
//     }
//   }

//   void jacobianImplImpl(const state_type & y, jacobian_type & J)
//   {
//     // compute the jacobian of the target model, however this happens
//     jacobian_obj_->compute(y, J, *model_, t_);

//     // this is bad, but for now leave it
//     jacobian_type A( J.rows(),J.cols() );
//     A.setIdentity();
//     J.scale(-c3_*dt_);
//     J += A;
//   }
