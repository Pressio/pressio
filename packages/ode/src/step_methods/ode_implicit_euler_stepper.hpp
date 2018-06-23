
#ifndef ODE_IMPLICIT_EULER_STEPPER_HPP_
#define ODE_IMPLICIT_EULER_STEPPER_HPP_

#include "./base/ode_implicit_stepper_base.hpp"

namespace ode{

template<typename state_type,
	 typename residual_type,
	 typename jacobian_type,
	 typename scalar_type,
	 typename state_resizer_fnctor_type,
	 typename model_type,
	 typename time_type,
	 typename solver_policy_type,
	 typename residual_policy_type,
	 typename jacobian_policy_type>
class implicitEulerStepper<state_type,residual_type,jacobian_type,scalar_type,
			   state_resizer_fnctor_type,model_type,time_type,
			   solver_policy_type,residual_policy_type,jacobian_policy_type,
			   typename 
			   std::enable_if<
			     !std::is_same<state_type,void>::value &&
			     core::meta::is_default_constructible<state_resizer_fnctor_type
								  >::value
			     >::type
			   >
  : public implicitStepperBase<
  implicitEulerStepper<state_type,residual_type,jacobian_type,scalar_type,
		       state_resizer_fnctor_type,model_type,time_type,
		       solver_policy_type,residual_policy_type,jacobian_policy_type>>
{

  static_assert(meta::derivesFromImplicitEulerResidualPolicyBase<residual_policy_type>::value,
		"RESIDUAL_POLICY_TYPE DOES NOT INHERIT FROM IMPLICIT EULER RESIDUAL POLICY BASE");
  static_assert(meta::derivesFromImplicitEulerJacobianPolicyBase<jacobian_policy_type>::value,
		"JACOBIAN_POLICY_TYPE DOES NOT INHERIT FROM IMPLICIT EULER JACOBIAN POLICY BASE");

public :
  using stepper_t = implicitEulerStepper<state_type,residual_type,jacobian_type,scalar_type,
					 state_resizer_fnctor_type,model_type,time_type,
					 solver_policy_type,residual_policy_type,
					 jacobian_policy_type>;
  using stepper_base_t = implicitStepperBase<stepper_t>;

  //*********************************************************
  // residual policy = NOT standard 
  // jacobian policy = NOT standard 
  //*********************************************************
  template < typename M = model_type,
	     typename S = solver_policy_type,
	     typename U = residual_policy_type,
	     typename T = jacobian_policy_type,
	     typename
	     std::enable_if<
	       !meta::isImplicitEulerResidualStandardPolicy<U>::value &&
	       !meta::isImplicitEulerJacobianStandardPolicy<T>::value
	       >::type * = nullptr
	     >
  implicitEulerStepper(M & model,
		       S & solver,
		       U & res_policy_obj,
		       T & jac_policy_obj)
    : stepper_base_t(model, res_policy_obj, jac_policy_obj), solver_(&solver){}

  //*********************************************************
  // residual policy = standard 
  // jacobian policy = NOT standard 
  //*********************************************************
  template < typename M = model_type,
	     typename S = solver_policy_type,
	     typename U = residual_policy_type,
	     typename T = jacobian_policy_type,
	     typename std::enable_if<
	       meta::isImplicitEulerResidualStandardPolicy<U>::value &&
	       !meta::isImplicitEulerJacobianStandardPolicy<T>::value
	      >::type * = nullptr>
  implicitEulerStepper(M & model,
		       S & solver,
		       T & jac_policy_obj)
    : stepper_base_t(model, jac_policy_obj), solver_(&solver){}

  //*********************************************************
  // residual policy = NOT standard 
  // jacobian policy = standard 
  //*********************************************************
  template < typename M = model_type,
	     typename S = solver_policy_type,
	     typename U = residual_policy_type,
	     typename T = jacobian_policy_type,
	     typename
	     std::enable_if<
	       !meta::isImplicitEulerResidualStandardPolicy<U>::value &&
	       meta::isImplicitEulerJacobianStandardPolicy<T>::value
	       >::type * = nullptr
	     >
  implicitEulerStepper(M & model,
		       S & solver,
		       U & res_policy_obj)
    : stepper_base_t(model, res_policy_obj), solver_(&solver){}

  //*********************************************************
  // residual policy = standard 
  // jacobian policy = standard 
  //*********************************************************
  template < typename M = model_type,
	     typename S = solver_policy_type,
	     typename U = residual_policy_type,
	     typename T = jacobian_policy_type,
	     typename
	     std::enable_if<
	       meta::isImplicitEulerResidualStandardPolicy<U>::value &&
	       meta::isImplicitEulerJacobianStandardPolicy<T>::value
	       >::type * = nullptr
	     >
  implicitEulerStepper(M & model,
		       S & solver)
    : stepper_base_t(model), solver_(&solver){}

  //*********************************************************
    
  implicitEulerStepper() = delete;
  ~implicitEulerStepper() = default;
  
private:
  void doStepImpl(state_type & y_inout, time_type t, time_type dt ){
    y_nm1_ = y_inout;
    dt_ = dt;
    t_ = t;

    std::cout << " doStepImpl: right " << std::endl;
    for (int i=0; i<y_inout.size(); ++i)
      std::cout << std::setprecision(14) << y_inout[i]  << " ";
    std::cout << "------------------ " << std::endl;
    solver_->solve(y_inout, *this);

    for (int i=0; i<y_inout.size(); ++i)
      std::cout << std::setprecision(14) << y_inout[i]  << " ";
    std::cout << std::endl;    
    std::cout << std::endl;    
  }
  
  void residualImpl(const state_type & y, state_type & R){
    this->residual_policy_obj_->compute(y, y_nm1_, R, *(this->model_), t_, dt_);
  }
  void jacobianImpl(const state_type & y, jacobian_type & J){
    this->jacobian_policy_obj_->compute(y, J, *(this->model_), t_, dt_);
  }

private:
  friend stepper_base_t;
  time_type t_;
  time_type dt_;
  state_type y_nm1_;
  solver_policy_type * solver_;

}; //end class
}//end namespace
#endif 






  // void operator()(const state_type & y, state_type & R, jacobian_type & J)
  // {
  //   R.resize(y.size());
  //   assert(y.size() == yOld_.size());

  //   if (romOn_ == false)
  //     doMathFullCase(y, R, J);
  //   else
  //     doMathReducedCase(y, R, J);
  // }

  // void doMathReducedCase(const state_type & yRed, state_type & R, jacobian_type & J)
  // {    
  //   state_type Vdoty;
  //   (*sysFunctor_).rescaleState(yRed, Vdoty);
  //   state_type Vdotyold;
  //   (*sysFunctor_).rescaleState(yOld_, Vdotyold);
    
  //   // compute RHS(y)
  //   jacobian_type Jfull; //Jfull is resized inside functor
  //   (*sysFunctor_)(yRed, R, Jfull, t_);
    
  //   // R = y - yOld - dt*rhs(y)
  //   for (decltype(R.size()) i=0; i < R.size(); i++){
  //     R[i] = Vdoty[i] - Vdotyold[i] -dt_*R[i];
  //   }

  //   auto & jac = Jfull.getNonConstRefToData();
  //   jac(0,0) = 1.0 - dt_ * jac(0,0);
  //   for (size_t i=1; i < J.rows(); ++i){
  //     jac(i,i-1) = - dt_ * jac(i,i-1);
  //     jac(i,i) = 1.0 - dt_ * jac(i,i);
  //   }
  //   // J should have here reduced size, this is done inside function 
  //   (*sysFunctor_).rescaleJacobian(Jfull, J);
  // }


  // void doMathFullCase(const state_type & y, state_type & R, jacobian_type & J)
  // {
  //   R.resize(y.size());
  //   assert(y.size() == yOld_.size());
    
  //   // compute RHS(y), this is bad because we should not know that y,R and J are
  //   // my vectors types, these could be anything. Needs to be fixed.
  //   (*sysFunctor_)(y, R, J, t_);

  //   // R = y - yOld - dt*rhs(y)
  //   for (decltype(R.size()) i=0; i < R.size(); i++){
  //     R[i] = y[i] - yOld_[i] -dt_*R[i];
  //   }

  //   auto & jac = J.getNonConstRefToData();
  //   jac(0,0) = 1.0 - dt_ * jac(0,0);
  //   for (size_t i=1; i < J.rows(); ++i){
  //     jac(i,i-1) = - dt_ * jac(i,i-1);
  //     jac(i,i) = 1.0 - dt_ * jac(i,i);
  //   }
  // }



  // void doStepImpl(state_type & y_inout,
  // 		  time_type t,
  // 		  time_type dt )
  // {
  //   yOld_ = y_inout;
  //   dt_ = dt;
  //   t_ = t;

  //   // std::cout << "STEP IMP EULER " << t_ << std::endl;
  //   // state_type Vdotyold;
  //   // (*sysFunctor_).rescaleState(yOld_, Vdotyold);
  //   // for (int i=0; i < Vdotyold.size(); ++i)
  //   //   std::cout << std::setprecision(10) << Vdotyold[i]  << " ";
  //   // std::cout << std::endl;

  //   //algebra::newtonRaph<stepper_t,state_type,jacobian_type>(*this, y_inout);
  //   algebra::nonLinearLstsq<stepper_t,state_type,jacobian_type>(*this, y_inout);

  //   // std::cout << "AFTER SOLVE " << std::endl;
  //   // state_type Vdotynew;
  //   // (*sysFunctor_).rescaleState(y_inout, Vdotynew);
  //   // for (int i=0; i < Vdotynew.size(); ++i)
  //   //   std::cout << std::setprecision(10) << Vdotynew[i]  << " ";
  //   // std::cout << std::endl;

  //   // std::cout << "AFTER SOLVE at t = " << t_+dt_ << std::endl;
  //   // for (int i=0; i < y_inout.size(); ++i)
  //   //   std::cout <<  << y_inout[i]  << " ";
  //   // std::cout << std::endl;    
  // }
