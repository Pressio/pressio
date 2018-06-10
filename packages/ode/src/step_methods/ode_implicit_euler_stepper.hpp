
#ifndef ODE_IMPLICIT_EULER_STEPPER_HPP_
#define ODE_IMPLICIT_EULER_STEPPER_HPP_

#include "./base/ode_implicit_stepper_base.hpp"
//#include "vector/core_vector_traits.hpp"
//#include "least_squares/algebra_leastsquares.hpp"
//#include "algebra_newtonraphson.hpp"

namespace ode{

template<typename state_type, typename rhs_type, typename jacobian_type,
	 typename scalar_type, typename state_resizer_fnctor_type,
	 typename model_type, typename residual_policy_type, typename jacobian_policy_type,
	 typename nonlinearsolver_policy_type>
class implicitEulerStepper<state_type, rhs_type, jacobian_type,
			   scalar_type, state_resizer_fnctor_type, model_type,
			   residual_policy_type, jacobian_policy_type,
			   nonlinearsolver_policy_type,
			   typename 
			   std::enable_if< !std::is_same<state_type,void>::value &&
					   core::meta::is_default_constructible<state_resizer_fnctor_type>::value
					   >::type
			   >
  : public implicitStepperBase<implicitEulerStepper<state_type, rhs_type, jacobian_type,
						    scalar_type, state_resizer_fnctor_type,model_type,
						    residual_policy_type, jacobian_policy_type,
						    nonlinearsolver_policy_type> >
{
public :
  using stepper_t = implicitEulerStepper<state_type, rhs_type, jacobian_type, scalar_type,
					 state_resizer_fnctor_type, model_type, residual_policy_type,
					 jacobian_policy_type, nonlinearsolver_policy_type>;
  using stepper_base_t = implicitStepperBase<stepper_t>;
  

  // (de)constructors


  // residual policy is the standard one, but jacobian is not
  template < typename U = residual_policy_type, typename T = jacobian_policy_type>
  implicitEulerStepper(model_type & model,
		       U & res_policy_obj,
		       T & jac_policy_obj,
		       typename
		       std::enable_if< std::is_same<U,
		                                    ode::policy::implicitEulerStandardResidual<state_type, rhs_type,
		                                                                               model_type, details::time_type
		                                                                               >
		                                    >::value &&
		                      !std::is_same<T,
		                                    ode::policy::implicitEulerStandardJacobian<state_type, jacobian_type,
		                                                                               model_type, details::time_type
		                                                                               >
		                                    >::value
		                     >::type * = 0)
  : stepper_base_t(model, U(), jac_policy_obj)
  {}

  // residual policy is not standard , and jacobian is standard
  template < typename U = residual_policy_type, typename T = jacobian_policy_type>
  implicitEulerStepper(model_type & model,
		       U & res_policy_obj,
		       T & jac_policy_obj,
		       typename
		       std::enable_if< !std::is_same<U,
		                                    ode::policy::implicitEulerStandardResidual<state_type, rhs_type,
		                                                                               model_type, details::time_type
		                                                                               >
		                                    >::value &&
		                      std::is_same<T,
		                                    ode::policy::implicitEulerStandardJacobian<state_type, jacobian_type,
		                                                                               model_type, details::time_type
		                                                                               >
		                                    >::value
		                     >::type * = 0)
  : stepper_base_t(model, res_policy_obj, T() )
  {}

  // residual policy is standard , and jacobian is standard
  template < typename U = residual_policy_type, typename T = jacobian_policy_type>
  implicitEulerStepper(model_type & model,
		       U & res_policy_obj,
		       T & jac_policy_obj,
		       typename
		       std::enable_if< std::is_same<U,
		                                    ode::policy::implicitEulerStandardResidual<state_type, rhs_type,
		                                                                               model_type, details::time_type
		                                                                               >
		                                    >::value &&
		                      std::is_same<T,
		                                    ode::policy::implicitEulerStandardJacobian<state_type, jacobian_type,
		                                                                               model_type, details::time_type
		                                                                               >
		                                    >::value
		                     >::type * = 0)
  : stepper_base_t( model, U(), T() )
  {}


  // residual policy is not standard , and jacobian is not standard
  template < typename U = residual_policy_type, typename T = jacobian_policy_type>
  implicitEulerStepper(model_type & model,
		       U & res_policy_obj,
		       T & jac_policy_obj,
		       typename
		       std::enable_if< !std::is_same<U,
		                                    ode::policy::implicitEulerStandardResidual<state_type, rhs_type,
		                                                                               model_type, details::time_type
		                                                                               >
		                                    >::value &&
		                      !std::is_same<T,
		                                    ode::policy::implicitEulerStandardJacobian<state_type, jacobian_type,
		                                                                               model_type, details::time_type
		                                                                               >
		                                    >::value
		                     >::type * = 0)
  : stepper_base_t( model, res_policy_obj, jac_policy_obj )
  {}
  
  
  ~implicitEulerStepper(){}


  void doStepImpl(state_type & y_inout,
  		  ode::details::time_type t,
  		  ode::details::time_type dt )
  {
    y_nm1_ = y_inout;
    dt_ = dt;
    t_ = t;

    nonlinearsolver_policy_type::compute(*this, y_inout);
    //algebra::newtonRaph<stepper_t,state_type,jacobian_type>(*this, y_inout);
    //algebra::nonLinearLstsq<stepper_t,state_type,jacobian_type>(*this, y_inout);
  }

   
  void residualImpl(const state_type & y, state_type & R){
    this->res_policy_obj_.compute(y, y_nm1_, R, dt_);
  }

  void jacobianImpl(const state_type & y, jacobian_type & J){
    this->jac_policy_obj_.compute(y, J, dt_);
  }

private:
  ode::details::time_type t_;
  ode::details::time_type dt_;
  state_type y_nm1_;

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
  // 		  ode::details::time_type t,
  // 		  ode::details::time_type dt )
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
