
#ifndef ODE_IMPLICIT_EULER_STEPPER_HPP_
#define ODE_IMPLICIT_EULER_STEPPER_HPP_

#include "ode_implicit_stepper_base.hpp"
#include "vector/core_vector_traits.hpp"
#include "least_squares/algebra_leastsquares.hpp"
#include "algebra_newtonraphson.hpp"

namespace ode{

template<typename state_type,
	 typename rhs_type,
	 typename jacobian_type,
	 typename functor_type,
	 typename scalar_type,
	 typename state_resizer_fnctor_type
	 >
class implicitEulerStepper<state_type, rhs_type, jacobian_type, functor_type,
			   scalar_type, state_resizer_fnctor_type,
			   typename 
			   std::enable_if< !std::is_same<state_type,void>::value &&
					   core::meta::is_default_constructible<state_resizer_fnctor_type>::value
					   >::type
			   >
  : public implicit_stepper_base<implicitEulerStepper<state_type,rhs_type,jacobian_type,functor_type,
						      scalar_type, state_resizer_fnctor_type> >
{
public :
  using stepper_t = implicitEulerStepper<state_type,rhs_type,jacobian_type,functor_type,
					 scalar_type,state_resizer_fnctor_type>;
  using stepper_base_t = implicit_stepper_base<stepper_t>;
  
  // (de)constructors
  implicitEulerStepper(functor_type & sysFunctor)
    : stepper_base_t(), sysFunctor_(&sysFunctor)
  {}
  ~implicitEulerStepper(){}

  
  void operator()(const state_type & y, state_type & R, jacobian_type & J)
  {
    R.resize(y.size());
    assert(y.size() == yOld_.size());
    
    // compute RHS(y)
    (*sysFunctor_)(*y.view(), R.getNonConstRefToData(), J.getNonConstRefToData(), t_);

    // R = y - yOld - dt*rhs(y)
    for (decltype(y.size()) i=0; i < y.size(); i++){
      R[i] = y[i] - yOld_[i] -dt_*R[i];
    }

    auto & jac = J.getNonConstRefToData();
    jac(0,0) = 1.0 - dt_ * jac(0,0);
    for (size_t i=1; i < J.rows(); ++i){
      jac(i,i-1) = - dt_ * jac(i,i-1);
      jac(i,i) = 1.0 - dt_ * jac(i,i);
    }    
  }

  
  void doStepImpl(state_type & y_inout,
		  ode::details::time_type t,
		  ode::details::time_type dt )
  {
    yOld_ = y_inout;
    dt_ = dt;
    t_ = t;

    std::cout << "STEP IMP EULER " << t_ << std::endl;
    for (int i=0; i < yOld_.size(); ++i)
      std::cout << yOld_[i]  << " ";
    std::cout << std::endl;

    algebra::newtonRaph<stepper_t,state_type,jacobian_type>(*this, y_inout);
    //algebra::nonLinearLstsq<stepper_t,state_type,jacobian_type>(*this, y_inout);

    std::cout << "AFTER SOLVE at t = " << t_+dt_ << std::endl;
    for (int i=0; i < y_inout.size(); ++i)
      std::cout << std::setprecision(10) << y_inout[i]  << " ";
    std::cout << std::endl;
    
    // this->myResizer_(y_inout, RHS_);
    // //eval RHS
    // functor(y_inout, RHS_, t);    
    // // TODO: if possible, change to use native operations of the target data types
    // // out = in + dt * rhs    
    // for (decltype(y_inout.size()) i=0; i < y_inout.size(); i++){
    //   y_inout[i] += dt*RHS_[i];
    // }
  }

private:
  ode::details::time_type t_;
  ode::details::time_type dt_;
  rhs_type RHS_;
  functor_type * sysFunctor_;
  state_type yOld_;

}; //end class


}//end namespace
#endif 




  
  // #ifndef DOXYGEN_SKIP
    // #else
    // typedef explicit_stepper_base< euler< ... > , ... > stepper_base_type;
    // #endif
    // typedef typename stepper_base_type::state_type state_type;
    // typedef typename stepper_base_type::value_type value_type;
    // typedef typename stepper_base_type::deriv_type deriv_type;
    // typedef typename stepper_base_type::time_type time_type;
    // typedef typename stepper_base_type::algebra_type algebra_type;
    // typedef typename stepper_base_type::operations_type operations_type;
    // typedef typename stepper_base_type::resizer_type resizer_type;

    // #ifndef DOXYGEN_SKIP
    // typedef typename stepper_base_type::stepper_type stepper_type;
    // typedef typename stepper_base_type::wrapped_state_type wrapped_state_type;
    // typedef typename stepper_base_type::wrapped_deriv_type wrapped_deriv_type;
    // #endif 
