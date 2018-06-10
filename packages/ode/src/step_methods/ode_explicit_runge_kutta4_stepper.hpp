
#ifndef ODE_EXPLICIT_RUNGEKUTTA4_STEPPER_HPP_
#define ODE_EXPLICIT_RUNGEKUTTA4_STEPPER_HPP_

#include "./base/ode_explicit_stepper_base.hpp"


namespace ode{

template<typename state_type,
	 typename rhs_type,
	 typename scalar_type,
	 typename state_resizer_fnctor_type,
	 typename model_type,
	 typename residual_policy_type
	 >
class explicitRungeKutta4Stepper<state_type, rhs_type, scalar_type,
			 state_resizer_fnctor_type, model_type,
			 residual_policy_type, 
			 typename
			 std::enable_if<!std::is_void<state_type>::value &&
					core::meta::is_default_constructible<state_resizer_fnctor_type>::value
					>::type
			 >
  : public explicitStepperBase< explicitRungeKutta4Stepper<state_type,rhs_type,scalar_type,
						   state_resizer_fnctor_type, model_type,
						   residual_policy_type
						   >
				>
{
public :
  using stepper_t = explicitRungeKutta4Stepper<state_type,rhs_type,scalar_type,
				       state_resizer_fnctor_type, model_type,
				       residual_policy_type>;
  using stepper_base_t = explicitStepperBase<stepper_t>;
  
  //(de)constructors
  template < typename U = residual_policy_type>
  explicitRungeKutta4Stepper(model_type & model,
  		       U & res_policy_obj, 
  		       typename std::enable_if< !std::is_same<U,
  		       ode::policy::explicitRungeKutta4StandardResidual<state_type, rhs_type,
		                                                model_type, details::time_type>
		       >::value >::type * = 0)
    : stepper_base_t(model, res_policy_obj)
  {}

  template < typename U = residual_policy_type>
  explicitRungeKutta4Stepper(model_type & model,
  		       typename
  		       std::enable_if< std::is_same<U,
  		       ode::policy::explicitRungeKutta4StandardResidual<state_type, rhs_type,
		                                                model_type, details::time_type>
  		       >::value >::type * = 0)
    : stepper_base_t( model, U() ) 
  {}

  ~explicitRungeKutta4Stepper()
  {}

  
  // methods
  void doStepImpl(state_type & y_inout,
		  ode::details::time_type t,
		  ode::details::time_type dt )
  {
    static const scalar_type val1 = static_cast< scalar_type >( 1 );
    const ode::details::time_type dt_half = dt / static_cast< scalar_type >(2);
    const ode::details::time_type t_phalf = t + dt_half;
    
    this->myResizer_(y_inout, rhs1_);
    this->myResizer_(y_inout, rhs2_);
    this->myResizer_(y_inout, rhs3_);
    this->myResizer_(y_inout, rhs4_);
    this->myResizer_(y_inout, y_tmp_);
    
    // ----------
    // stage 1: 
    // ----------
    // rhs1_(y_n,t)
    this->residual_policy_obj_.compute(y_inout, rhs1_, *(this->model_), t);
    // y_tmp_ = y_n + rhs1_*dt/2
    for (size_t i=0; i<y_inout.size(); i++){
      y_tmp_[i] = y_inout[i] + dt_half*rhs1_[i];
    }

    // ----------
    // stage 2: 
    // ----------
    // rhs2_
    this->residual_policy_obj_.compute(y_tmp_, rhs2_, *(this->model_), t_phalf);
    // y_tmp_ = y_n + rhs2_*dt/2
    for (size_t i=0; i<y_inout.size(); i++){
      y_tmp_[i] = y_inout[i] + dt_half*rhs2_[i];
    }

    // ----------
    // stage 3: 
    // ----------
    // rhs3_
    this->residual_policy_obj_.compute(y_tmp_, rhs3_, *(this->model_), t_phalf);
    //y_tmp_ = y_n + rhs3_*dt/2
    for (size_t i=0; i<y_inout.size(); i++){
      y_tmp_[i] = y_inout[i] + dt*rhs3_[i];
    }

    // ----------
    // stage 4: 
    // ----------
    // rhs4_
    this->residual_policy_obj_.compute(y_tmp_, rhs4_, *(this->model_), t + dt);
    //x += dt/6 * ( k1 + 2 * k2 + 2 * k3 + k4 )
    ode::details::time_type dt6 = dt / static_cast< scalar_type >( 6.0 );
    ode::details::time_type dt3 = dt / static_cast< scalar_type >( 3.0 );
    for (size_t i=0; i < y_inout.size(); i++){
      y_inout[i] = y_inout[i] + dt6*rhs1_[i] + dt3*rhs2_[i] + dt3*rhs3_[i] + dt6*rhs4_[i];
    }
    
  }//end doStepImpl

private:
  rhs_type rhs1_;
  rhs_type rhs2_;
  rhs_type rhs3_;
  rhs_type rhs4_;
  state_type y_tmp_;

  // inherited from base: model_, myResizer_, residual_policy_t
  
}; //end class
  
}//end namespace
#endif 
