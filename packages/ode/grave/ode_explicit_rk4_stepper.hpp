
#ifndef ODE_RK4_STEPPER_HPP_
#define ODE_RK4_STEPPER_HPP_

#include "ode_explicit_stepper_base.hpp"


namespace ode{

template<typename state_type,
	 typename rhs_type,
	 typename scalar_type,
	 typename state_resizer_fnctor_type
	 >
class rungeKutta4Stepper<state_type, rhs_type, scalar_type, state_resizer_fnctor_type,
			 typename
			 std::enable_if< !std::is_same<state_type,void>::value &&
					 core::meta::is_default_constructible<state_resizer_fnctor_type>::value
					 >::type
			 >
  : public explicit_stepper_base< rungeKutta4Stepper<state_type,rhs_type,
						     scalar_type,state_resizer_fnctor_type> >
{
public :
  using stepper_t = rungeKutta4Stepper<state_type,rhs_type,scalar_type,state_resizer_fnctor_type>;
  using stepper_base_t = explicit_stepper_base<stepper_t>;
  using time_t = ode::details::time_type;
  
  // (de)constructors
  rungeKutta4Stepper() : stepper_base_t(){}
  virtual ~rungeKutta4Stepper(){}

    // methods
  template< typename functor_type>
  void do_step_impl( functor_type & functor,
		     state_type & y_inout,
		     time_t t,
		     time_t dt )
  {
    static const scalar_type val1 = static_cast< scalar_type >( 1 );
    const time_t dt_half = dt / static_cast< scalar_type >(2);
    const time_t t_phalf = t + dt_half;

    myResizer_(y_inout, rhs1_);
    myResizer_(y_inout, rhs2_);
    myResizer_(y_inout, rhs3_);
    myResizer_(y_inout, rhs4_);

    // ----------
    // stage 1: 
    // ----------
    // rhs1_(y_n,t)
    functor( y_inout, rhs1_, t );
    // y_tmp_ = y_n + rhs1_*dt/2
    for (int i=0; i<y_inout.size(); i++){
      y_tmp_[i] = y_inout[i] + dt_half*rhs1_[i];
    }
    // ----------
    // stage 2: 
    // ----------
    // rhs2_
    functor( y_tmp_, rhs2_, t_phalf );
    // y_tmp_ = y_n + rhs2_*dt/2
    for (int i=0; i<y_inout.size(); i++){
      y_tmp_[i] = y_inout[i] + dt_half*rhs2_[i];
    }
    // ----------
    // stage 3: 
    // ----------
    // rhs3_
    functor( y_tmp_, rhs3_, t_phalf );
    //y_tmp_ = y_n + rhs3_*dt/2
    for (int i=0; i<y_inout.size(); i++){
      y_tmp_[i] = y_inout[i] + dt*rhs3_[i];
    }

    // ----------
    // stage 4: 
    // ----------
    // rhs4_
    functor( y_tmp_, rhs4_, t + dt );

    //x += dt/6 * ( k1 + 2 * k2 + 2 * k3 + k4 )
    time_t dt6 = dt / static_cast< scalar_type >( 6.0 );
    time_t dt3 = dt / static_cast< scalar_type >( 3.0 );
    for (int i=0; i < y_inout.size(); i++){
      y_inout[i] = y_inout[i] + dt6*rhs1_[i] + dt3*rhs2_[i] + dt3*rhs3_[i] + dt6*rhs4_[i];
    }
  }//end step_impl
  
private:
    rhs_type rhs1_;
    rhs_type rhs2_;
    rhs_type rhs3_;
    rhs_type rhs4_;
    state_type y_tmp_;
};

}//end namespace
#endif
