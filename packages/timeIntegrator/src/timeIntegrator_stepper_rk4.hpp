
#ifndef TIMEINTEGRATOR_STEPPER_RK4_HPP
#define TIMEINTEGRATOR_STEPPER_RK4_HPP

#include "timeIntegrator_explicit_stepper_base.hpp"


template<
  class state_type,
  class deriv_type = state_type,
  class scalar_type = typename timeIntegrator::details::defaultTypes::scalar_t,
  class time_type = typename timeIntegrator::details::defaultTypes::scalar_t
>
class runge_kutta4
: public explicit_stepper_base<
  runge_kutta4< state_type, deriv_type, scalar_type, time_type>,
  4, state_type, deriv_type, scalar_type, time_type>
{
public :
  using stepper_base_t = explicit_stepper_base<
  runge_kutta4< state_type, deriv_type, scalar_type, time_type>,
  4, state_type, deriv_type, scalar_type, time_type>;

  // (de)constructors
  runge_kutta4() : stepper_base_t(){}
  virtual ~runge_kutta4(){}

    // methods
  template< class system_type>
  void step_impl( system_type & func,
		  const state_type & y_n,
		  const deriv_type & rhs,
		  state_type & y_next,
		  time_type t,
		  time_type dt )
  {
    cassert( in.size()==out.size() && in.size()==rhs.size() );
    static const scalar_type val1 = static_cast< scalar_type >( 1 );

    const time_type dh = dt / static_cast< scalar_type >( 2 );
    const time_type th = t + dh;

    // k1 = dt * rhs (rhs already calculated)

    // y_tmp = y_n + dh*dxdt
    for (int i=0; i<y_n.size(); i++){
      y_tmp[i] = y_n[i] + dh*rhs[i];
    }
    // k2 = dt * rhs2
    func( y_tmp, rhs2.data_, th );

    // y_tmp = y_n + dh*rhs2
    for (int i=0; i<y_n.size(); i++){
      y_tmp[i] = y_n[i] + dh*rhs2[i];
    }
    // k3 = dt * rhs3
    func( y_tmp, rhs3.data_, th );

    //y_tmp = y_n + dt*rhs3
    for (int i=0; i<y_n.size(); i++){
      y_tmp[i] = y_n[i] + dt*rhs3[i];
    }
    // k4 = dt * rhs4
    func( y_tmp, rhs4.data_, t + dt );

    //x += dt/6 * ( k1 + 2 * k2 + 2 * k3 + k4 )
    time_type dt6 = dt / static_cast< scalar_type >( 6.0 );
    time_type dt3 = dt / static_cast< scalar_type >( 3.0 )
    for (int i=0; i<in.size(); i++){
      y_next[i] = y_n[i] + dt6*rhs[i] + dt3*rhs2[i] + dt3*rhs3[i] + dt6*rhs4[i];
    }
  }//end step_impl
  
private:
    deriv_type rhs2;
    deriv_type rhs3;
    deriv_type rhs4;
    state_type y_tmp;
};
