
#ifndef ODE_EULER_IMPLICIT_IMPL_HPP_
#define ODE_EULER_IMPLICIT_IMPL_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{
namespace impl{

template<typename state_type, typename residual_type, typename time_type>
void implicit_euler_residual_impl(const state_type & yn, const state_type & ynm1,
				  residual_type & R, time_type dt)
{
  // on input, R contains the application RHS, i.e. if
  // dudt = f(x,u,...), R contains f(...)
  for (decltype(R.size()) i=0; i < R.size(); i++){
    R[i] = yn[i] - ynm1[i] - dt*R[i];
  }
}
//----------------------------------------------------

template<typename jacobian_type, typename time_type>
void implicit_euler_jacobian_impl(jacobian_type & jac,
				  time_type dt)
{
  // obviously this needs to be fixed to use operator []
  // of the type coming in 
  // auto & jac = J.getNonConstRefToData();
 jac[0,0] = 1.0 - dt * jac[0,0];
 for (size_t i=1; i < jac.rows(); ++i){
   jac[i,i-1] = - dt * jac[i,i-1];
   jac[i,i] = 1.0 - dt * jac[i,i];
 }
}
//----------------------------------------------------



}//end namespace impl
}//end namespace ode
#endif 

