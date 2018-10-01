
#ifndef ODE_RESIDUAL_IMPL_HPP_
#define ODE_RESIDUAL_IMPL_HPP_

#include "ode_ConfigDefs.hpp"

namespace rompp{
namespace ode{
namespace impl{

template<typename state_type, typename time_type>
void implicit_euler_time_discrete_residual(const state_type & yn,
					   const state_type & ynm1,
					   state_type & R,
					   time_type dt){
  // On input: R contains the application RHS, i.e. if
  //           dudt = f(x,u,...), R contains f(...)
  // On output, it contains the residual
  R = yn - ynm1 - dt*R;
  // for (decltype(R.size()) i=0; i < R.size(); i++){
  //   R[i] = yn[i] - ynm1[i] - dt*R[i];
  // }
}


// template<typename state_type,
// 	 typename residual_type,
// 	 typename time_type>
// void implicit_bdf2_time_discrete_residual(const state_type & yn,
// 					  const state_type & ynm1,
// 					  const state_type & ynm2,
// 					  residual_type & R,
// 					  time_type dt)
// {
//   using sc_t = typename core::details::traits<state_type>::scalar_t;
//   const sc_t c1 = static_cast<sc_t>(4)/3;
//   const sc_t c2 = static_cast<sc_t>(1)/3;
//   const sc_t c3 = static_cast<sc_t>(2)/3;
  
//   // On input: R contains the application RHS, i.e. if
//   //           dudt = f(x,u,...), R contains f(...)
//   // On output, it contains the residual
//   for (decltype(R.size()) i=0; i < R.size(); i++){
//     R[i] = yn[i] - c1*ynm1[i] +c2*ynm2[i] - c3*dt*R[i];
//   }
//   // std::cout << "\ndoImpl res euler " << std::endl;
//   // std::cout << *R.data();
// }


}//end namespace impl
}//end namespace ode
}//end namespace rompp
#endif 

