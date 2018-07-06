
#ifndef ODE_EULER_IMPLICIT_RESIDUAL_IMPL_HPP_
#define ODE_EULER_IMPLICIT_RESIDUAL_IMPL_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{
namespace impl{

template<typename state_type,
	 typename residual_type,
	 typename time_type>
void implicit_euler_residual_impl(const state_type & yn,
				  const state_type & ynm1,
				  residual_type & R,
				  time_type dt)
{
  // On input: R contains the application RHS, i.e. if
  //           dudt = f(x,u,...), R contains f(...)
  // On output, it contains the residual
  for (decltype(R.size()) i=0; i < R.size(); i++){
    R[i] = yn[i] - ynm1[i] - dt*R[i];
  }
  // std::cout << "\ndoImpl res euler " << std::endl;
  // std::cout << *R.data();
}


}//end namespace impl
}//end namespace ode
#endif 

