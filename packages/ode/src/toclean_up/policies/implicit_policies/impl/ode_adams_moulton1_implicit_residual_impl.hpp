
#ifndef ODE_ADAMS_MOULTON1_IMPLICIT_RESIDUAL_IMPL_HPP_
#define ODE_ADAMS_MOULTON1_IMPLICIT_RESIDUAL_IMPL_HPP_

#include "ode_ConfigDefs.hpp"
#include "vector/core_vector_traits.hpp"

namespace ode{
namespace impl{

template<typename state_type,
	 typename residual_type,
	 typename time_type>
void implicit_adams_moulton1_residual_impl(const state_type & yn,
				 const state_type & ynm1,
				 const state_type & RHSnm1,
				 residual_type & R,
				 time_type dt)
{
  using sc_t = typename core::details::traits<state_type>::scalar_t;
  const sc_t c1 = static_cast<sc_t>(1)/2;
  
  // On input: R contains the application RHS, i.e. if
  //           dudt = f(x,u,...), R contains f(...)
  // On output, it contains the residual
  for (decltype(R.size()) i=0; i < R.size(); i++){
    R[i] = yn[i] - ynm1[i] - c1*dt*R[i] - c1*dt*RHSnm1[i];
  }
  // std::cout << "\ndoImpl res euler " << std::endl;
  // std::cout << *R.data();
}


}//end namespace impl
}//end namespace ode
#endif 

