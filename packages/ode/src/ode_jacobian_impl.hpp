
#ifndef ODE_JACOBIAN_IMPL_HPP_
#define ODE_JACOBIAN_IMPL_HPP_

#include "ode_ConfigDefs.hpp"

namespace rompp{
namespace ode{
namespace impl{

template <typename jacobian_type, typename time_type,
	  typename
	  std::enable_if<
	    core::details::traits<jacobian_type>::isMatrix==1 &&
	    core::details::traits<jacobian_type>::isSparse==1 &&
	    core::details::traits<jacobian_type>::isEigen==1
	    >::type * = nullptr
	  >
void implicit_euler_time_discrete_jacobian(jacobian_type & jac,
					   time_type dt)
{
  jac.scale(-dt);
  // auwo II(jac);
  // II.setIdentity();
  //jac += II;  
  jac.addToDiagonal(static_cast<time_type>(1));
}


//need to add also other overloads  

}//end namespace impl
}//end namespace ode
}//end namespace rompp
#endif 

