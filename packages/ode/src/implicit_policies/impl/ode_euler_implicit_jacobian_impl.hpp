
#ifndef ODE_EULER_IMPLICIT_JACOBIAN_IMPL_HPP_
#define ODE_EULER_IMPLICIT_JACOBIAN_IMPL_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{
namespace impl{
 
//******************************************************
// enable for SPARSE serial eigen matrix
//******************************************************
template <typename jacobian_type,
	  typename time_type,
	  typename
	  std::enable_if<
	    core::details::traits<jacobian_type>::isMatrix==1 &&
	    core::details::traits<jacobian_type>::isSparse==1 &&
	    core::details::traits<jacobian_type>::isEigen==1
	    >::type * = nullptr
	  >
void implicit_euler_jacobian_impl(jacobian_type & jac,
				  jacobian_type & A,
				  time_type dt)
{
  jac.scale(-dt);
  jac += A;
}

  
}//end namespace impl
}//end namespace ode
#endif 
