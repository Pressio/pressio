
#ifndef ODE_ADAMS_MOULTON1_IMPLICIT_JACOBIAN_IMPL_HPP_
#define ODE_ADAMS_MOULTON1_IMPLICIT_JACOBIAN_IMPL_HPP_

#include "ode_ConfigDefs.hpp"
#include "matrix/core_matrix_traits.hpp"

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
void implicit_adams_moulton1_jacobian_impl(jacobian_type & jac,
					   jacobian_type & A,
					   time_type dt)
{
  using sc_t = typename core::details::traits<jacobian_type>::scalar_t;
  const sc_t c1 = static_cast<sc_t>(1)/2;
  jac.scale(-c1*dt);
  jac += A;
}

  
}//end namespace impl
}//end namespace ode
#endif 
