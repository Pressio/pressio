
#ifndef ODE_JACOBIAN_IMPL_HPP_
#define ODE_JACOBIAN_IMPL_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../../../core/src/matrix/core_matrix_meta.hpp"

namespace rompp{ namespace ode{ namespace impl{

template <typename jacobian_type,
	  typename time_type,
	  typename std::enable_if<
	    core::meta::is_eigen_sparse_matrix_wrapper<jacobian_type>::value
	    >::type * = nullptr
	  >
void implicit_euler_time_discrete_jacobian(jacobian_type & jac,
					   time_type dt){
  
  jac.scale(-dt);
  jac.addToDiagonal(static_cast<time_type>(1));
}

template <typename jacobian_type,
	  typename time_type,
	  typename std::enable_if<
	    core::meta::is_eigen_dense_matrix_wrapper<jacobian_type>::value
	    >::type * = nullptr
	  >
void implicit_euler_time_discrete_jacobian(jacobian_type & jac,
					   time_type dt){
  
  jac.scale(-dt);
  jac.addToDiagonal(static_cast<time_type>(1));
}

      
//need to add also other overloads  

      
}}}//end namespace rompp::ode::impl
#endif 
