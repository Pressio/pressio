#include "./impl/wls_policies_impl.hpp"

namespace pressio{ namespace rom{ namespace wls{ 
template<typename fom_type,
	 typename decoder_type,
	 typename local_residual_policy_type, 
	 typename local_jacobian_policy_type,
	 typename ode_tag>
using hessian_gradient_policy = pressio::rom::wls::impl::hessian_gradient_policy<fom_type,decoder_type,local_residual_policy_type,local_jacobian_policy_type,ode_tag>;
}}}
