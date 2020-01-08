#include "./impl/wls_hess_and_grad_policies_impl.hpp"

namespace pressio{ namespace rom{ namespace wls{ 
template<typename fom_type,
	 typename decoder_type,
	 typename ode_tag>
using hessian_gradient_policy = pressio::rom::wls::impl::hessian_gradient_policy<fom_type,decoder_type,ode_tag>;
}}}
