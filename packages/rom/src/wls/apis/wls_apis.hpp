#include "./impl/wls_apis_impl.hpp"
// basic alias to impl namespace
namespace pressio{ namespace rom{ namespace wls{ 
template<typename fom_type, 
	typename wls_state_type, 
	typename decoder_type,
  typename ode_type, 
	typename hessian_type>
using  WlsSystemHessianAndGradientApi = pressio::rom::wls::impl::WlsSystemHessianAndGradientApi<fom_type,wls_state_type,decoder_type,ode_type,hessian_type>;
}}}





