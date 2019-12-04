namespace pressio{ namespace rom{ namespace wls{


template<typename api_tag, typename fom_type, typename wls_state_type,typename ... rest>
struct Specializer{
  using type = void;
};


/// Specializer for routunes that compute the hessian and gradient at the same time
template<
  typename fom_type, typename wls_state_type,
  typename decoder_type, typename hessian_gradient_pol_type
>
struct Specializer<pressio::rom::wls::WlsSystemHessianAndGradientApi<fom_type,wls_state_type,decoder_type, hessian_gradient_pol_type>, fom_type, wls_state_type, decoder_type, hessian_gradient_pol_type>
{
  using type = pressio::rom::wls::WlsSystemHessianAndGradientApi<fom_type, wls_state_type,decoder_type,hessian_gradient_pol_type>;
};


/// System type
template<typename api_type,typename fom_type,typename wls_state_type,typename ... rest>
using WlsSystem = typename pressio::rom::wls::Specializer<api_type, fom_type, wls_state_type, rest...>::type;





//struct DefaultApi{};
/*
template<
  typename fom_type, typename wls_state_type,
  typename residual_pol_type, typename jac_pol_type, typename decoder_type
>
struct Specializer<pressio::rom::wls::WlsSystemDefaultApi<fom_type,wls_state_type,residual_pol_type,jac_pol_type,decoder_type>, fom_type, wls_state_type, residual_pol_type, jac_pol_type,decoder_type>
{ 
//  static_assert( ::pressio::rom::meta::is_legitimate_residual_policy_for_wls<residual_pol_t>,
//     "You are trying to use Wls with default api but the residual policy passed \
is not admissible for this: maybe you have the wrong api? blas blas");
  using type = pressio::rom::wls::WlsSystemDefaultApi<fom_type, wls_state_type, residual_pol_type, jac_pol_type,decoder_type>;
};
*/

}}} 


