
#ifndef ROM_WLS_HESSIAN_GRADIENT_SYSTEM_API_HPP_
#define ROM_WLS_HESSIAN_GRADIENT_SYSTEM_API_HPP_

#include "rom_wls_hessian_gradient_system_api_impl.hpp"

namespace pressio{ namespace rom{ namespace wls{

template<typename ... Args>
using SystemHessianGradientApi = impl::SystemHessianGradientApi<Args...>;

template<typename ... Args>
using SystemHessianAndGradientApi = impl::SystemHessianGradientApi<Args...>;

}}}
#endif
