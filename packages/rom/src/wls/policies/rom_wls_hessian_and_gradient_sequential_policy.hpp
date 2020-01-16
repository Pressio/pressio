
#ifndef ROM_WLS_HESSIAN_AND_GRADIENT_SEQUENTIAL_POLICY_HPP_
#define ROM_WLS_HESSIAN_AND_GRADIENT_SEQUENTIAL_POLICY_HPP_

#include "rom_wls_hessian_and_gradient_sequential_policy_impl.hpp"

namespace pressio{ namespace rom{ namespace wls{

template<typename ... Args>
using HessianGradientSequentialPolicy = ::pressio::rom::wls::impl::HessianGradientSequentialPolicy<Args...>;

}}} // end namespace pressio::rom::wls
#endif
