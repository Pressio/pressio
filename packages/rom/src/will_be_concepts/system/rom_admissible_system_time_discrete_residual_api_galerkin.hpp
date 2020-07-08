
#ifndef rom_admissible_system_time_discrete_residual_api_galerkin_hpp_
#define rom_admissible_system_time_discrete_residual_api_galerkin_hpp_

namespace pressio{ namespace rom{ namespace meta {

template<typename T, typename ...Args>
using admissible_system_time_discrete_residual_api_galerkin
 = admissible_system_time_discrete_residual_api_unsteady_lspg<T, Args...>;

}}} // namespace pressio::rom::meta
#endif
