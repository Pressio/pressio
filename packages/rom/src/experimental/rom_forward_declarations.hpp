
#ifndef ROM_FORWARD_DECLARATIONS_HPP_
#define ROM_FORWARD_DECLARATIONS_HPP_

#include "../rom_ConfigDefs.hpp"
#include "../../../ode/src/ode_enum_steppers.hpp"

namespace rompp{ namespace rom{ namespace exp{
    
template<::rompp::ode::ImplicitSteppersEnum,
	 typename app_state_w_type,
	 typename app_res_w_type,
	 // basis operator type
	 typename phi_op_type,  
	 // weighting matrix by default is = void
	 typename A_type = void>
class RomLSPGResidualPolicy;

template<::rompp::ode::ImplicitSteppersEnum,
	 typename app_state_w_type,
	 typename app_jac_w_type,
	 // basis operator type
	 typename phi_op_type,  
	 // weighting matrix by default is = void
	 typename A_type = void>
class RomLSPGJacobianPolicy;
      
    
}}} // end namespace rompp::rom::exp
#endif
