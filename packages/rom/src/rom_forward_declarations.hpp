
#ifndef ROM_FORWARD_DECLARATIONS_HPP_
#define ROM_FORWARD_DECLARATIONS_HPP_

#include "rom_ConfigDefs.hpp"
#include "../../ode/src/ode_enum_steppers.hpp"

namespace rompp{ namespace rom{
    
template<::rompp::ode::ImplicitSteppersEnum,
	  typename app_state_w_type,
	  typename app_res_w_type,
	  typename phi_op_type,  
	  typename A_type = void>
class RomLSPGResidualPolicy;

template<::rompp::ode::ImplicitSteppersEnum,
	  typename app_state_w_type,
	  typename app_jac_w_type,
	  typename phi_op_type,  
	  typename A_type = void>
class RomLSPGJacobianPolicy;



template<::rompp::ode::ImplicitSteppersEnum,
	  typename app_state_w_type,
	  typename app_res_w_type,
	  typename phi_op_type,
	  typename ode_state_w_type,
	  typename ode_res_w_type,
	  typename A_type = phi_op_type>
class RomGalerkinImplicitResidualPolicy;

template<::rompp::ode::ImplicitSteppersEnum,
	  typename app_state_w_type,
	  typename app_jac_w_type,
	  typename phi_op_type,  
	  typename ode_state_w_type,
	  typename ode_jac_w_type,
	  typename A_type = phi_op_type>
class RomGalerkinImplicitJacobianPolicy;
    
    
}} // end namespace rompp::rom
#endif
