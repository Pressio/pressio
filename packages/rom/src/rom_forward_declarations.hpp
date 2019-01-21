
#ifndef ROM_FORWARD_DECLARATIONS_HPP_
#define ROM_FORWARD_DECLARATIONS_HPP_

#include "rom_ConfigDefs.hpp"
#include "../../ode/src/ode_enum.hpp"

namespace rompp{ namespace rom{

template <typename preconditionable,
	  typename enable = void>
class Preconditioned;


template<typename app_state_w_type,
	 typename app_res_w_type,
	 typename phi_op_type,
	 int maxNstates,
	 int maxNrhs,
	 typename A_type = void>
class RomLSPGResidualPolicy;


template<typename app_state_w_type,
	 typename app_jac_w_type,
	 typename phi_op_type,
	 int maxNstates,
	 typename A_type = void>
class RomLSPGJacobianPolicy;


template<typename operator_type,
	 typename enable = void>
class MultiVectorOperator;

template<::rompp::ode::ImplicitEnum,
	  typename app_state_w_type,
	  typename app_res_w_type,
	  typename phi_op_type,
	  typename ode_state_w_type,
	  typename ode_res_w_type,
	  typename A_type = phi_op_type>
class RomGalerkinImplicitResidualPolicy;


template<::rompp::ode::ImplicitEnum,
	  typename app_state_w_type,
	  typename app_jac_w_type,
	  typename phi_op_type,
	  typename ode_state_w_type,
	  typename ode_jac_w_type,
	  typename A_type = phi_op_type>
class RomGalerkinImplicitJacobianPolicy;


}} // end namespace rompp::rom
#endif
