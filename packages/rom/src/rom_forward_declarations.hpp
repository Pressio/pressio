
#ifndef ROM_FORWARD_DECLARATIONS_HPP_
#define ROM_FORWARD_DECLARATIONS_HPP_

#include "../../ode/src/ode_enum.hpp"

namespace rompp{ namespace rom{

template <typename fom_states_data_t,
	  typename fom_rhs_data_t,
	  typename fom_rhs_eval_policy>
class LSPGResidualPolicy;

template< typename fom_states_data_t,
	  typename jac_type,
	  typename fom_apply_jac_policy>
class LSPGJacobianPolicy;

namespace decorator{

// decorators usable to decorate policies above
template <typename preconditionable,
	  typename enable = void>
class Preconditioned;

}//end namespace rompp::rom::decorator
//-------------------------------------------


template<typename operator_type,
	 typename enable = void>
class MultiVectorOperator;

//-------------------------------------------

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
