
#ifndef ROM_FORWARD_DECLARATIONS_HPP_
#define ROM_FORWARD_DECLARATIONS_HPP_

#include "../../ode/src/ode_enum.hpp"

namespace rompp{ namespace rom{

/* forward declare all decorators */
namespace decorator{

template <typename preconditionable,
	  typename enable = void>
class Preconditioned;

template <typename maskable,
	  typename enable = void>
class Masked;

}// namespace rompp::rom::decorator
//-------------------------------------------

/* LSPG policies */
template <typename fom_states_data_t,
	  typename fom_rhs_data_t,
	  typename fom_rhs_eval_policy>
class LSPGResidualPolicy;

template< typename fom_states_data_t,
	  typename jac_type,
	  typename fom_apply_jac_policy>
class LSPGJacobianPolicy;
//-------------------------------------------

/* operators */
template<typename wrapped_type,
	 typename enable = void>
class MultiVectorOperator;

template<typename wrapped_type,
	 typename enable = void>
class MatrixOperator;
//-------------------------------------------


template <typename lspg_type_generator_t, typename enable = void>
struct StepperObjectGenerator;
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
