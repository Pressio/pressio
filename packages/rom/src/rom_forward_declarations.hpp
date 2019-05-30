
#ifndef ROM_FORWARD_DECLARATIONS_HPP_
#define ROM_FORWARD_DECLARATIONS_HPP_

#include "../../ode/src/ode_enum.hpp"

namespace rompp{ namespace rom{

/* forward declare all decorators */
namespace decorator{

template <
  typename preconditionable,
  typename enable = void>
class Preconditioned;

template <
  typename maskable,
  typename enable = void>
class Masked;

}// namespace rompp::rom::decorator
//-------------------------------------------

namespace policy{

template <bool is_steady_problem>
struct EvaluateFomRhsDefault;

template <bool is_steady_problem>
struct ApplyFomJacobianDefault;

}// namespace rompp::rom::policy
//-------------------------------------------

/* operators */
template<
  typename wrapped_type,
  typename enable = void
  >
class MultiVectorOperator;

template<
  typename wrapped_type,
  typename enable = void
  >
class MatrixOperator;
//-------------------------------------------

/* LSPG policies */
template <
  typename fom_states_data_t,
  typename fom_rhs_data_t,
  typename fom_rhs_eval_policy,
  typename td_ud_ops = void
  >
class LSPGResidualPolicy;

template<
  typename fom_states_data_t,
  typename jac_type,
  typename fom_apply_jac_policy,
  typename decoder_t,
  typename td_ud_ops = void
  >
class LSPGJacobianPolicy;
//-------------------------------------------

/* steady LSPG policies */
template <
  typename fom_states_data_t,
  typename fom_rhs_data_t,
  typename fom_rhs_eval_policy
  >
class LSPGSteadyResidualPolicy;

template<
  typename fom_states_data_t,
  typename jac_type,
  typename fom_apply_jac_policy,
  typename decoder_t
  >
class LSPGSteadyJacobianPolicy;
//-------------------------------------------

/* Explicit Galerkin policies */
template <
  typename fom_states_data_t,
  typename fom_rhs_data_t,
  typename decoder_jac_t
  >
class DefaultGalerkinExplicitResidualPolicy;
//-------------------------------------------

template <
  typename type_generator_t,
  typename enable = void
  >
struct GalerkinStepperObjectGenerator;

template <
  typename type_generator_t,
  typename enable = void
  >
struct LSPGUnsteadyProblemGenerator;

template <
  typename type_generator_t,
  typename enable = void
  >
struct LSPGSteadyProblemGenerator;
//-------------------------------------------

template<
  typename app_type,
  typename lspg_state_type,
  typename lspg_residual_type,
  typename lspg_jacobian_type,
  typename residual_policy_type,
  typename jacobian_policy_type,
  typename enable = void
  >
class LSPGSteadySystem;

}} // end namespace rompp::rom
#endif
