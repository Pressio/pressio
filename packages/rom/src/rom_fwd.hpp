
#ifndef ROM_FORWARD_DECLARATIONS_HPP_
#define ROM_FORWARD_DECLARATIONS_HPP_

#include "../../ode/src/ode_enum.hpp"

namespace pressio{ namespace rom{

/* decorators */
namespace decorator{

template <
  typename preconditionable,
  typename enable = void>
class Preconditioned;

template <
  typename maskable,
  typename enable = void>
class Masked;

}// namespace pressio::rom::decorator
//---------------------------------

namespace policy{

template <bool is_steady_problem>
struct EvaluateFomVelocityDefault;

template <bool is_steady_problem>
struct ApplyFomJacobianDefault;

}// namespace pressio::rom::policy
//---------------------------------

#ifdef HAVE_PYBIND11
template <
  typename matrix_type,
  typename ops_t,
  typename enable = void
  >
struct PyLinearDecoder;
#endif

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


/* ------------------
 * steady LSPG
 ------------------ */

/* policies */
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

/* problem */
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

template <
  typename type_generator_t,
  typename enable = void
  >
struct LSPGSteadyProblemGenerator;


/* ------------------
 * UNsteady LSPG
 ------------------ */

/* policies */
template <
  typename fom_states_data_t,
  typename fom_rhs_data_t,
  typename fom_rhs_eval_policy,
  typename ud_ops = void
  >
class LSPGResidualPolicy;

template<
  typename fom_states_data_t,
  typename jac_type,
  typename fom_apply_jac_policy,
  typename decoder_t,
  typename ud_ops = void
  >
class LSPGJacobianPolicy;

template <
  typename type_generator_t,
  typename enable = void
  >
struct LSPGUnsteadyProblemGenerator;
//-----------------------------------

/* Explicit Galerkin policies */
template <
  typename fom_states_data_t,
  typename fom_rhs_data_t,
  typename decoder_jac_t
  >
class DefaultGalerkinExplicitVelocityPolicy;

template <
  typename type_generator_t,
  typename enable = void
  >
struct GalerkinProblemGenerator;


}} // end namespace pressio::rom
#endif
