
#ifndef ROM_LSPG_STEADY_TYPE_GENERATOR_DEFAULT_HPP_
#define ROM_LSPG_STEADY_TYPE_GENERATOR_DEFAULT_HPP_

#include "../../rom_lspg_type_generator_common.hpp"

namespace rompp{ namespace rom{

template <
  typename fom_type,
  typename decoder_type,
  typename lspg_state_type
  >
struct DefaultLSPGSteadyTypeGenerator
  : LSPGCommonTypes<
  fom_type, decoder_type, lspg_state_type
  >{

  using this_t = DefaultLSPGSteadyTypeGenerator
    <fom_type, decoder_type, lspg_state_type>;

  using base_t = LSPGCommonTypes
    <fom_type, decoder_type, lspg_state_type>;

  using typename base_t::fom_t;
  using typename base_t::scalar_t;
  using typename base_t::fom_native_state_t;
  using typename base_t::fom_state_t;
  using typename base_t::fom_rhs_t;
  using typename base_t::lspg_state_t;
  using typename base_t::lspg_residual_t;
  using typename base_t::decoder_t;
  using typename base_t::decoder_jac_t;
  using typename base_t::fom_state_reconstr_t;
  using typename base_t::fom_states_data;
  using typename base_t::fom_rhs_data;

  static constexpr bool is_steady = true;

  /* lspg_matrix_t is type of J*decoder_jac_t (in the most basic case) where
   * * J is the jacobian of the fom rhs
   * * decoder_jac_t is the type of the decoder jacobian
   * In more complex cases, we might have (something)*J*decoder_jac_t,
   * where (something) is product of few matrices.
   * For now, set lspg_matrix_t to be of same type as decoder_jac_t
   * if phi is MV<>, then lspg_matrix_t = containers::MV<>
   * if phi is Matrix<>, then we have containers::Matrix<>
   * not a bad assumption since all matrices are left-applied to decoder_jac_t
   */
  using lspg_matrix_t		= decoder_jac_t;

  // policy for evaluating the rhs of the fom object (<true> for steady overload)
  using fom_eval_rhs_policy_t	= ::rompp::rom::policy::EvaluateFomRhsDefault<this_t::is_steady>;

  // policy for left multiplying the fom jacobian with decoder_jac_t
  // possibly involving other stuff like explained above (<true> for steady overload
  using fom_apply_jac_policy_t	= ::rompp::rom::policy::ApplyFomJacobianDefault<this_t::is_steady>;

  // Policy defining how to compute the LSPG residual
  using lspg_residual_policy_t	= ::rompp::rom::LSPGSteadyResidualPolicy<
	fom_states_data, fom_rhs_data, fom_eval_rhs_policy_t>;

  // policy defining how to compute the LSPG jacobian
  using lspg_jacobian_policy_t	= ::rompp::rom::LSPGSteadyJacobianPolicy<
    fom_states_data, lspg_matrix_t, fom_apply_jac_policy_t, decoder_t>;

  // system's type
  using lspg_system_t		= ::rompp::rom::LSPGSteadySystem<
    fom_t, lspg_state_type, lspg_residual_t, lspg_matrix_t,
    lspg_residual_policy_t, lspg_jacobian_policy_t>;

};//end class

}}//end  namespace rompp::rom
#endif
