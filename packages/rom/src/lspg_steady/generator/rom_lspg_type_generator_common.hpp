
#ifndef ROM_LSPG_STEADY_TYPE_GENERATOR_COMMON_HPP_
#define ROM_LSPG_STEADY_TYPE_GENERATOR_COMMON_HPP_

#include "../../rom_ConfigDefs.hpp"
#include "../../rom_forward_declarations.hpp"
#include "../../rom_data_fom_rhs.hpp"
#include "../../rom_data_fom_states.hpp"
#include "../../rom_reconstructor_fom_state.hpp"
#include "../../policies/rom_evaluate_fom_rhs_policy.hpp"
#include "../../policies/rom_apply_fom_jacobian_policy.hpp"
#include "../../../../ode/src/ode_forward_declarations.hpp"

namespace rompp{ namespace rom{

template <
  typename fom_type,
  typename decoder_type,
  typename lspg_state_type
  >
struct LSPGSteadyCommonTypes{
  // these are native types of the full-order model (fom)
  using fom_t			= fom_type;
  using scalar_t		= typename fom_t::scalar_type;
  using fom_state_t		= typename fom_t::state_type;
  using fom_rhs_t		= typename fom_t::residual_type;

  // fom wrapper types
  using fom_state_w_t		= ::rompp::core::Vector<fom_state_t>;
  using fom_rhs_w_t		= ::rompp::core::Vector<fom_rhs_t>;

  // rom state type (passed in)
  using lspg_state_t		= lspg_state_type;

  // for LSPG, the rom residual type = core::wrapper of application rhs
  // i.e. the wrapped fom rhs type
  using lspg_residual_t		= fom_rhs_w_t;

  // decoder types (passed in)
  using decoder_t		= decoder_type;
  using decoder_jac_t		= typename decoder_t::jacobian_t;

  // fom state reconstructor type
  using fom_state_reconstr_t	= FomStateReconstructor<
    fom_state_w_t, decoder_t>;

  // class type holding fom states data
  using fom_states_data = ::rompp::rom::FomStatesData<
	fom_state_w_t, 1, fom_state_reconstr_t>;

  // class type holding fom rhs data
  using fom_rhs_data = ::rompp::rom::FomRhsData<fom_rhs_w_t>;
};

}}//end  namespace rompp::rom
#endif
