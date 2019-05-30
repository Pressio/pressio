
#ifndef ROM_LSPG_TYPE_GENERATOR_COMMON_HPP_
#define ROM_LSPG_TYPE_GENERATOR_COMMON_HPP_

#include "../rom_ConfigDefs.hpp"
#include "../rom_forward_declarations.hpp"
#include "../rom_data_fom_rhs.hpp"
#include "../rom_data_fom_states.hpp"
#include "../rom_reconstructor_fom_state.hpp"
#include "../policies/rom_evaluate_fom_rhs_policy.hpp"
#include "../policies/rom_apply_fom_jacobian_policy.hpp"
#include "../../../ode/src/ode_forward_declarations.hpp"

namespace rompp{ namespace rom{

template <ode::ImplicitEnum odeName>
struct statesStorageCapacityHelper{
  static constexpr int maxAuxStates_ = 1;
};

template <>
struct statesStorageCapacityHelper<ode::ImplicitEnum::Euler>{
  static constexpr int maxAuxStates_ = 1;
};

template <>
struct statesStorageCapacityHelper<ode::ImplicitEnum::BDF2>{
  static constexpr int maxAuxStates_ = 2;
};
//-------------------------------------------------------


template <ode::ImplicitEnum odeName, typename ... Rest>
struct auxStepperHelper{
  using type = void;
};

template <
  typename lspg_state_type,
  typename lspg_residual_t,
  typename lspg_matrix_t,
  typename fom_type,
  typename lspg_residual_policy_t,
  typename lspg_jacobian_policy_t
  >
struct auxStepperHelper<
  ode::ImplicitEnum::BDF2,
  lspg_state_type,
  lspg_residual_t,
  lspg_matrix_t,
  fom_type,
  lspg_residual_policy_t,
  lspg_jacobian_policy_t
  >{
  using type = ode::ImplicitStepper<
    ode::ImplicitEnum::Euler, lspg_state_type, lspg_residual_t,
    lspg_matrix_t, fom_type, void, lspg_residual_policy_t,
    lspg_jacobian_policy_t>;
};
//-------------------------------------------------------


template <
  typename fom_type,
  typename decoder_type,
  typename lspg_state_type,
  ode::ImplicitEnum odeName = ode::ImplicitEnum::Undefined,
  typename time_discrete_ud_ops = void
  >
struct LSPGCommonTypes{
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

  // max num of states needed for time integration.
  // this is deduced based on the integrator, i.e. odeName
  static constexpr auto maxAuxStates =
    statesStorageCapacityHelper<odeName>::maxAuxStates_;

  // class type holding fom states data
  using fom_states_data = ::rompp::rom::FomStatesData<
	fom_state_w_t, maxAuxStates, fom_state_reconstr_t>;

  // class type holding fom rhs data
  using fom_rhs_data = ::rompp::rom::FomRhsData<fom_rhs_w_t>;

  // type for user-defined time discrete ops
  using td_ud_ops = time_discrete_ud_ops;
};

}}//end  namespace rompp::rom
#endif
