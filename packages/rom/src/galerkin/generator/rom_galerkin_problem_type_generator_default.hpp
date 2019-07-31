
#ifndef ROM_GALERKIN_PROBLEM_TYPE_GENERATOR_DEFAULT_HPP_
#define ROM_GALERKIN_PROBLEM_TYPE_GENERATOR_DEFAULT_HPP_

#include "rom_galerkin_type_generator_common.hpp"

namespace pressio{ namespace rom{

template <
  ode::ExplicitEnum odeName,
  typename galerkin_state_type,
  typename ...Args
  >
struct DefaultGalerkinExplicitTypeGenerator
  : GalerkinCommonTypes<galerkin_state_type, Args...>
{

  using base_t = GalerkinCommonTypes<galerkin_state_type, Args...>;

  static constexpr ode::ExplicitEnum odeName_ = odeName;

  using typename base_t::fom_t;
  using typename base_t::scalar_t;
  using typename base_t::fom_native_state_t;
  using typename base_t::fom_state_t;
  using typename base_t::fom_velocity_t;
  using typename base_t::galerkin_state_t;
  using typename base_t::galerkin_residual_t;
  using typename base_t::decoder_t;
  using typename base_t::decoder_jac_t;
  using typename base_t::fom_state_reconstr_t;
  using typename base_t::fom_states_data;
  using typename base_t::fom_velocity_data;

  // policy for evaluating the ode residual
  using galerkin_residual_policy_t =
    ::pressio::rom::DefaultGalerkinExplicitVelocityPolicy<
    fom_states_data, fom_velocity_data, decoder_t>;

  // declare type of stepper object
  using galerkin_stepper_t = ::pressio::ode::ExplicitStepper<
    odeName, galerkin_state_type, fom_t,
    galerkin_residual_t, galerkin_residual_policy_t, scalar_t>;

};//end class

}}//end  namespace pressio::rom
#endif
