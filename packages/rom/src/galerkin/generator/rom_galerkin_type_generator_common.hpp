
#ifndef ROM_GALERKIN_TYPE_GENERATOR_COMMON_HPP_
#define ROM_GALERKIN_TYPE_GENERATOR_COMMON_HPP_

#include "../../rom_ConfigDefs.hpp"
#include "../../rom_fwd.hpp"
#include "../../rom_data_fom_rhs.hpp"
#include "../../rom_data_fom_states.hpp"
#include "../../policies/rom_evaluate_fom_velocity_unsteady_policy.hpp"
#include "../../policies/rom_apply_fom_jacobian_unsteady_policy.hpp"
#include "../../../../ode/src/ode_fwd.hpp"
#include "../../../../ode/src/meta/ode_is_legitimate_model_for_explicit_ode.hpp"
#include "../../meta/rom_is_legitimate_decoder_type.hpp"

namespace pressio{ namespace rom{

template < typename galerkin_state_type, typename ...Args >
struct GalerkinCommonTypes{

  // verify that args contains a valid fom/adapter type
  // a valid fom/adapter type is one that is valid for explicit ode
  using ic1 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::ode::meta::is_legitimate_model_for_explicit_ode, Args...>;
  using fom_t = ::pressio::mpl::variadic::at_or_t<void, ic1::value, Args...>;
  static_assert(!std::is_void<fom_t>::value and ic1::value < sizeof... (Args),
		"A valid adapter/fom type must be passed to define a ROM Galerkin problem");

  // get the native types from the full-order model (fom)
  using scalar_t		= typename fom_t::scalar_type;
  using fom_native_state_t	= typename fom_t::state_type;
  using fom_native_velocity_t	= typename fom_t::velocity_type;

  // declare fom wrapper types
  using fom_state_t		= ::pressio::containers::Vector<fom_native_state_t>;
  using fom_velocity_t		= ::pressio::containers::Vector<fom_native_velocity_t>;

  // rom state type (passed in)
  using galerkin_state_t	= galerkin_state_type;

  // the GALERKIN residual type is (for now) same as state type
  using galerkin_residual_t	= galerkin_state_type;

  // verify the sequence contains a valid decoder type
  using ic2 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::rom::meta::is_legitimate_decoder_type, Args...>;
  using decoder_t = ::pressio::mpl::variadic::at_or_t<void, ic2::value, Args...>;
  static_assert(!std::is_void<decoder_t>::value and ic2::value < sizeof... (Args),
		"A valid decoder type must be passed to define a ROM Galerkin problem");
  using decoder_jac_t		= typename decoder_t::jacobian_t;

  // fom state reconstructor type
  using fom_state_reconstr_t	= FomStateReconstructor<fom_state_t, decoder_t>;

  // class type holding fom states data
  using fom_states_data = ::pressio::rom::FomStatesData<
	fom_state_t, 0, fom_state_reconstr_t>;

  // class type holding fom rhs data
  using fom_velocity_data = ::pressio::rom::FomRhsData<fom_velocity_t>;
};

}}//end  namespace pressio::rom
#endif
