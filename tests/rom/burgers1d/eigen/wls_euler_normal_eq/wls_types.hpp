
#ifndef ROM_LSPG_STEADY_TYPE_GENERATOR_COMMON_HPP_
#define ROM_LSPG_STEADY_TYPE_GENERATOR_COMMON_HPP_

#include "../../../../../packages/rom/src/rom_fwd.hpp"
#include "../../../../../packages/rom/src/rom_static_container_fom_states.hpp"
#include "../../../../../packages/rom/src/rom_reconstructor_fom_state.hpp"

namespace pressio{ namespace rom{ namespace wls{ 

template <
  typename fom_type,
  typename decoder_type,
  typename wls_state_type
  >
struct CommonTypes<
  fom_type, decoder_type, wls_state_type,
  mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper<wls_state_type>::value
    >
  >
{
  // these are native types of the full-order model (fom)
  using fom_t                   = fom_type;
  using scalar_t                = typename fom_t::scalar_type;
  using fom_native_state_t      = typename fom_t::state_type;
//  using fom_native_velocity_t   = typename fom_t::velocity_type;

  // fom wrapper types
  using fom_state_t     = ::pressio::containers::Vector<fom_native_state_t>;
//  using fom_velocity_t  = ::pressio::containers::Vector<fom_native_velocity_t>;

  // rom state type (passed in)
  using wls_state_t            = wls_state_type;

  // for LSPG, the rom residual type = containers::wrapper of application rhs
  // i.e. the wrapped fom rhs type
//  using lspg_residual_t         = fom_velocity_t;

  // decoder types (passed in)
  using decoder_t               = decoder_type;
  using decoder_jac_t           = typename decoder_t::jacobian_t;

  // fom state reconstructor type
  using fom_state_reconstr_t    = FomStateReconstructor<fom_state_t, decoder_t>;

  // class type holding fom states data: we only need to store one FOM state
  using fom_states_data = ::pressio::rom::FomStatesStaticContainer<fom_state_t, 1, fom_state_reconstr_t>;
};

}}}//end  namespace pressio::rom::lspg::steady
#endif


