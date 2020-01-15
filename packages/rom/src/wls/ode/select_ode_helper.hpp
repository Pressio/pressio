#include "./impl/implicit_euler.hpp"
#include "./impl/bdf2.hpp"
#include "./impl/explicit_euler.hpp"

namespace pressio{ namespace rom{ namespace wls{ namespace ode{ namespace helpers{
namespace impl{
template <typename ode_tag, typename fom_state_t, typename wls_state_t>
struct DefaultHelper{
  using ode_type = void;
};

//Explicit Euler specialization
template <typename fom_state_t, typename wls_state_t>
struct DefaultHelper<::pressio::ode::implicitmethods::Euler, fom_state_t, wls_state_t>
{
  using ode_policies_t = ::pressio::rom::wls::ode::impl::ImplicitEuler<fom_state_t, wls_state_t>;
};

//BDF2 specialization
template <typename fom_state_t, typename wls_state_t>
struct DefaultHelper<::pressio::ode::implicitmethods::BDF2, fom_state_t, wls_state_t>
{
  using ode_policies_t = ::pressio::rom::wls::ode::impl::BDF2<fom_state_t, wls_state_t>;
};

//Explicit Euler specialization
template <typename fom_state_t, typename wls_state_t>
struct DefaultHelper<::pressio::ode::explicitmethods::Euler, fom_state_t, wls_state_t>
{
  using ode_policies_t = ::pressio::rom::wls::ode::impl::ExplicitEuler<fom_state_t, wls_state_t>;
};

}

template <typename ode_tag, typename fom_state_t, typename wls_state_t>
using ode_policies_t = typename impl::DefaultHelper<ode_tag,fom_state_t, wls_state_t>::ode_policies_t;

}}}}}

