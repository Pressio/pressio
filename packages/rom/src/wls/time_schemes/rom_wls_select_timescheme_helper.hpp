
#ifndef ROM_WLS_SELECT_ODE_HELPER_HPP_
#define ROM_WLS_SELECT_ODE_HELPER_HPP_

#include "rom_wls_implicit_euler.hpp"
#include "rom_wls_bdf2.hpp"
#include "rom_wls_explicit_euler.hpp"

namespace pressio{ namespace rom{ namespace wls{ namespace timeschemes{

namespace impl{

template <typename ode_tag, typename fom_state_t, typename wls_state_t>
struct DefaultHelper{
  using type = void;
};

//Explicit Euler specialization
template <typename fom_state_t, typename wls_state_t>
struct DefaultHelper<::pressio::ode::implicitmethods::Euler, fom_state_t, wls_state_t>{
  using type = ::pressio::rom::wls::timeschemes::impl::ImplicitEuler<fom_state_t, wls_state_t>;
};

//BDF2 specialization
template <typename fom_state_t, typename wls_state_t>
struct DefaultHelper<::pressio::ode::implicitmethods::BDF2, fom_state_t, wls_state_t>{
  using type = ::pressio::rom::wls::timeschemes::impl::BDF2<fom_state_t, wls_state_t>;
};

//Explicit Euler specialization
template <typename fom_state_t, typename wls_state_t>
struct DefaultHelper<::pressio::ode::explicitmethods::Euler, fom_state_t, wls_state_t>{
  using type = ::pressio::rom::wls::timeschemes::impl::ExplicitEuler<fom_state_t, wls_state_t>;
};

} //end namespace pressio::rom::wls::timeschemes::helpers::impl

template <typename ode_tag, typename fom_state_t, typename wls_state_t>
using timescheme_t = typename impl::DefaultHelper<ode_tag,fom_state_t, wls_state_t>::type;

}}}} //ena namespace pressio::rom::wls::timeschemes::helpers
#endif
