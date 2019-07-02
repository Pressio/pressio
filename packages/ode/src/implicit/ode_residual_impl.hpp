
#ifndef ODE_RESIDUAL_IMPL_HPP_
#define ODE_RESIDUAL_IMPL_HPP_

#include "../ode_ConfigDefs.hpp"
#include "ode_implicit_constants.hpp"

namespace pressio{ namespace ode{ namespace impl{

// in all functions below, on input R should contain
// the target application RHS, i.e. if the system is
// expressed as dudt = f(x,u,...), then on input R contains f(...)
// On output, R contains the time-discrete residual

template < typename ode::ImplicitEnum method,
	   int n,
	   typename state_type,
	   typename residual_type,
	   typename scalar_type,
	   mpl::enable_if_t<
	     method == ode::ImplicitEnum::Euler
	     > * = nullptr
	   >
void time_discrete_residual(const state_type & y,
			    residual_type & R,
			    const std::array<state_type, n> & oldYs,
			    scalar_type dt) {
  constexpr auto one = ::pressio::utils::constants::one<scalar_type>();
  constexpr auto negOne = ::pressio::utils::constants::negOne<scalar_type>();
  const scalar_type negDt = -dt;
  ::pressio::containers::ops::do_update(R, negDt, y, one, oldYs[0], negOne);
}


template < typename ode::ImplicitEnum method,
	   int n,
	   typename state_type,
	   typename residual_type,
	   typename scalar_type,
	   mpl::enable_if_t<
	     method == ode::ImplicitEnum::BDF2
	     > * = nullptr
	   >
void time_discrete_residual(const state_type & y,
			    residual_type & R,
			    const std::array<state_type, n> & oldYs,
			    scalar_type dt) {

  constexpr auto one = ::pressio::utils::constants::one<scalar_type>();
  constexpr auto negOne = ::pressio::utils::constants::negOne<scalar_type>();

  constexpr auto a = ::pressio::ode::coeffs::bdf2<scalar_type>::c1_*negOne;
  constexpr auto b = ::pressio::ode::coeffs::bdf2<scalar_type>::c2_;
  const auto c = ::pressio::ode::coeffs::bdf2<scalar_type>::c3_*dt*negOne;

  ::pressio::containers::ops::do_update(R, c, y, one, oldYs[1], a, oldYs[0], b);
}

}}}//end namespace pressio::ode::impl
#endif
