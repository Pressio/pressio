
#ifndef ODE_IS_VALID_USER_DEFINED_OPS_EXPLICIT_ODE_HPP_
#define ODE_IS_VALID_USER_DEFINED_OPS_EXPLICIT_ODE_HPP_

#include "ode_is_valid_user_defined_ops_for_explicit_euler.hpp"
#include "ode_is_valid_user_defined_ops_for_explicit_rk4.hpp"

namespace pressio{ namespace ode{ namespace meta {

template<typename T,
	 typename scalar_t,
	 typename state_t,
	 typename residual_t,
	 typename enable = void>
struct is_valid_user_defined_ops_for_explicit_ode : std::false_type{};

template<typename T,
	 typename scalar_t,
	 typename state_t,
	 typename residual_t>
struct is_valid_user_defined_ops_for_explicit_ode<
  T, scalar_t, state_t, residual_t,
    mpl::enable_if_t<
      is_valid_user_defined_ops_for_explicit_euler<
	T, scalar_t, state_t, residual_t
      	>::value and
      is_valid_user_defined_ops_for_explicit_rk4<
      	T, scalar_t, state_t, residual_t
      	>::value
      >
  > : std::true_type{};

}}} // namespace pressio::ode::meta
#endif
