
#ifndef ODE_admissible_USER_DEFINED_OPS_EXPLICIT_RK4_HPP_
#define ODE_admissible_USER_DEFINED_OPS_EXPLICIT_RK4_HPP_

namespace pressio{ namespace ode{ namespace meta {

template<typename T,
	 typename scalar_t,
	 typename state_t,
	 typename residual_t,
	 typename enable = void>
struct admissible_user_defined_ops_for_explicit_rk4
  : std::false_type{};

template<typename T,
	 typename scalar_t,
	 typename state_t,
	 typename residual_t>
struct admissible_user_defined_ops_for_explicit_rk4<
  T, scalar_t, state_t, residual_t,
    mpl::enable_if_t<
      ::pressio::ops::meta::has_method_do_update_two_terms<
	T,
	scalar_t,
	typename containers::details::traits<state_t>::wrapped_t,
	typename containers::details::traits<residual_t>::wrapped_t,
	typename containers::details::traits<residual_t>::wrapped_t
	>::value and
      ::pressio::ops::meta::has_method_do_update_four_terms<
	T,
	scalar_t,
	typename containers::details::traits<state_t>::wrapped_t,
	typename containers::details::traits<residual_t>::wrapped_t,
	typename containers::details::traits<residual_t>::wrapped_t,
	typename containers::details::traits<residual_t>::wrapped_t,
	typename containers::details::traits<residual_t>::wrapped_t
	>::value
      >
  > : std::true_type{};

}}} // namespace pressio::ode::meta
#endif
