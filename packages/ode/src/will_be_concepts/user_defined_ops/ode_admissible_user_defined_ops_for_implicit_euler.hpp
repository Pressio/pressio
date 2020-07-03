
#ifndef ODE_IS_VALID_USER_DEFINED_OPS_IMPLICIT_EULER_HPP_
#define ODE_IS_VALID_USER_DEFINED_OPS_IMPLICIT_EULER_HPP_

namespace pressio{ namespace ode{ namespace meta {

template<
  typename T,
  typename scalar_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t,
  typename enable = void
  >
struct admissible_user_defined_ops_for_implicit_euler : std::false_type{};

template<
  typename T,
  typename scalar_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t
  >
struct admissible_user_defined_ops_for_implicit_euler<
  T, scalar_t, state_t, residual_t, jacobian_t,
    mpl::enable_if_t<
      ::pressio::ops::meta::has_method_do_update_two_terms<
	typename T::update_op,
	scalar_t,
	typename containers::details::traits<residual_t>::wrapped_t,
	typename containers::details::traits<state_t>::wrapped_t,
	typename containers::details::traits<state_t>::wrapped_t
	>::value
      and
      ::pressio::ops::meta::has_method_scale<
	typename T::update_op,
	typename containers::details::traits<jacobian_t>::wrapped_t,
	scalar_t
	>::value
      and
      ::pressio::ops::meta::has_method_add_to_diagonal<
	typename T::update_op,
	typename containers::details::traits<jacobian_t>::wrapped_t,
	scalar_t
	>::value
      >
  > : std::true_type{};

}}} // namespace pressio::ode::meta
#endif
