
#ifndef ODE_ADMISSIBLE_USER_DEFINED_OPS_EXPLICIT_EULER_HPP_
#define ODE_ADMISSIBLE_USER_DEFINED_OPS_EXPLICIT_EULER_HPP_

namespace pressio{ namespace ode{ namespace concepts {

template<typename T,
	 typename scalar_t,
	 typename state_t,
	 typename residual_t,
	 typename enable = void>
struct user_defined_ops_for_explicit_euler
  : std::false_type{};


template<typename T,
	 typename scalar_t,
	 typename state_t,
	 typename residual_t>
struct user_defined_ops_for_explicit_euler<
  T, scalar_t, state_t, residual_t,
    mpl::enable_if_t<
      ::pressio::ops::meta::has_method_do_update_one_term<
	T,
	scalar_t,
	typename containers::details::traits<state_t>::wrapped_t,
	typename containers::details::traits<residual_t>::wrapped_t
	>::value
      >
  > : std::true_type{};

}}} // namespace pressio::ode::concepts
#endif
