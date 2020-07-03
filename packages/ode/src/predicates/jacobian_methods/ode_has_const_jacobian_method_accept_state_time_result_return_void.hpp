
#ifndef ODE_HAS_CONST_JACOBIAN_METHOD_ACCEPT_STATE_TIME_RESULT_RETURN_VOID_HPP_
#define ODE_HAS_CONST_JACOBIAN_METHOD_ACCEPT_STATE_TIME_RESULT_RETURN_VOID_HPP_

namespace pressio{ namespace ode{ namespace meta {

template <
  typename T,
  typename state_t,
  typename time_type,
  typename jac_t,
  typename = void
  >
struct has_const_jacobian_method_accept_state_time_result_return_void
  : std::false_type{};

template <
  typename T,
  typename state_t,
  typename time_type,
  typename jac_t
  >
struct has_const_jacobian_method_accept_state_time_result_return_void<
  T, state_t, time_type, jac_t,
  ::pressio::mpl::void_t<
  decltype(
	   std::declval<T const>().jacobian(
					    std::declval<state_t const&>(),
					    std::declval<time_type const &>(),
					    std::declval<jac_t &>()
					    )
	   )
    >
  >: std::true_type{};

}}} // namespace pressio::ode::meta
#endif
