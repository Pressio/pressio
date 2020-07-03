
#ifndef ODE_HAS_CONST_JACOBIAN_METHOD_ACCEPT_STATE_TIME_RETURN_RESULT_HPP_
#define ODE_HAS_CONST_JACOBIAN_METHOD_ACCEPT_STATE_TIME_RETURN_RESULT_HPP_

namespace pressio{ namespace ode{ namespace meta {
  

template <
  typename T,
  typename state_t,
  typename time_type,
  typename jacobian_t,
  typename = void
  >
struct has_const_jacobian_method_accept_state_time_return_result
  : std::false_type{};

template <
  typename T,
  typename state_t,
  typename time_type,
  typename jacobian_t
  >
struct has_const_jacobian_method_accept_state_time_return_result<
  T, state_t, time_type, jacobian_t,
  mpl::enable_if_t<
    !std::is_void<jacobian_t>::value and
    std::is_same<
      jacobian_t,
      decltype(
         std::declval<T const>().jacobian(
            std::declval<state_t const &>(),
            std::declval<time_type const &>()
            )
         )
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::ode::meta
#endif
