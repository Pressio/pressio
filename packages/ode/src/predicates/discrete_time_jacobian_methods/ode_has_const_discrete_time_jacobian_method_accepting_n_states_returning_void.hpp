
#ifndef ODE_HAS_CONST_TD_JACOBIAN_METHOD_ACCEPT_STEP_TIME_DT_RESULT_STATES_RETURN_VOID_HPP_
#define ODE_HAS_CONST_TD_JACOBIAN_METHOD_ACCEPT_STEP_TIME_DT_RESULT_STATES_RETURN_VOID_HPP_

namespace pressio{ namespace ode{ namespace predicates {

template <
  typename T, int n, typename step_t, typename sc_t, typename state_t, typename result_t,
  typename = void
  >
struct has_const_discrete_time_jacobian_method_accept_step_time_dt_result_n_states_return_void
  : std::false_type{};


template <typename T, typename step_t, typename sc_t, typename state_t, typename result_t>
struct has_const_discrete_time_jacobian_method_accept_step_time_dt_result_n_states_return_void<
  T, 1, step_t, sc_t, state_t, result_t,
  ::pressio::mpl::enable_if_t<
    std::is_void<
      decltype(
	       std::declval<T const>().discreteTimeJacobian(
							    std::declval<step_t const &>(),
							    std::declval<sc_t const &>(),
							    std::declval<sc_t const &>(),
							    std::declval<result_t &>(),
							    std::declval<state_t const&>()
							    )
	       )
      >::value
    >
  > : std::true_type{};


template <typename T, typename step_t, typename sc_t, typename state_t, typename result_t>
struct has_const_discrete_time_jacobian_method_accept_step_time_dt_result_n_states_return_void<
  T, 2, step_t, sc_t, state_t, result_t,
  ::pressio::mpl::enable_if_t<
    std::is_void<
      decltype(
	       std::declval<T const>().discreteTimeJacobian(
							    std::declval<step_t const &>(),
							    std::declval<sc_t const &>(),
							    std::declval<sc_t const &>(),
							    std::declval<result_t &>(),
							    std::declval<state_t const&>(),
							    std::declval<state_t const&>()
							    )
	       )
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::ode::predicates
#endif
