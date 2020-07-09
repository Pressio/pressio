
#ifndef has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void_hpp_
#define has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void_hpp_

namespace pressio{ namespace rom{ namespace predicates {

template <
  typename T, int n, typename step_t, typename sc_t, typename state_t, typename operand_t, typename result_t,
  typename = void
  >
struct has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void
  : std::false_type{};


template <typename T, typename step_t, typename sc_t, typename state_t, typename operand_t, typename result_t>
struct has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void<
  T, 1, step_t, sc_t, state_t, operand_t, result_t,
  ::pressio::mpl::enable_if_t<
    std::is_void<
      decltype(
	       std::declval<T const>().applyDiscreteTimeJacobian(
							    std::declval<step_t const &>(),
							    std::declval<sc_t const &>(),
							    std::declval<sc_t const &>(),
							    std::declval<operand_t const &>(),
							    std::declval<result_t &>(),
							    std::declval<state_t const&>()
							    )
	       )
      >::value
    >
  > : std::true_type{};


template <typename T, typename step_t, typename sc_t, typename state_t, typename operand_t, typename result_t>
struct has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void<
  T, 2, step_t, sc_t, state_t, operand_t, result_t,
  ::pressio::mpl::enable_if_t<
    std::is_void<
      decltype(
	       std::declval<T const>().applyDiscreteTimeJacobian(
							    std::declval<step_t const &>(),
							    std::declval<sc_t const &>(),
							    std::declval<sc_t const &>(),
							    std::declval<operand_t const &>(),
							    std::declval<result_t &>(),
							    std::declval<state_t const&>(),
							    std::declval<state_t const&>()
							    )
	       )
      >::value
    >
  > : std::true_type{};


template <typename T, typename step_t, typename sc_t, typename state_t, typename operand_t, typename result_t>
struct has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void<
  T, 3, step_t, sc_t, state_t, operand_t, result_t,
  ::pressio::mpl::enable_if_t<
    std::is_void<
      decltype(
	       std::declval<T const>().applyDiscreteTimeJacobian(
							    std::declval<step_t const &>(),
							    std::declval<sc_t const &>(),
							    std::declval<sc_t const &>(),
							    std::declval<operand_t const &>(),
							    std::declval<result_t &>(),
							    std::declval<state_t const&>(),
							    std::declval<state_t const&>(),
							    std::declval<state_t const&>()
							    )
	       )
      >::value
    >
  > : std::true_type{};

}}} 
#endif
