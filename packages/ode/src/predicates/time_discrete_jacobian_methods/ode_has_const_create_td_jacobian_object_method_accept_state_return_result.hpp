
#ifndef ODE_HAS_CONST_CREATE_TD_JAC_OBJECT_METHOD_ACCEPT_STATE_RETURN_RESULT_HPP_
#define ODE_HAS_CONST_CREATE_TD_JAC_OBJECT_METHOD_ACCEPT_STATE_RETURN_RESULT_HPP_

namespace pressio{ namespace ode{ namespace meta {

template <
  typename T, typename state_t, typename jacobian_t,
  typename = void
  >
struct has_const_create_td_jacobian_object_method_accept_state_return_result
  : std::false_type{};


template <typename T, typename state_t, typename jacobian_t>
struct has_const_create_td_jacobian_object_method_accept_state_return_result<
  T, state_t, jacobian_t,
  ::pressio::mpl::enable_if_t<
    !std::is_void<jacobian_t>::value and
    mpl::is_same<
      jacobian_t,
      decltype(
	       std::declval<T const>().createTimeDiscreteJacobianObject(
              std::declval<state_t const&>()
            )
	       )
      >::value
    >
  > : std::true_type{};


}}} // namespace pressio::ode::meta
#endif
