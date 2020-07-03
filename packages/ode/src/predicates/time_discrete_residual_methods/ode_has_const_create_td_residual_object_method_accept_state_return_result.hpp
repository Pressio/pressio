
#ifndef ODE_HAS_CONST_CREATE_TD_RESIDUAL_OBJECT_METHOD_ACCEPT_STATE_RETURN_RESULT_HPP_
#define ODE_HAS_CONST_CREATE_TD_RESIDUAL_OBJECT_METHOD_ACCEPT_STATE_RETURN_RESULT_HPP_

namespace pressio{ namespace ode{ namespace meta {

template <
  typename T, typename state_t, typename residual_t,
  typename = void
  >
struct has_const_create_td_residual_object_method_accept_state_return_result
  : std::false_type{};


template <typename T, typename state_t, typename residual_t>
struct has_const_create_td_residual_object_method_accept_state_return_result<
  T, state_t, residual_t,
  ::pressio::mpl::enable_if_t<
    !std::is_void<residual_t>::value and
    mpl::is_same<
      residual_t,
      decltype(
	       std::declval<T const>().createTimeDiscreteResidualObject(
              std::declval<state_t const&>()
              )
	       )
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::ode::meta
#endif
