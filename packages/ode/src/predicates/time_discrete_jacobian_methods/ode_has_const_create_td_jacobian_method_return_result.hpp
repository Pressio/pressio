
#ifndef ODE_HAS_CONST_CREATE_TD_JAC_OBJECT_METHOD_RETURN_RESULT_HPP_
#define ODE_HAS_CONST_CREATE_TD_JAC_OBJECT_METHOD_RETURN_RESULT_HPP_

namespace pressio{ namespace ode{ namespace meta {

template <
  typename T, typename jacobian_t,
  typename = void
  >
struct has_const_create_td_jacobian_method_return_result
  : std::false_type{};


template <typename T, typename jacobian_t>
struct has_const_create_td_jacobian_method_return_result<
  T, jacobian_t,
  ::pressio::mpl::enable_if_t<
    !std::is_void<jacobian_t>::value and
    mpl::is_same<
      jacobian_t,
      decltype(
	       std::declval<T const>().createTimeDiscreteJacobian()
	       )
      >::value
    >
  > : std::true_type{};


}}} // namespace pressio::ode::meta
#endif
