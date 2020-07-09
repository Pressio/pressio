
#ifndef ODE_HAS_CONST_CREATE_TD_RESIDUAL_OBJECT_METHOD_RETURN_RESULT_HPP_
#define ODE_HAS_CONST_CREATE_TD_RESIDUAL_OBJECT_METHOD_RETURN_RESULT_HPP_

namespace pressio{ namespace ode{ namespace predicates {

template <typename T, typename result_t, typename = void>
struct has_const_create_discrete_time_residual_method_return_result
  : std::false_type{};


template <typename T, typename result_t>
struct has_const_create_discrete_time_residual_method_return_result<
  T, result_t,
  ::pressio::mpl::enable_if_t<
    !std::is_void<result_t>::value and
    mpl::is_same<
      result_t,
      decltype(
	       std::declval<T const>().createDiscreteTimeResidual()
	       )
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::ode::predicates
#endif
