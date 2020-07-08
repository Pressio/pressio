
#ifndef ode_has_const_create_jacobian_method_return_result_HPP_
#define ode_has_const_create_jacobian_method_return_result_HPP_

namespace pressio{ namespace ode{ namespace predicates {
  
template <
  typename T,
  typename jacobian_t,
  typename = void
  >
struct has_const_create_jacobian_method_return_result
  : std::false_type{};

template <
  typename T,
  typename jacobian_t
  >
struct has_const_create_jacobian_method_return_result<
  T, jacobian_t,
  mpl::enable_if_t<
    !std::is_void<jacobian_t>::value and
    std::is_same<
      jacobian_t,
      decltype(
         std::declval<T const>().createJacobian()
         )
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::ode::predicates
#endif
