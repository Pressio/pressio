
#ifndef solvers_has_const_jacobian_method_accept_state_result_return_void_HPP_
#define solvers_has_const_jacobian_method_accept_state_result_return_void_HPP_

namespace pressio{ namespace solvers{ namespace meta {
  
template <
  typename T,
  typename state_t,
  typename jac_t,
  typename = void
  >
struct has_const_jacobian_method_accept_state_result_return_void
  : std::false_type{};

template <
  typename T,
  typename state_t,
  typename jac_t
  >
struct has_const_jacobian_method_accept_state_result_return_void<
  T, state_t, jac_t,
  mpl::enable_if_t<
    std::is_void<
      decltype(
         std::declval<T const>().jacobian
            (
              std::declval<state_t const &>(),
              std::declval<jac_t &>()
            )
         )
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::solvers::meta
#endif
