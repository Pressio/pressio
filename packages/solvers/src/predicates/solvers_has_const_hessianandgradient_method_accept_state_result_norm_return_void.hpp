
#ifndef solvers_has_const_hessianandgradient_method_accept_state_result_return_void_hpp_
#define solvers_has_const_hessianandgradient_method_accept_state_result_return_void_hpp_

namespace pressio{ namespace solvers{ namespace predicates {
  
template <
  typename T,
  typename state_t,
  typename hess_t,
  typename grad_t,
  typename norm_t,
  typename = void
  >
struct has_const_hessianandgradient_method_accept_state_result_norm_return_void
  : std::false_type{};

template <
  typename T,
  typename state_t,
  typename hess_t,
  typename grad_t,
  typename norm_t
  >
struct has_const_hessianandgradient_method_accept_state_result_norm_return_void<
  T, state_t, hess_t, grad_t, norm_t,
  mpl::enable_if_t<
    std::is_void<
      decltype(
         std::declval<T const>().hessianAndGradient
          (
            std::declval<state_t const &>(),
            std::declval<hess_t &>(),
            std::declval<grad_t &>(),
            /* does not matter here what we pass, just to test */
            ::pressio::Norm::Undefined,
            std::declval<norm_t &>()
          )
         )
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::solvers::predicates
#endif
