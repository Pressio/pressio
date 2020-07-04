
#ifndef solvers_has_const_residualnorm_method_accept_state_norm_return_void_hpp_
#define solvers_has_const_residualnorm_method_accept_state_norm_return_void_hpp_

namespace pressio{ namespace solvers{ namespace meta {
  

template <
  typename T,
  typename state_t,
  typename norm_t,
  typename = void
  >
struct has_const_residualnorm_method_accept_state_norm_return_void
  : std::false_type{};

template <
  typename T,
  typename state_t,
  typename norm_t
  >
struct has_const_residualnorm_method_accept_state_norm_return_void<
  T, state_t, norm_t, 
  mpl::enable_if_t<
    std::is_void<
      decltype(
         std::declval<T const>().residualNorm
            (
              std::declval<state_t const &>(),
              ::pressio::Norm::Undefined,
              std::declval<norm_t &>()
            )
         )
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::solvers::meta
#endif
