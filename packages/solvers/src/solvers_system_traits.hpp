
#ifndef SOLVERS_SYSTEM_TRAITS_HPP
#define SOLVERS_SYSTEM_TRAITS_HPP

#include "../../core/src/meta/core_meta_detection_idiom.hpp"


namespace rompp{ namespace solvers{ namespace details {

template <typename T>
using has_public_state_type = typename T::state_type;

template <typename T>
using has_public_residual_type = typename T::residual_type;

template <typename T>
using has_public_jacobian_type = typename T::jacobian_type;


template <typename T, typename Arg>
using has_residual_callable_with_one_arg =
  decltype(std::declval<T>().residual(std::declval<Arg const&>()));

template <typename T, typename FirstArg, typename SecondArg>
using has_residual_callable_with_two_args =
  decltype(std::declval<T>().residual(std::declval<FirstArg const&>(),
				      std::declval<SecondArg&>()));

template <typename T, typename Arg>
using has_jacobian_callable_with_one_arg =
  decltype(std::declval<T>().jacobian(std::declval<Arg const&>()));

template <typename T, typename FirstArg, typename SecondArg>
using has_jacobian_callable_with_two_args =
  decltype(std::declval<T>().jacobian(std::declval<FirstArg const&>(),
				      std::declval<SecondArg&>()));

template <typename T>
struct system_traits {
  using state_type = typename core::meta::detected_t<has_public_state_type, T>;
  using residual_type = typename core::meta::detected_t<has_public_residual_type, T>;
  using jacobian_type = typename core::meta::detected_t<has_public_jacobian_type, T>;

  static constexpr bool has_public_state_type =
    core::meta::is_detected<solvers::details::has_public_state_type, T>::value;

  static constexpr bool has_public_residual_type =
    core::meta::is_detected<solvers::details::has_public_residual_type, T>::value;

  static constexpr bool has_public_jacobian_type =
    core::meta::is_detected<solvers::details::has_public_jacobian_type, T>::value;


  static constexpr bool has_residual_callable_with_one_arg =
    core::meta::is_detected<solvers::details::has_residual_callable_with_one_arg,
			    T, state_type>::value;
  static constexpr bool has_residual_callable_with_two_args =
    core::meta::is_detected<solvers::details::has_residual_callable_with_two_args,
			    T, state_type, residual_type>::value;
  static constexpr bool has_residual_methods =
    has_residual_callable_with_one_arg || has_residual_callable_with_two_args;


  static constexpr bool has_jacobian_callable_with_one_arg =
    core::meta::is_detected<solvers::details::has_jacobian_callable_with_one_arg,
			    T, state_type>::value;
  static constexpr bool has_jacobian_callable_with_two_args =
    core::meta::is_detected<solvers::details::has_jacobian_callable_with_two_args,
			    T, state_type, jacobian_type>::value;
  static constexpr bool has_jacobian_methods =
    has_jacobian_callable_with_one_arg || has_jacobian_callable_with_two_args;

  static constexpr bool is_system = has_residual_methods and  has_jacobian_methods;
};


}}}//end namespace rompp::solvers::details
#endif
