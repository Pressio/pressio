
#ifndef SOLVERS_SYSTEM_HAS_ALL_NEEDED_RESIDUAL_METHODS_HPP_
#define SOLVERS_SYSTEM_HAS_ALL_NEEDED_RESIDUAL_METHODS_HPP_

#include "../solvers_ConfigDefs.hpp"
#include "../../../mpl/src/detection_idiom.hpp"

namespace rompp{ namespace solvers{ namespace meta {

template <typename T, typename Arg, typename enable = void>
struct has_residual_method_callable_with_one_arg{
  using type = void;
  static constexpr bool value = false;
};

template <typename T, typename Arg>
struct has_residual_method_callable_with_one_arg<
  T, Arg,
  ::rompp::mpl::void_t<
    decltype(std::declval<T>().residual(std::declval<Arg const&>()))
    >
  >{
  using type = decltype(std::declval<T>().residual(std::declval<Arg const&>()));
  static constexpr bool value = true;
};


template <typename T, typename FirstArg, typename SecondArg, typename = void>
struct has_residual_method_callable_with_two_args : std::false_type{};

template <typename T, typename FirstArg, typename SecondArg>
struct has_residual_method_callable_with_two_args<
  T, FirstArg, SecondArg,
  ::rompp::mpl::void_t<
    decltype(std::declval<T>().residual(std::declval<FirstArg const&>(),
					std::declval<SecondArg&>()))
    >
  > : std::true_type{};



template<
  typename system_type,
  typename state_type,
  typename residual_type,
  typename enable = void
  >
struct system_has_needed_residual_methods : std::false_type{};

template<
  typename system_type,
  typename state_type,
  typename residual_type
  >
struct system_has_needed_residual_methods
< system_type, state_type, residual_type,
  ::rompp::mpl::enable_if_t<
      has_residual_method_callable_with_one_arg<
	system_type, state_type
      >::value and
      has_residual_method_callable_with_two_args<
	system_type, state_type, residual_type
      >::value and
    // residual method with 1 argument returns a residual_type
    std::is_same<
	typename has_residual_method_callable_with_one_arg<
	  system_type, state_type>::type,
      residual_type
      >::value
  >
  > : std::true_type{};

}}} // namespace rompp::solvers::meta
#endif
