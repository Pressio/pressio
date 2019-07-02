
#ifndef SOLVERS_SYSTEM_HAS_ALL_NEEDED_JACOBIAN_METHODS_HPP_
#define SOLVERS_SYSTEM_HAS_ALL_NEEDED_JACOBIAN_METHODS_HPP_

#include "../solvers_ConfigDefs.hpp"
#include "../../../mpl/src/detection_idiom.hpp"

namespace pressio{ namespace solvers{ namespace meta {

template <typename T, typename Arg, typename enable = void>
struct has_jacobian_method_callable_with_one_arg{
  using type = void;
  static constexpr bool value = false;
};

template <typename T, typename Arg>
struct has_jacobian_method_callable_with_one_arg<
  T, Arg,
  ::pressio::mpl::void_t<
    decltype(std::declval<T>().jacobian(std::declval<Arg const&>()))
    >
  >{
  using type = decltype(std::declval<T>().jacobian(std::declval<Arg const&>()));
  static constexpr bool value = true;
};


template <typename T, typename FirstArg, typename SecondArg, typename = void>
struct has_jacobian_method_callable_with_two_args : std::false_type{};

template <typename T, typename FirstArg, typename SecondArg>
struct has_jacobian_method_callable_with_two_args<
  T, FirstArg, SecondArg,
  ::pressio::mpl::void_t<
    decltype(std::declval<T>().jacobian(std::declval<FirstArg const&>(),
					std::declval<SecondArg&>()))
    >
  > : std::true_type{};



template<
  typename system_type,
  typename state_type,
  typename jacobian_type,
  typename enable = void
  >
struct system_has_needed_jacobian_methods : std::false_type{};

template<
  typename system_type,
  typename state_type,
  typename jacobian_type
  >
struct system_has_needed_jacobian_methods
< system_type, state_type, jacobian_type,
  ::pressio::mpl::enable_if_t<
      has_jacobian_method_callable_with_one_arg<
	system_type, state_type
      >::value and
   // has method with 2 arguments,
      has_jacobian_method_callable_with_two_args<
	system_type, state_type, jacobian_type
      >::value and
    // method with 1 argument returns a jacobian_type
    std::is_same<
	typename has_jacobian_method_callable_with_one_arg<
	  system_type, state_type>::type,
      jacobian_type
      >::value
  >
  > : std::true_type{};

}}} // namespace pressio::solvers::meta
#endif
