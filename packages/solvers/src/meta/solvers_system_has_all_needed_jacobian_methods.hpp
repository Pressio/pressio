
#ifndef SOLVERS_SYSTEM_HAS_ALL_NEEDED_JACOBIAN_METHODS_HPP_
#define SOLVERS_SYSTEM_HAS_ALL_NEEDED_JACOBIAN_METHODS_HPP_

#include "../solvers_ConfigDefs.hpp"
#include "../../../mpl/src/detection_idiom.hpp"

namespace rompp{ namespace solvers{ namespace meta {

template <typename T, typename Arg>
using has_jacobian_method_callable_with_one_arg =
  decltype(std::declval<T>().jacobian(std::declval<Arg const&>()));

template <typename T, typename FirstArg, typename SecondArg>
using has_jacobian_method_callable_with_two_args =
  decltype(std::declval<T>().jacobian(std::declval<FirstArg const&>(),
				      std::declval<SecondArg&>()));

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
  ::rompp::mpl::enable_if_t<
   // has method with 1 arguments,
    ::rompp::mpl::is_detected<
      has_jacobian_method_callable_with_one_arg,
      system_type, state_type
      >::value and
   // has method with 2 arguments,
    ::rompp::mpl::is_detected<
      has_jacobian_method_callable_with_two_args,
      system_type, state_type, jacobian_type
      >::value and
    // method with 1 argument returns a jacobian_type
    std::is_same<
      ::rompp::mpl::detected_t<
	has_jacobian_method_callable_with_one_arg,
	system_type, state_type>,
      jacobian_type
      >::value
  >
  > : std::true_type{};

}}} // namespace rompp::solvers::meta
#endif
