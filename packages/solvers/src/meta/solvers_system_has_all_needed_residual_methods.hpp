
#ifndef SOLVERS_SYSTEM_HAS_ALL_NEEDED_RESIDUAL_METHODS_HPP_
#define SOLVERS_SYSTEM_HAS_ALL_NEEDED_RESIDUAL_METHODS_HPP_

#include "../solvers_ConfigDefs.hpp"
#include "../../../mpl/src/detection_idiom.hpp"

namespace rompp{ namespace solvers{ namespace meta {

template <typename T, typename Arg>
using has_residual_method_callable_with_one_arg =
  decltype(std::declval<T>().residual(std::declval<Arg const&>()));

template <typename T, typename FirstArg, typename SecondArg>
using has_residual_method_callable_with_two_args =
  decltype(std::declval<T>().residual(std::declval<FirstArg const&>(),
				      std::declval<SecondArg&>()));

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
    // has residual method with 1 arguments,
    ::rompp::mpl::is_detected<
      has_residual_method_callable_with_one_arg,
      system_type, state_type
      >::value and
    // has residual method with 2 arguments,
    ::rompp::mpl::is_detected<
      has_residual_method_callable_with_two_args,
      system_type, state_type, residual_type
      >::value and
    // residual method with 1 argument returns a residual_type
    std::is_same<
      ::rompp::mpl::detected_t<
	has_residual_method_callable_with_one_arg,
	system_type, state_type>,
      residual_type
      >::value
  >
  > : std::true_type{};

}}} // namespace rompp::solvers::meta
#endif
