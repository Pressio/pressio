
#ifndef SOLVERS_EXPERIMENTAL_SYSTEM_TRAITS_HPP
#define SOLVERS_EXPERIMENTAL_SYSTEM_TRAITS_HPP

#include <type_traits>

#include "../../../algebra/src/meta/algebra_meta_detection_idiom.hpp"


namespace rompp{
namespace solvers{
namespace details {


template <typename T>
using has_public_matrix_type = typename T::matrix_type;


template <typename T>
using has_public_vector_type = typename T::vector_type;


template <
  typename T,
  typename Arg
>
using has_residual_callable_with_one_arg =
  decltype(std::declval<T>().residual(std::declval<Arg const&>()));


template <
  typename T,
  typename FirstArg,
  typename SecondArg
>
using has_residual_callable_with_two_args =
  decltype(std::declval<T>().residual(std::declval<FirstArg const&>(),
				      std::declval<SecondArg&>()));


template <
  typename T,
  typename Arg
>
using has_jacobian_callable_with_one_arg =
  decltype(std::declval<T>().jacobian(std::declval<Arg const&>()));


template <
  typename T,
  typename FirstArg,
  typename SecondArg
>
using has_jacobian_callable_with_two_args =
  decltype(std::declval<T>().jacobian(std::declval<FirstArg const&>(),
				      std::declval<SecondArg&>()));


template <typename T>
struct system_traits {
  typedef typename algebra::meta::detected_t<has_public_vector_type, T> vector_type;
  typedef typename algebra::meta::detected_t<has_public_matrix_type, T> matrix_type;

  static constexpr bool has_public_vector_type =
    algebra::meta::is_detected<has_public_vector_type, T>::value;
  static constexpr bool has_public_matrix_type =
    algebra::meta::is_detected<has_public_matrix_type, T>::value;

  static constexpr bool has_residual_callable_with_one_arg =
    algebra::meta::is_detected<has_residual_callable_with_one_arg,
			    T, vector_type>::value;
  static constexpr bool has_residual_callable_with_two_args =
    algebra::meta::is_detected<has_residual_callable_with_two_args,
			    T, vector_type, vector_type>::value;
  static constexpr bool has_residual_methods =
    has_residual_callable_with_one_arg || has_residual_callable_with_two_args;

  static constexpr bool has_jacobian_callable_with_one_arg =
    algebra::meta::is_detected<has_jacobian_callable_with_one_arg,
			    T, vector_type>::value;
  static constexpr bool has_jacobian_callable_with_two_args =
    algebra::meta::is_detected<has_jacobian_callable_with_two_args,
			    T, vector_type, matrix_type>::value;
  static constexpr bool has_jacobian_methods =
    has_jacobian_callable_with_one_arg || has_jacobian_callable_with_two_args;

  static constexpr bool is_system = has_residual_methods and  has_jacobian_methods;
};

} // end namespace details
} // end namespace solvers

}//end namespace rompp
#endif
