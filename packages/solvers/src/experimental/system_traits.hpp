
#ifndef SYSTEM_TRAITS_HPP_
#define SYSTEM_TRAITS_HPP_

#include <type_traits>

#include <Eigen/Core>

#include "meta/core_meta_detection_idiom.hpp"
#include "matrix/core_matrix_traits_exp.hpp"
#include "vector/core_vector_traits_exp.hpp"


namespace solvers {
namespace details {


template <typename T>
using has_public_matrix_type = typename T::matrix_type;


template <typename T>
using has_public_state_type = typename T::state_type;


template <typename T, typename Arg>
using has_residual_callable_with_one_arg = decltype(std::declval<T>().residual(std::declval<Arg const&>()));


template <typename T, typename Arg> 
using has_residual_callable_with_two_args = decltype(std::declval<T>().residual(std::declval<Arg const&>(), std::declval<Arg&>()));


template <typename T, typename Arg>
using has_jacobian_callable_with_one_arg = decltype(std::declval<T>().jacobian(std::declval<Arg const&>()));


template <typename T, typename FirstArg, typename SecondArg>
using has_jacobian_callable_with_two_args = decltype(std::declval<T>().jacobian(std::declval<FirstArg const&>(), std::declval<SecondArg&>()));


template <typename T>
struct system_traits {
  typedef typename core::meta::detected_t<has_public_state_type, T> state_type;
  typedef typename core::meta::detected_t<has_public_matrix_type, T> matrix_type;

  static constexpr bool has_public_state_type = core::meta::is_detected<has_public_state_type, T>::value;
  static constexpr bool has_public_matrix_type = core::meta::is_detected<has_public_matrix_type, T>::value;

  static constexpr bool has_residual_methods = 
    core::meta::is_detected<has_residual_callable_with_one_arg, T, state_type>::value &&
    core::meta::is_detected<has_residual_callable_with_two_args, T, state_type>::value;
  
  static constexpr bool has_jacobian_methods = 
    core::meta::is_detected<has_jacobian_callable_with_one_arg, T, state_type>::value &&
    core::meta::is_detected<has_jacobian_callable_with_two_args, T, state_type, matrix_type>::value;

  static constexpr bool is_system = has_residual_methods && has_jacobian_methods;
};


template <typename T, typename U>
struct are_system_compatible {
  typedef system_traits<T> t_traits_type;
  typedef system_traits<U> u_traits_type;

  static constexpr bool are_system = t_traits_type::is_system && u_traits_type::is_system;
  static constexpr bool are_vector_compatible = core::meta::are_vector_compatible<typename t_traits_type::vector_type, typename u_traits_type::vector_type>::value;
  static constexpr bool are_matrix_compatible = core::meta::are_matrix_compatible<typename t_traits_type::matrix_type, typename u_traits_type::matrix_type>::value;
  static constexpr bool value = are_system && are_vector_compatible && are_matrix_compatible;
};


} // end namespace details
} // end namespace solvers

#endif
