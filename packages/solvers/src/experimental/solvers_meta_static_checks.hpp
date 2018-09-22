
#ifndef SOLVERS_META_META_STATIC_CHECKS_HPP
#define SOLVERS_META_META_STATIC_CHECKS_HPP

#include <type_traits>

#include "../../../core/src/matrix/core_matrix_traits.hpp"
#include "../../../core/src/vector/core_vector_traits.hpp"


namespace solvers{
namespace meta {

/**
 * @brief Check whether two matrices are compatible.
 *
 * @section DESCRIPTION
 *
 * Two matrices are compatible iff the following two conditions are satisfied:
 *   1. the underlying matrix representation is the same; this implies that the
 *      underlying linear algebra library is also the same;
 *   2. the matrices have the same structure (sparse or dense).
 */
template <typename T, typename U>
struct are_matrix_compatible {
  static constexpr bool valid_matrix = core::details::traits<T>::wrapped_package_identifier != core::details::WrappedPackageIdentifier::Undefined;
  static constexpr bool value = valid_matrix && (core::details::traits<T>::wrapped_matrix_identifier == core::details::traits<U>::wrapped_matrix_identifier);
};


/**
 * @brief Check whether a vector and a matrix are compatible.
 *
 * @section DESCRIPTION
 *
 * A vector and a matrix are compatible iff the underlying linear algebra package
 * used to represent them is the same.
 */
template <typename T, typename U>
struct are_vector_matrix_compatible {
  static constexpr bool valid_vector = core::details::traits<T>::wrapped_package_identifier != core::details::WrappedPackageIdentifier::Undefined;
  static constexpr bool valid_matrix = core::details::traits<U>::wrapped_package_identifier != core::details::WrappedPackageIdentifier::Undefined;
  static constexpr bool value = valid_vector && valid_matrix && (core::details::traits<T>::wrapped_package_identifier == core::details::traits<U>::wrapped_package_identifier);
};


/**
 * @brief Check whether two vectors are compatible.
 *
 * @section DESCRIPTION
 *
 * Two vectors are compatible iff their underlying structure is the same and they
 * are represented using the same underlying linear algebra package.
 */
template <
  typename T,
  typename U,
  typename Enable = void
>
struct are_vector_compatible;


/**
 * Either one or both arguments are not vectors.
 */
template <
  typename T,
  typename U
>
struct are_vector_compatible<
  T,
  U,
  typename std::enable_if<
    !core::details::traits<T>::is_vector || !core::details::traits<U>::is_vector,
    void
  >::type
> {
  static constexpr bool value = false;
};


/**
 * Both arguments are vectors.
 */
template <
  typename T,
  typename U
>
struct are_vector_compatible<
  T,
  U,
  typename std::enable_if<
    core::details::traits<T>::is_vector && core::details::traits<U>::is_vector,
    void
  >::type
> {

  static constexpr bool valid_vector =
    core::details::traits<T>::wrapped_package_identifier !=
    core::details::WrappedPackageIdentifier::Undefined;

  static constexpr bool same_type =
    valid_vector &&
    core::details::traits<T>::wrapped_package_identifier ==
    core::details::traits<U>::wrapped_package_identifier;

  static constexpr bool value = same_type &&
    (core::details::traits<T>::is_dynamic ||
     core::details::traits<U>::is_dynamic ||
     core::details::traits<T>::rows == core::details::traits<U>::rows);
};


} // end namespace meta
} // end namespace solvers

#endif
