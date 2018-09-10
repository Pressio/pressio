
#ifndef SOLVERS_META_META_STATIC_CHECKS_HPP
#define SOLVERS_META_META_STATIC_CHECKS_HPP

#include "matrix/core_matrix_traits_exp.hpp"
#include "vector/core_vector_traits_exp.hpp"

namespace solvers{
namespace meta {

template <typename T, typename U>
struct are_matrix_compatible {
  static constexpr bool valid_matrix = core::details::matrix_traits<T>::wrapped_package_identifier != core::details::WrappedPackageIdentifier::Undefined;
  static constexpr bool same_type = core::details::matrix_traits<T>::wrapped_package_identifier == core::details::matrix_traits<U>::wrapped_package_identifier;
  static constexpr bool same_structure = core::details::matrix_traits<T>::is_sparse ? core::details::matrix_traits<U>::is_sparse : !core::details::matrix_traits<U>::is_sparse;
  static constexpr bool value = valid_matrix && same_type && same_structure;
};
//----------------------------------------------------


// template <typename T, typename U>
// struct are_vector_compatible {
//   static constexpr bool valid_vector = core::details::vector_traits<T>::wrapped_package_identifier != core::details::WrappedPackageIdentifier::Undefined;
//   static constexpr bool same_type = core::details::vector_traits<T>::wrapped_package_identifier == core::details::vector_traits<U>::wrapped_package_identifier;
//   static constexpr bool same_structure = core::details::vector_traits<T>::is_dynamic || core::details::vector_traits<U>::is_dynamic || core::details::vector_traits<T>::rows == core::details::vector_traits<U>::rows;
//   static constexpr bool value = valid_vector && same_type && same_structure;
// };
// //----------------------------------------------------


template <typename T, typename U>
struct are_vector_matrix_compatible {
  static constexpr bool valid_vector = core::details::vector_traits<T>::wrapped_package_identifier != core::details::WrappedPackageIdentifier::Undefined;
  static constexpr bool valid_matrix = core::details::matrix_traits<U>::wrapped_package_identifier != core::details::WrappedPackageIdentifier::Undefined; 
  static constexpr bool same_type = core::details::vector_traits<T>::wrapped_package_identifier == core::details::matrix_traits<U>::wrapped_package_identifier;
  static constexpr bool value = valid_matrix && valid_vector && same_type;
};
//----------------------------------------------------


template <typename T, typename U>
struct same_vector_structure {
  static constexpr bool value = core::details::vector_traits<T>::is_dynamic ||
    core::details::vector_traits<U>::is_dynamic ||
    core::details::vector_traits<T>::rows == core::details::vector_traits<U>::rows;
};
//----------------------------------------------------


template <typename T, typename U>
struct are_vector_compatible {

  static constexpr bool valid_vector =
    core::details::vector_traits<T>::wrapped_package_identifier !=
    core::details::WrappedPackageIdentifier::Undefined;

  static constexpr bool same_type =
    core::details::vector_traits<T>::wrapped_package_identifier ==
    core::details::vector_traits<U>::wrapped_package_identifier;

  static constexpr bool same_structure =
    core::details::vector_traits<T>::is_dynamic ||
    core::details::vector_traits<U>::is_dynamic ||
    core::details::vector_traits<T>::rows == core::details::vector_traits<U>::rows;

  static constexpr bool value =
    valid_vector && same_type && same_structure;
};
//----------------------------------------------------
  

  
} // end namespace meta	
} // end namespace solvers

#endif
