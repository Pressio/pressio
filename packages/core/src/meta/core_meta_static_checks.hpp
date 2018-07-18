
#ifndef CORE_META_STATIC_CHECKS_HPP
#define CORE_META_STATIC_CHECKS_HPP

#include "matrix/core_matrix_traits_exp.hpp"
#include "vector/core_vector_traits_exp.hpp"

namespace core {
namespace meta {

template <typename T, typename U>
struct are_matrix_compatible {
  static constexpr bool valid_matrix = details::matrix_traits<T>::matrix_class != details::WrappedClass::Undefined;
  static constexpr bool same_type = details::matrix_traits<T>::matrix_class == details::matrix_traits<U>::matrix_class;
  static constexpr bool same_structure = details::matrix_traits<T>::is_sparse ? details::matrix_traits<U>::is_sparse : !details::matrix_traits<U>::is_sparse;
  static constexpr bool value = valid_matrix && same_type && same_structure;
};


template <typename T, typename U>
struct are_vector_compatible {
  static constexpr bool valid_vector = details::vector_traits<T>::vector_class != details::WrappedClass::Undefined;
  static constexpr bool same_type = details::vector_traits<T>::vector_class == details::vector_traits<U>::vector_class;
  static constexpr bool same_structure = details::vector_traits<T>::is_dynamic || details::vector_traits<U>::is_dynamic || details::vector_traits<T>::rows == details::vector_traits<U>::rows;
  static constexpr bool value = valid_vector && same_type && same_structure;
};


template <typename T, typename U>
struct are_vector_matrix_compatible {
  static constexpr bool valid_vector = details::vector_traits<T>::vector_class != details::WrappedClass::Undefined;
  static constexpr bool valid_matrix = details::matrix_traits<U>::matrix_class != details::WrappedClass::Undefined; 
  static constexpr bool same_type = details::vector_traits<T>::vector_class == details::matrix_traits<U>::matrix_class;
  static constexpr bool value = valid_matrix && valid_vector && same_type;
};


} // end namespace meta	
} // end namespace core

#endif