
#ifndef CORE_NATIVE_EIGEN_MULTIVECTOR_META_HPP_
#define CORE_NATIVE_EIGEN_MULTIVECTOR_META_HPP_

#include "../core_meta_basic.hpp"
#include "../native_matrix/core_native_eigen_matrix_meta.hpp"

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_dynamic_multi_vector_eigen : std::false_type {};

template <typename T>
struct is_dynamic_multi_vector_eigen<T,
  typename
   std::enable_if<
    is_dense_dynamic_matrix_eigen<T>::value
   >::type
  > : std::true_type{};


}}}//end namespace rompp::core::meta
#endif
