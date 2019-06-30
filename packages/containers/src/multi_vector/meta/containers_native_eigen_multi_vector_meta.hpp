
#ifndef CONTAINERS_NATIVE_EIGEN_MULTIVECTOR_META_HPP_
#define CONTAINERS_NATIVE_EIGEN_MULTIVECTOR_META_HPP_

#include "../../meta/containers_meta_basic.hpp"
#include "../../matrix/meta/containers_native_eigen_matrix_meta.hpp"

namespace rompp{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_dynamic_multi_vector_eigen : std::false_type {};

template <typename T>
struct is_dynamic_multi_vector_eigen<T,
  typename
   std::enable_if<
    is_dense_dynamic_matrix_eigen<T>::value
   >::type
  > : std::true_type{};


}}}//end namespace rompp::containers::meta
#endif
