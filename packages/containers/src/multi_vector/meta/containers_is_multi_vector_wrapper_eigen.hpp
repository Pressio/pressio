
#ifndef CONTAINERS_IS_MULTI_VECTOR_WRAPPER_EIGEN_HPP_
#define CONTAINERS_IS_MULTI_VECTOR_WRAPPER_EIGEN_HPP_

#include "../containers_multi_vector_traits.hpp"

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_multi_vector_wrapper_eigen : std::false_type {};

template <typename T>
struct is_multi_vector_wrapper_eigen<
  T, ::pressio::mpl::enable_if_t<
       containers::details::traits<T>::is_multi_vector &&
       containers::details::traits<T>::wrapped_multi_vector_identifier==
       containers::details::WrappedMultiVectorIdentifier::Eigen
       >
  >
  : std::true_type{};

}}}//end namespace pressio::containers::meta
#endif
