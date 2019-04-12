
#ifndef CORE_IS_VECTOR_WRAPPER_EIGEN_HPP_
#define CORE_IS_VECTOR_WRAPPER_EIGEN_HPP_

#include "../core_vector_traits.hpp"

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_vector_wrapper_eigen : std::false_type {};

template <typename T>
struct is_vector_wrapper_eigen<
  T, ::rompp::mpl::enable_if_t<
       core::details::traits<T>::is_vector &&
       (core::details::traits<T>::wrapped_vector_identifier==
        core::details::WrappedVectorIdentifier::EigenColStatic or
        core::details::traits<T>::wrapped_vector_identifier==
        core::details::WrappedVectorIdentifier::EigenColDynamic or
        core::details::traits<T>::wrapped_vector_identifier==
        core::details::WrappedVectorIdentifier::EigenRowStatic or
        core::details::traits<T>::wrapped_vector_identifier==
        core::details::WrappedVectorIdentifier::EigenRowDynamic)
       >
  > : std::true_type{};

}}}//end namespace rompp::core::meta
#endif
