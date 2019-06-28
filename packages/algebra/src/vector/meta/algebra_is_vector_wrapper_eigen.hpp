
#ifndef ALGEBRA_IS_VECTOR_WRAPPER_EIGEN_HPP_
#define ALGEBRA_IS_VECTOR_WRAPPER_EIGEN_HPP_

#include "../algebra_vector_traits.hpp"

namespace rompp{ namespace algebra{ namespace meta {

template <typename T, typename enable = void>
struct is_vector_wrapper_eigen : std::false_type {};

template <typename T>
struct is_vector_wrapper_eigen<
  T, ::rompp::mpl::enable_if_t<
       algebra::details::traits<T>::is_vector &&
       (algebra::details::traits<T>::wrapped_vector_identifier==
        algebra::details::WrappedVectorIdentifier::EigenColStatic or
        algebra::details::traits<T>::wrapped_vector_identifier==
        algebra::details::WrappedVectorIdentifier::EigenColDynamic or
        algebra::details::traits<T>::wrapped_vector_identifier==
        algebra::details::WrappedVectorIdentifier::EigenRowStatic or
        algebra::details::traits<T>::wrapped_vector_identifier==
        algebra::details::WrappedVectorIdentifier::EigenRowDynamic)
       >
  > : std::true_type{};

}}}//end namespace rompp::algebra::meta
#endif
