
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_IS_MULTI_VECTOR_WRAPPER_TPETRA_BLOCK_HPP_
#define CONTAINERS_IS_MULTI_VECTOR_WRAPPER_TPETRA_BLOCK_HPP_

#include "../containers_multi_vector_traits.hpp"

namespace rompp{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_multi_vector_wrapper_tpetra_block : std::false_type {};

template <typename T>
struct is_multi_vector_wrapper_tpetra_block<
  T, ::rompp::mpl::enable_if_t<
       containers::details::traits<T>::is_multi_vector &&
       containers::details::traits<T>::wrapped_multi_vector_identifier==
       containers::details::WrappedMultiVectorIdentifier::TpetraBlock
       >
  >
  : std::true_type{};

}}}//end namespace rompp::containers::meta
#endif
#endif
