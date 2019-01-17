
#ifndef CORE_MULTI_VECTOR_META_HPP_
#define CORE_MULTI_VECTOR_META_HPP_

#include "core_multi_vector_traits.hpp"

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_core_multi_vector_wrapper : std::false_type {};

template <typename T>
struct is_core_multi_vector_wrapper<T,
	   typename
	   std::enable_if<
	     core::details::traits<T>::is_multi_vector
	     >::type
	   > : std::true_type{};
//------------------------------------------------------------

#define STATIC_ASSERT_IS_CORE_MULTI_VECTOR_WRAPPER(TYPE) \
  static_assert(core::meta::is_core_multi_vector_wrapper<TYPE>::value,\
		"THIS_IS_NOT_A_CORE_MULTI_VECTOR_WRAPPER")

//------------------------------------------------------------


#ifdef HAVE_TRILINOS
template <typename T, typename enable = void>
struct is_multi_vector_wrapper_epetra : std::false_type {};

template <typename T>
struct is_multi_vector_wrapper_epetra<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_multi_vector &&
       core::details::traits<T>::wrapped_multi_vector_identifier==
       core::details::WrappedMultiVectorIdentifier::Epetra
       >
  >
  : std::true_type{};
#endif
//------------------------------------------------------------


#ifdef HAVE_TRILINOS
template <typename T, typename enable = void>
struct is_multi_vector_wrapper_tpetra : std::false_type {};

template <typename T>
struct is_multi_vector_wrapper_tpetra<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_multi_vector &&
       core::details::traits<T>::wrapped_multi_vector_identifier==
       core::details::WrappedMultiVectorIdentifier::Tpetra
       >
  >
  : std::true_type{};
#endif
//------------------------------------------------------------


template <typename T, typename enable = void>
struct is_multi_vector_wrapper_eigen : std::false_type {};

template <typename T>
struct is_multi_vector_wrapper_eigen<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_multi_vector &&
       core::details::traits<T>::wrapped_multi_vector_identifier==
       core::details::WrappedMultiVectorIdentifier::Eigen
       >
  >
  : std::true_type{};

}}}//end namespace rompp::core::meta
#endif
