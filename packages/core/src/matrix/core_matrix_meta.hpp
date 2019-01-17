\
#ifndef CORE_MATRIX_META_HPP_
#define CORE_MATRIX_META_HPP_

#include "core_matrix_traits.hpp"

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_core_matrix_wrapper : std::false_type {};

template <typename T>
struct is_core_matrix_wrapper< T,
		       typename
		       std::enable_if<
			 core::details::traits<T>::is_matrix
			 >::type
		       > : std::true_type{};


#define STATIC_ASSERT_IS_CORE_MATRIX_WRAPPER(TYPE) \
  static_assert( core::meta::is_core_matrix_wrapper<TYPE>::value, \
		 "THIS_IS_NOT_A_CORE_MATRIX_WRAPPER")
//------------------------------------------------------------


#ifdef HAVE_TRILINOS
template <typename T, typename enable = void>
struct is_dense_matrix_wrapper_teuchos : std::false_type {};

template <typename T>
struct is_dense_matrix_wrapper_teuchos<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_matrix &&
       core::details::traits<T>::wrapped_matrix_identifier==
       core::details::WrappedMatrixIdentifier::TeuchosSerialDense
       >
  > : std::true_type{};
#endif

//------------------------------------------------------------

#ifdef HAVE_TRILINOS
template <typename T, typename enable = void>
struct is_dense_matrix_wrapper_epetra : std::false_type {};

template <typename T>
struct is_dense_matrix_wrapper_epetra<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_matrix &&
       core::details::traits<T>::wrapped_matrix_identifier==
       core::details::WrappedMatrixIdentifier::DenseEpetra
       >
  >
  : std::true_type{};
#endif
//------------------------------------------------------------


#ifdef HAVE_TRILINOS
template <typename T, typename enable = void>
struct is_sparse_matrix_wrapper_epetra : std::false_type {};

template <typename T>
struct is_sparse_matrix_wrapper_epetra<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_matrix &&
       (core::details::traits<T>::wrapped_matrix_identifier==
	core::details::WrappedMatrixIdentifier::CrsEpetra)
       >
  >
  : std::true_type{};
  #endif
//------------------------------------------------------------


#ifdef HAVE_TRILINOS
template <typename T, typename enable = void>
struct is_sparse_matrix_wrapper_tpetra : std::false_type {};

template <typename T>
struct is_sparse_matrix_wrapper_tpetra<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_matrix &&
       core::details::traits<T>::wrapped_matrix_identifier==
       core::details::WrappedMatrixIdentifier::SparseTpetra
       >
  >
  : std::true_type{};
#endif
//------------------------------------------------------------


template <typename T, typename enable = void>
struct is_dense_matrix_wrapper_eigen : std::false_type {};

template <typename T>
struct is_dense_matrix_wrapper_eigen<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_matrix &&
       (core::details::traits<T>::wrapped_matrix_identifier==
	core::details::WrappedMatrixIdentifier::DenseEigen)
       >
  >
  : std::true_type{};
//------------------------------------------------------------

template <typename T, typename enable = void>
struct is_sparse_matrix_wrapper_eigen : std::false_type {};

template <typename T>
struct is_sparse_matrix_wrapper_eigen<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_matrix &&
       (core::details::traits<T>::wrapped_matrix_identifier==
       core::details::WrappedMatrixIdentifier::SparseEigen)
       >
  >
  : std::true_type{};
//------------------------------------------------------------


template <typename T, typename enable = void>
struct is_matrix_wrapper_eigen : std::false_type {};

template <typename T>
struct is_matrix_wrapper_eigen<
  T, core::meta::enable_if_t<
       is_sparse_matrix_wrapper_eigen<T>::value or
       is_dense_matrix_wrapper_eigen<T>::value
       >
  >
  : std::true_type{};
//------------------------------------------------------------



}}}//end namespace rompp::core::meta
#endif
