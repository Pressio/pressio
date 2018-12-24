
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
struct is_teuchos_serial_dense_matrix_wrapper : std::false_type {};

template <typename T>
struct is_teuchos_serial_dense_matrix_wrapper<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_matrix &&
       core::details::traits<T>::wrapped_matrix_identifier==
       core::details::WrappedVectorIdentifier::TeuchosSerialDense
       >
  > : std::true_type{};
#endif

//------------------------------------------------------------

#ifdef HAVE_TRILINOS
template <typename T, typename enable = void>
struct is_epetra_dense_matrix_wrapper : std::false_type {};

template <typename T>
struct is_epetra_dense_matrix_wrapper<
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
struct is_epetra_sparse_matrix_wrapper : std::false_type {};

template <typename T>
struct is_epetra_sparse_matrix_wrapper<
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
struct is_tpetra_sparse_matrix_wrapper : std::false_type {};

template <typename T>
struct is_tpetra_sparse_matrix_wrapper<
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
struct is_eigen_dense_matrix_wrapper : std::false_type {};

template <typename T>
struct is_eigen_dense_matrix_wrapper<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_matrix &&
       (core::details::traits<T>::wrapped_matrix_identifier==
	core::details::WrappedMatrixIdentifier::DenseEigen)
       >
  >
  : std::true_type{};
//------------------------------------------------------------

template <typename T, typename enable = void>
struct is_eigen_sparse_matrix_wrapper : std::false_type {};

template <typename T>
struct is_eigen_sparse_matrix_wrapper<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_matrix &&
       (core::details::traits<T>::wrapped_matrix_identifier==
       core::details::WrappedMatrixIdentifier::SparseEigen)
       >
  >
  : std::true_type{};
//------------------------------------------------------------


}}}//end namespace rompp::core::meta
#endif
