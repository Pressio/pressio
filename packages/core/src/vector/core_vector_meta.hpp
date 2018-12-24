
#ifndef CORE_VECTOR_META_HPP_
#define CORE_VECTOR_META_HPP_

#include "core_vector_traits.hpp"

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_core_vector_wrapper : std::false_type {};

template <typename T>
struct is_core_vector_wrapper<T,
	   typename
	   std::enable_if<
	      core::details::traits<T>::is_vector
	     >::type
	   > : std::true_type{};

//------------------------------------------------------------

#define STATIC_ASSERT_IS_CORE_VECTOR_WRAPPER(TYPE) \
  static_assert( core::meta::is_core_vector_wrapper<TYPE>::value, \
		 "THIS_IS_NOT_A_CORE_VECTOR_WRAPPER")

//------------------------------------------------------------


#ifdef HAVE_TRILINOS
template <typename T, typename enable = void>
struct is_teuchos_serial_dense_vector_wrapper : std::false_type {};

template <typename T>
struct is_teuchos_serial_dense_vector_wrapper<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_vector &&
       core::details::traits<T>::wrapped_vector_identifier==
       core::details::WrappedVectorIdentifier::TeuchosSerialDense
       >
  > : std::true_type{};
#endif

//------------------------------------------------------------

#ifdef HAVE_TRILINOS
template <typename T, typename enable = void>
struct is_epetra_vector_wrapper : std::false_type {};

template <typename T>
struct is_epetra_vector_wrapper<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_vector &&
       core::details::traits<T>::wrapped_vector_identifier==
       core::details::WrappedVectorIdentifier::Epetra
       >
  > : std::true_type{};
#endif

//------------------------------------------------------------

#ifdef HAVE_TRILINOS
template <typename T, typename enable = void>
struct is_tpetra_vector_wrapper : std::false_type {};

template <typename T>
struct is_tpetra_vector_wrapper<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_vector &&
       core::details::traits<T>::wrapped_vector_identifier==
       core::details::WrappedVectorIdentifier::Tpetra
       >
  > : std::true_type{};
#endif

//------------------------------------------------------------


template <typename T, typename enable = void>
struct is_eigen_vector_wrapper : std::false_type {};

template <typename T>
struct is_eigen_vector_wrapper<
  T, core::meta::enable_if_t<
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
//------------------------------------------------------------


#ifdef HAVE_ARMADILLO
template <typename T, typename enable = void>
struct is_armadillo_row_vector_wrapper : std::false_type {};

template <typename T>
struct is_armadillo_row_vector_wrapper<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_vector &&
       core::details::traits<T>::wrapped_vector_identifier==
       core::details::WrappedVectorIdentifier::ArmadilloRow
       >
  > : std::true_type{};


template <typename T, typename enable = void>
struct is_armadillo_column_vector_wrapper : std::false_type {};

template <typename T>
struct is_armadillo_column_vector_wrapper<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_vector &&
       core::details::traits<T>::wrapped_vector_identifier==
       core::details::WrappedVectorIdentifier::ArmadilloCol
       >
  > : std::true_type{};
#endif // HAVE_ARMADILLO
//------------------------------------------------------------


#ifdef HAVE_BLAZE
template <typename T, typename enable = void>
struct is_blaze_dynamic_vector_wrapper : std::false_type {};

template <typename T>
struct is_blaze_dynamic_vector_wrapper<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_vector &&
       core::details::traits<T>::wrapped_vector_identifier==
       core::details::WrappedVectorIdentifier::BlazeDynamic
       >
  > : std::true_type{};


template <typename T, typename enable = void>
struct is_blaze_static_vector_wrapper : std::false_type {};

template <typename T>
struct is_blaze_static_vector_wrapper<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_vector &&
       core::details::traits<T>::wrapped_vector_identifier==
       core::details::WrappedVectorIdentifier::BlazeStatic
       >
  > : std::true_type{};
#endif // HAVE_BLAZE

//------------------------------------------------------------


template <typename T,
    typename enable = void>
struct is_admissible_vec_for_dist_expression : std::false_type{};

template <typename T>
struct is_admissible_vec_for_dist_expression<T,
      core::meta::enable_if_t<
  core::meta::is_core_vector_wrapper<T>::value &&
  !core::details::traits<T>::is_shared_mem
      >> : std::true_type{};

///-----------------------------------------------------

template <typename T,
    typename enable = void>
struct is_admissible_vec_for_sharedmem_expression : std::false_type{};

template <typename T>
struct is_admissible_vec_for_sharedmem_expression<T,
      core::meta::enable_if_t<
  core::meta::is_core_vector_wrapper<T>::value &&
  core::details::traits<T>::is_shared_mem
      >> : std::true_type{};

///-----------------------------------------------------



}}}//end namespace rompp::core::meta
#endif
