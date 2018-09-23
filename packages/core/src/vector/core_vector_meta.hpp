
#ifndef CORE_VECTOR_META_HPP_
#define CORE_VECTOR_META_HPP_

#include "core_vector_traits.hpp"

namespace rompp{
namespace core{
namespace meta {

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

template <typename T, typename enable = void>
struct is_epetra_vector_wrapper : std::false_type {};

template <typename T>
struct is_epetra_vector_wrapper<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_vector &&
       core::details::traits<T>::wrapped_vector_identifier==
       core::details::WrappedVectorIdentifier::Epetra
       >
  >
  : std::true_type{};
//------------------------------------------------------------

  
template <typename T, typename enable = void>
struct is_eigen_vector_wrapper : std::false_type {};

template <typename T>
struct is_eigen_vector_wrapper<
  T, core::meta::enable_if_t<
       core::details::traits<T>::is_vector &&
       core::details::traits<T>::wrapped_vector_identifier==
       core::details::WrappedVectorIdentifier::Eigen
       >
  >
  : std::true_type{};

} // namespace meta
} // namespace core

}//end namespace rompp
#endif
