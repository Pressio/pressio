
#ifndef CORE_NATIVE_MULTIVECTOR_META_MULTIVECTOR_META_HPP_
#define CORE_NATIVE_MULTIVECTOR_META_MULTIVECTOR_META_HPP_

#include "core_meta_basic.hpp"
#include "core_native_vector_meta.hpp"
#include "core_native_matrix_meta.hpp"

namespace rompp{ namespace core{ namespace meta {
 
#ifdef HAVE_TRILINOS
template <typename T, typename enable = void>
struct is_multi_vector_epetra : std::false_type {};

template <typename T>
struct is_multi_vector_epetra<T,
  typename
   std::enable_if<
    std::is_same<T,Epetra_MultiVector>::value
   >::type
  > : std::true_type{};
#endif
//-------------------------------------------------


#ifdef HAVE_TRILINOS  
template <typename T, typename enable = void>
struct is_multi_vector_tpetra : std::false_type {};

template <typename T>
struct is_multi_vector_tpetra<T,
      typename
      std::enable_if<
	std::is_same<T,
		     Tpetra::MultiVector<
		       typename T::impl_scalar_type,
		       typename T::local_ordinal_type,
		       typename T::global_ordinal_type,
		       typename T::node_type
		       >
		     >::value 
	>::type
      > : std::true_type{};
#endif 
//--------------------------------------------
      

template <typename T, typename enable = void>
struct is_multi_vector_eigen_dynamic : std::false_type {};

template <typename T>
struct is_multi_vector_eigen_dynamic<T,
  typename
   std::enable_if<
    is_matrix_dense_sharedmem_eigen_dynamic<T>::value
   >::type
  > : std::true_type{};
      

}}}//end namespace rompp::core::meta
#endif
