
#ifndef CONTAINERS_NATIVE_MULTI_VECTOR_STATIC_ASSERTS_HPP_
#define CONTAINERS_NATIVE_MULTI_VECTOR_STATIC_ASSERTS_HPP_

#include "./meta/containers_native_eigen_multi_vector_meta.hpp"
#include "./meta/containers_native_epetra_multi_vector_meta.hpp"
#include "./meta/containers_native_tpetra_multi_vector_meta.hpp"
#include "./meta/containers_native_tpetra_block_multi_vector_meta.hpp"

namespace rompp{ namespace containers{

#ifdef HAVE_TRILINOS
#define STATIC_ASSERT_IS_MULTIVECTOR_EPETRA(TYPE) \
  static_assert( containers::meta::is_multi_vector_epetra<TYPE>::value, \
		 "THIS_IS_NOT_A_MULTIVECTOR_EPETRA")
#define STATIC_ASSERT_IS_NOT_MULTIVECTOR_EPETRA(TYPE) \
  static_assert( !containers::meta::is_multi_vector_epetra<TYPE>::value, \
		 "THIS_IS_A_MULTIVECTOR_EPETRA")

#define STATIC_ASSERT_IS_MULTIVECTOR_TPETRA(TYPE)		  \
  static_assert( containers::meta::is_multi_vector_tpetra<TYPE>::value, \
		 "THIS_IS_NOT_A_MULTIVECTOR_TPETRA")
#define STATIC_ASSERT_IS_NOT_MULTIVECTOR_TPETRA(TYPE) \
  static_assert( !containers::meta::is_multi_vector_tpetra<TYPE>::value, \
		 "THIS_IS_A_MULTIVECTOR_TPETRA")

#define STATIC_ASSERT_IS_MULTIVECTOR_TPETRA_BLOCK(TYPE)		  \
  static_assert( containers::meta::is_multi_vector_tpetra_block<TYPE>::value, \
		 "THIS_IS_NOT_A_MULTIVECTOR_TPETRA_BLOCK")
#define STATIC_ASSERT_IS_NOT_MULTIVECTOR_TPETRA_BLOCK(TYPE) \
  static_assert( !containers::meta::is_multi_vector_tpetra_block<TYPE>::value, \
		 "THIS_IS_A_MULTIVECTOR_TPETRA_BLOCK")

#endif

}}//end namespace rompp::containers
#endif
