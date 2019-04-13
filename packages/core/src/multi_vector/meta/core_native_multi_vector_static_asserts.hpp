
#ifndef CORE_NATIVE_MULTI_VECTOR_STATIC_ASSERTS_HPP_
#define CORE_NATIVE_MULTI_VECTOR_STATIC_ASSERTS_HPP_

#include "core_native_eigen_multi_vector_meta.hpp"
#include "core_native_trilinos_multi_vector_meta.hpp"

namespace rompp{ namespace core{

#ifdef HAVE_TRILINOS
#define STATIC_ASSERT_IS_MULTIVECTOR_EPETRA(TYPE) \
  static_assert( core::meta::is_multi_vector_epetra<TYPE>::value, \
		 "THIS_IS_NOT_A_MULTIVECTOR_EPETRA")
#define STATIC_ASSERT_IS_NOT_MULTIVECTOR_EPETRA(TYPE) \
  static_assert( !core::meta::is_multi_vector_epetra<TYPE>::value, \
		 "THIS_IS_A_MULTIVECTOR_EPETRA")

#define STATIC_ASSERT_IS_MULTIVECTOR_TPETRA(TYPE)		  \
  static_assert( core::meta::is_multi_vector_tpetra<TYPE>::value, \
		 "THIS_IS_NOT_A_MULTIVECTOR_TPETRA")
#define STATIC_ASSERT_IS_NOT_MULTIVECTOR_TPETRA(TYPE) \
  static_assert( !core::meta::is_multi_vector_tpetra<TYPE>::value, \
		 "THIS_IS_A_MULTIVECTOR_TPETRA")
#endif

}}//end namespace rompp::core
#endif
