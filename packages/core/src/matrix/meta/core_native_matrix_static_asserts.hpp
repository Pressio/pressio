
#ifndef CORE_NATIVE_MATRIX_STATIC_ASSERTS_HPP_
#define CORE_NATIVE_MATRIX_STATIC_ASSERTS_HPP_

#include "core_native_eigen_matrix_meta.hpp"
#include "core_native_stdlib_matrix_meta.hpp"
#include "core_native_trilinos_matrix_meta.hpp"

namespace rompp{ namespace core{

#define STATIC_ASSERT_IS_MATRIX_DENSE_SHAREDMEM_EIGEN(TYPE) \
  static_assert( core::meta::is_dense_matrix_eigen<TYPE>::value,	\
		 "THIS_IS_NOT_A_MATRIX_DENSE_SHAREDMEM_EIGEN")
#define STATIC_ASSERT_IS_NOT_MATRIX_DENSE_SHAREDMEM_EIGEN(TYPE) \
  static_assert( !core::meta::is_dense_matrix_eigen<TYPE>::value, \
		 "THIS_IS_A_MATRIX_DENSE_SHAREDMEM_EIGEN")


#define STATIC_ASSERT_IS_MATRIX_SPARSE_SHAREDMEM_EIGEN(TYPE) \
  static_assert( core::meta::is_sparse_matrix_eigen<TYPE>::value,	\
		 "THIS_IS_NOT_A_MATRIX_SPARSE_SHAREDMEM_EIGEN")
#define STATIC_ASSERT_IS_NOT_MATRIX_SPARSE_SHAREDMEM_EIGEN(TYPE) \
  static_assert( !core::meta::is_sparse_matrix_eigen<TYPE>::value, \
		 "THIS_IS_A_MATRIX_SPARSE_SHAREDMEM_EIGEN")


#define STATIC_ASSERT_IS_MATRIX_DENSE_SHAREDMEM_STDLIB(TYPE)	      \
  static_assert( core::meta::is_dense_matrix_stdlib<TYPE>::value, \
		 "THIS_IS_NOT_A_MATRIX_DENSE_SHAREDMEM_STDLIB")
#define STATIC_ASSERT_IS_NOT_MATRIX_DENSE_SHAREDMEM_STDLIB(TYPE) \
  static_assert( !core::meta::is_dense_matrix_stdlib<TYPE>::value, \
		 "THIS_IS_A_MATRIX_DENSE_SHAREDMEM_STDLIB")


#ifdef HAVE_TRILINOS
#define STATIC_ASSERT_IS_MATRIX_SPARSE_DISTRIBUTED_EPETRA(TYPE)	      \
  static_assert( core::meta::is_sparse_matrix_epetra<TYPE>::value, \
		 "THIS_IS_NOT_A_MATRIX_SPARSE_DIST_EPETRA")
#define STATIC_ASSERT_IS_NOT_MATRIX_SPARSE_DISTRIBUTED_EPETRA(TYPE) \
  static_assert( !core::meta::is_sparse_matrix_epetra<TYPE>::value, \
		 "THIS_IS_A_MATRIX_SPARSE_DIST_EPETRA")


#define STATIC_ASSERT_IS_MATRIX_SPARSE_DISTRIBUTED_TPETRA(TYPE)     \
  static_assert( core::meta::is_sparse_matrix_tpetra<TYPE>::value, \
     "THIS_IS_NOT_A_MATRIX_SPARSE_DIST_TPETRA")
#define STATIC_ASSERT_IS_NOT_MATRIX_SPARSE_DISTRIBUTED_TPETRA(TYPE) \
  static_assert( !core::meta::is_sparse_matrix_tpetra<TYPE>::value, \
     "THIS_IS_A_MATRIX_SPARSE_DIST_TPETRA")


#define STATIC_ASSERT_IS_MATRIX_DENSE_DISTRIBUTED_EPETRA(TYPE)	      \
  static_assert( core::meta::is_dense_matrix_epetra<TYPE>::value, \
		 "THIS_IS_NOT_A_MATRIX_DENSE_DIST_EPETRA")
#define STATIC_ASSERT_IS_NOT_MATRIX_DENSE_DISTRIBUTED_EPETRA(TYPE) \
  static_assert( !core::meta::is_dense_matrix_epetra<TYPE>::value, \
		 "THIS_IS_A_MATRIX_DENSE_DIST_EPETRA")
#endif


}}//end namespace rompp::core
#endif
