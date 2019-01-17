
#ifndef CORE_UTILS_STATIC_ASSERT_DEFINITIONS_HPP_
#define CORE_UTILS_STATIC_ASSERT_DEFINITIONS_HPP_

#include "../meta/core_native_vector_meta.hpp"
#include "../meta/core_native_multi_vector_meta.hpp"
#include "../meta/core_native_matrix_meta.hpp"

namespace rompp{ namespace core{

//////////////////////
// VECTOR
/////////////////////

#define STATIC_ASSERT_IS_VECTOR_EIGEN(TYPE) \
  static_assert( core::meta::is_vector_eigen<TYPE>::value, \
		 "THIS_IS_NOT_A_VECTOR_FROM_EIGEN")
#define STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(TYPE) \
  static_assert( !core::meta::is_vector_eigen<TYPE>::value, \
		 "THIS_IS_A_VECTOR_FROM_EIGEN")

#define STATIC_ASSERT_IS_VECTOR_STDLIB(TYPE) \
  static_assert( core::meta::is_vector_stdlib<TYPE>::value, \
		 "THIS_IS_NOT_A_STDLIB_VECTOR")
#define STATIC_ASSERT_IS_NOT_VECTOR_STDLIB(TYPE) \
  static_assert( !core::meta::is_vector_stdlib<TYPE>::value, \
		 "THIS_IS_A_STDLIB_VECTOR")

#ifdef HAVE_TRILINOS
#define STATIC_ASSERT_IS_VECTOR_EPETRA(TYPE) \
  static_assert( core::meta::is_vector_epetra<TYPE>::value, \
		 "THIS_IS_NOT_A_VECTOR_EPETRA")
#define STATIC_ASSERT_IS_NOT_VECTOR_EPETRA(TYPE) \
  static_assert( !core::meta::is_vector_epetra<TYPE>::value, \
		 "THIS_IS_A_VECTOR_EPETRA")

#define STATIC_ASSERT_IS_VECTOR_TPETRA(TYPE) \
  static_assert( core::meta::is_vector_tpetra<TYPE>::value, \
		 "THIS_IS_NOT_A_VECTOR_TPETRA")
#define STATIC_ASSERT_IS_NOT_VECTOR_TPETRA(TYPE) \
  static_assert( !core::meta::is_vector_tpetra<TYPE>::value, \
		 "THIS_IS_A_VECTOR_TPETRA")

#define STATIC_ASSERT_IS_VECTOR_KOKKOS(TYPE) \
  static_assert( core::meta::is_vector_kokkos<TYPE>::value, \
		 "THIS_IS_NOT_A_VECTOR_KOKKOS")
#define STATIC_ASSERT_IS_NOT_VECTOR_KOKKOS(TYPE)\
  static_assert( !core::meta::is_vector_kokkos<TYPE>::value, \
		 "THIS_IS_A_VECTOR_KOKKOS")
#endif


#ifdef HAVE_BLAZE
#define STATIC_ASSERT_IS_VECTOR_BLAZE(TYPE)				\
  static_assert( core::meta::is_static_vector_blaze<TYPE>::value ||	\
                 core::meta::is_dynamic_vector_blaze<TYPE>::value,	\
		 "THIS_IS_NOT_A_VECTOR_BLAZE")
#define STATIC_ASSERT_IS_NOT_VECTOR_BLAZE(TYPE)				\
  static_assert( !core::meta::is_static_vector_blaze<TYPE>::value &&	\
                 !core::meta::is_dynamic_vector_blaze<TYPE>::value,	\
		 "THIS_IS_A_VECTOR_BLAZE")
#endif

#ifdef HAVE_ARMADILLO
#define STATIC_ASSERT_IS_VECTOR_ARMADILLO(TYPE)			\
  static_assert( core::meta::is_vector_armadillo<TYPE>::value,	\
		 "THIS_IS_NOT_A_VECTOR_ARMADILLO")
#define STATIC_ASSERT_IS_NOT_VECTOR_ARMADILLO(TYPE)		\
  static_assert( !core::meta::is_vector_armadillo<TYPE>::value,	\
		 "THIS_IS_A_VECTOR_ARMADILLO")
#endif


////////////////////////
// MULTI VECTOR
///////////////////////

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


////////////////////////
// MATRIX
///////////////////////

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
