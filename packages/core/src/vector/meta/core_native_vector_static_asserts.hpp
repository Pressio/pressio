
#ifndef CORE_NATIVE_VECTOR_STATIC_ASSERTS_HPP_
#define CORE_NATIVE_VECTOR_STATIC_ASSERTS_HPP_

#include "core_native_blaze_vector_meta.hpp"
#include "core_native_eigen_vector_meta.hpp"
#include "core_native_stdlib_vector_meta.hpp"
#include "core_native_epetra_vector_meta.hpp"
#include "core_native_tpetra_vector_meta.hpp"
#include "core_native_teuchos_vector_meta.hpp"
#include "core_native_kokkos_vector_meta.hpp"

namespace rompp{ namespace core{

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

}}//end namespace rompp::core
#endif
