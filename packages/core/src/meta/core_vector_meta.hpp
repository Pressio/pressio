
#ifndef CORE_VECTOR_META_VECTOR_META_HPP_
#define CORE_VECTOR_META_VECTOR_META_HPP_

#include "core_meta_basic.hpp"
#include "core_meta_detection_idiom.hpp"
#include <vector>
#include <Eigen/Dense>
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include <Kokkos_Core.hpp>

namespace core{
namespace meta {

  
template <typename T, typename enable = void>
struct is_vector_eigen : std::false_type {};

template <typename T>
struct is_vector_eigen< T,
     typename
     std::enable_if<
       std::is_same<T,
	 Eigen::Matrix<typename T::Scalar,
		       1,
		       T::ColsAtCompileTime
		       >
	 >::value
       >::type
     > : std::true_type{};

template <typename T>
struct is_vector_eigen< T,
      typename
      std::enable_if<
	std::is_same<T,
	  Eigen::Matrix<typename T::Scalar,
			T::RowsAtCompileTime,
			1>
	  >::value
	>::type
      > : std::true_type{};
//----------------------------------------------------------------------

  
template <typename T, typename enable = void>
struct is_vector_stdlib : std::false_type {};

template <typename T>
struct is_vector_stdlib<T,
      typename
      std::enable_if<
	std::is_same<T,
	  std::vector<typename T::value_type>
	  >::value &&
	// we do not want to have Vector<Vector<...>>
	// so we need to check that the T::value_type is a
	// scalar type or integral type or complex
	(std::is_floating_point<typename T::value_type>::value ||
	 std::is_integral<typename T::value_type>::value ||
	 is_std_complex<typename T::value_type>::value
	 )
	>::type
      > : std::true_type{};
//----------------------------------------------------------------------

  
template <typename T, typename enable = void>
struct is_vector_epetra : std::false_type {};

template <typename T>
struct is_vector_epetra<T,
      typename
      std::enable_if<
	std::is_same<T,Epetra_Vector>::value 
	>::type
      > : std::true_type{};

//----------------------------------------------------------------------


template <typename T, typename enable = void>
struct is_vector_kokkos : std::false_type {};

template <typename T>
struct is_vector_kokkos<T,
	 core::meta::enable_if_t<
	   // kokkos vector is it is a view and has rank=1
	   Kokkos::is_view<T>::value && 
	   T::traits::rank==1>
      > : std::true_type{};

//----------------------------------------------------------------------
  
 
//////////////////////
} // namespace meta
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

#define STATIC_ASSERT_IS_VECTOR_EPETRA(TYPE) \
  static_assert( core::meta::is_vector_epetra<TYPE>::value, \
		 "THIS_IS_NOT_A_VECTOR_EPETRA")
#define STATIC_ASSERT_IS_NOT_VECTOR_EPETRA(TYPE) \
  static_assert( !core::meta::is_vector_epetra<TYPE>::value, \
		 "THIS_IS_A_VECTOR_EPETRA")

#define STATIC_ASSERT_IS_VECTOR_KOKKOS(TYPE) \
  static_assert( core::meta::is_vector_kokkos<TYPE>::value, \
		 "THIS_IS_NOT_A_VECTOR_KOKKOS")
#define STATIC_ASSERT_IS_NOT_VECTOR_KOKKOS(TYPE) \
  static_assert( !core::meta::is_vector_kokkos<TYPE>::value, \
		 "THIS_IS_A_VECTOR_KOKKOS")
  
  
/////////////////
} // namespace core
/////////////////

#endif
