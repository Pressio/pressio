
#ifndef CORE_NATIVE_VECTOR_META_VECTOR_META_HPP_
#define CORE_NATIVE_VECTOR_META_VECTOR_META_HPP_

#include "core_meta_basic.hpp"
#include "core_meta_detection_idiom.hpp"
//std::vector 
#include <vector>
// Eigen dynamic and static vectors
#include <Eigen/Dense>

#ifdef HAVE_BLAZE
// blaze offers few things but we include two for now
#include <blaze/math/DynamicVector.h>
#include <blaze/math/StaticVector.h>
#endif

#ifdef HAVE_ARMADILLO
#include <armadillo>
#endif

#ifdef HAVE_TRILINOS
// tpetra
#include <Tpetra_Vector.hpp>
// epetra 
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
// kokkos views
#include <Kokkos_Core.hpp>
#endif


namespace rompp{ namespace core{ namespace meta {

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
//----------------------------------------------
  
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
//--------------------------------------------
  
#ifdef HAVE_TRILINOS  
template <typename T, typename enable = void>
struct is_vector_epetra : std::false_type {};

template <typename T>
struct is_vector_epetra<T,
      typename
      std::enable_if<
	std::is_same<T,Epetra_Vector>::value 
	>::type
      > : std::true_type{};
#endif 
//--------------------------------------------

#ifdef HAVE_TRILINOS  
template <typename T, typename enable = void>
struct is_vector_tpetra : std::false_type {};

template <typename T>
struct is_vector_tpetra<T,
      typename
      std::enable_if<
	std::is_same<T,
		     Tpetra::Vector<
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
  
#ifdef HAVE_TRILINOS 
template <typename T, typename enable = void>
struct is_vector_kokkos : std::false_type {};

template <typename T>
struct is_vector_kokkos<T,
	 core::meta::enable_if_t<
	   // kokkos vector is it is a view and has rank=1
	   Kokkos::is_view<T>::value && 
	   T::traits::rank==1>
      > : std::true_type{};
#endif
//--------------------------------------------

  
#ifdef HAVE_BLAZE
template <typename T, typename enable = void>
struct is_static_vector_blaze : std::false_type {};

template <typename T>
struct is_static_vector_blaze<T,
	 core::meta::enable_if_t<
	   blaze::IsStatic<typename T::This>::value
	   >
      > : std::true_type{};

template <typename T, typename enable = void>
struct is_dynamic_row_vector_blaze : std::false_type {};

template <typename T>
struct is_dynamic_row_vector_blaze<T,
	 core::meta::enable_if_t<
	   std::is_same<T, blaze::DynamicVector<
				typename T::ElementType,
				blaze::rowVector>
			>::value
	   >
      > : std::true_type{};
  
template <typename T, typename enable = void>
struct is_dynamic_column_vector_blaze : std::false_type {};

template <typename T>
struct is_dynamic_column_vector_blaze<T,
	 core::meta::enable_if_t<
	   std::is_same<T, blaze::DynamicVector<
				typename T::ElementType,
				blaze::columnVector>
			>::value
	   >
      > : std::true_type{};

template <typename T, typename enable = void>
struct is_dynamic_vector_blaze : std::false_type {};

template <typename T>
struct is_dynamic_vector_blaze<T,
	   core::meta::enable_if_t<
	     is_dynamic_row_vector_blaze<T>::value ||
	     is_dynamic_column_vector_blaze<T>::value
	   >
      > : std::true_type{};

#endif
//--------------------------------------------


#ifdef HAVE_ARMADILLO
template <typename T, typename enable = void>
struct is_armadillo_column_vector : std::false_type {};

template <typename T>
struct is_armadillo_column_vector<T,
	 core::meta::enable_if_t<
	   std::is_same<T,
     	    arma::Col<typename T::elem_type>
			>::value
	   >
      > : std::true_type{};

template <typename T, typename enable = void>
struct is_armadillo_row_vector : std::false_type {};

template <typename T>
struct is_armadillo_row_vector<T,
	 core::meta::enable_if_t<
	   std::is_same<T,
     	    arma::Row<typename T::elem_type>
			>::value
	   >
      > : std::true_type{};
  
template <typename T, typename enable = void>
struct is_vector_armadillo : std::false_type {};

template <typename T>
struct is_vector_armadillo<T,
	 core::meta::enable_if_t<
	   is_armadillo_row_vector<T>::value or
	   is_armadillo_column_vector<T>::valu
	   >
      > : std::true_type{};

#endif
//--------------------------------------------



}}}//end namespace rompp::core::meta
#endif
