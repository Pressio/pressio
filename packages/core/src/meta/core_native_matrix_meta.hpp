
#ifndef CORE_NATIVE_MATRIX_MATRIX_META_HPP_
#define CORE_NATIVE_MATRIX_MATRIX_META_HPP_

#include "core_meta_basic.hpp"
#include "core_native_vector_meta.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"

#ifdef HAVE_TRILINOS
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_LocalMap.h>
#include <Tpetra_CrsMatrix_decl.hpp>
#endif

namespace rompp{ namespace core{ namespace meta {


template <typename T, typename enable = void>
struct is_matrix_dense_sharedmem_eigen_dynamic : std::false_type {};

/* a type is a dense eigen matrix if
   (a) it is not a eigen vector,
   (b) its rows and cols are fully dynamic*/
template<typename T>
struct is_matrix_dense_sharedmem_eigen_dynamic<
  T,
  typename
  std::enable_if<
    !is_vector_eigen<T>::value &&
    std::is_same<
      T,
      Eigen::Matrix<typename T::Scalar,Eigen::Dynamic, Eigen::Dynamic>
      >::value
    >::type
  > : std::true_type{};

//----------------------------------------------------------------------

template <typename T, typename enable = void>
struct is_matrix_dense_sharedmem_eigen_static : std::false_type {};

/* a type is a dense eigen matrix if
   (a) it is not a eigen vector,
   (b) matches a eigen matrix sized at compile time*/
template<typename T>
struct is_matrix_dense_sharedmem_eigen_static<
  T,
  typename
  std::enable_if<
    !is_vector_eigen<T>::value and
    !is_matrix_dense_sharedmem_eigen_dynamic<T>::value and
    std::is_same<
      T,
      Eigen::Matrix<typename T::Scalar,
		    T::RowsAtCompileTime,
		    T::ColsAtCompileTime
		    >
      >::value
    >::type
  > : std::true_type{};
//----------------------------------------------------------------------


template <typename T, typename enable = void>
struct is_matrix_dense_sharedmem_eigen : std::false_type {};

template<typename T>
struct is_matrix_dense_sharedmem_eigen<
  T, typename
  std::enable_if<
       is_matrix_dense_sharedmem_eigen_static<T>::value ||
       is_matrix_dense_sharedmem_eigen_dynamic<T>::value
    >::type
  > : std::true_type{};


//----------------------------------------------------------------------

template <typename T, typename enable = void>
struct is_matrix_sparse_sharedmem_eigen : std::false_type {};

/* a type is a sparse eigen matrix if
   (a) it is not a eigen vector,
   (b) matches a eigen sparse matrix sized at compile time
   (c) it is not a dense matrix
*/
template<typename T>
struct is_matrix_sparse_sharedmem_eigen<T, typename
  std::enable_if< !is_vector_eigen<T>::value &&
		  !is_matrix_dense_sharedmem_eigen<T>::value &&
		  std::is_same<T,
			       Eigen::SparseMatrix<
			       typename T::Scalar,
				 T::Options,
			       typename T::StorageIndex>
			      >::value
		  >::type
      > : std::true_type{};

//----------------------------------------------------------------------

template <typename T, typename enable = void>
struct is_matrix_dense_sharedmem_stdlib : std::false_type {};

template <typename T>
struct is_matrix_dense_sharedmem_stdlib<T,
     typename
     std::enable_if<
       std::is_same<T,std::vector<
			std::vector<typename
			T::value_type::value_type
			>
		     >
		    >::value
		    >::type
    > : std::true_type{};

//----------------------------------------------------------------------


template <typename T1, typename T2, typename enable = void>
struct sparse_sharedmem_eigen_same_storage : std::false_type{};

template <typename T1, typename T2>
struct sparse_sharedmem_eigen_same_storage<
  T1, T2, typename
  std::enable_if<
	    (T1::is_row_major && T2::is_row_major) ||
	    (T1::is_col_major && T2::is_col_major)
	    >::type
  > : std::true_type{};

//----------------------------------------------------------------------


#ifdef HAVE_TRILINOS
template <typename T, typename enable = void>
struct is_matrix_sparse_distributed_tpetra : std::false_type {};

template <typename T>
struct is_matrix_sparse_distributed_tpetra<T,
      typename
      std::enable_if<
  std::is_same<T,
         Tpetra::CrsMatrix<
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
struct is_matrix_sparse_distributed_epetra
  : std::false_type {};

template<typename T>
struct is_matrix_sparse_distributed_epetra<T,
    typename std::enable_if<
      std::is_same<T, Epetra_CrsMatrix>::value
      >::type >
  : std::true_type{};

//----------------------------------------------------------------------


template <typename T, typename enable = void>
struct is_matrix_dense_distributed_epetra
  : std::false_type {};

template<typename T>
struct is_matrix_dense_distributed_epetra<T,
    typename std::enable_if<
      std::is_same<T, Epetra_MultiVector>::value
      >::type >
  : std::true_type{};
#endif


}}}//end namespace rompp::core::meta
#endif
