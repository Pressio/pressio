
#ifndef CORE_MATRIX_MATRIX_META_HPP_
#define CORE_MATRIX_MATRIX_META_HPP_

#include "../../meta/core_meta_basic.hpp"
#include "../../vector/meta/core_vector_meta.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <Epetra_CrsMatrix.h>


namespace core{
namespace meta {

template <typename T, typename enable = void>
struct is_matrix_dense_serial_eigen : std::false_type {};

/* a type is a dense eigen matrix if
   (a) it is not a eigen vector, 
   (b) matches a eigen matrix sized at compile time
   (c) its rows and cols are fully dynamic
*/
template<typename T>
struct is_matrix_dense_serial_eigen<T,
    typename
    std::enable_if< !is_vector_eigen<T>::value &&
	(std::is_same<T,Eigen::Matrix<typename T::Scalar,
	 T::RowsAtCompileTime,
	 T::ColsAtCompileTime
	 >
	 >::value ||
	 std::is_same<T, Eigen::Matrix<typename T::Scalar,
	 Eigen::Dynamic,
	 Eigen::Dynamic
	 >
	 >::value
	 )
	>::type
		   > : std::true_type{};

//----------------------------------------------------------------------
  
template <typename T, typename enable = void>
struct is_matrix_sparse_serial_eigen : std::false_type {};

/* a type is a sparse eigen matrix if
   (a) it is not a eigen vector, 
   (b) matches a eigen sparse matrix sized at compile time
   (c) it is not a dense matrix
*/
template<typename T>
struct is_matrix_sparse_serial_eigen<T, typename
  std::enable_if< !is_vector_eigen<T>::value &&
		  !is_matrix_dense_serial_eigen<T>::value &&
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
struct is_matrix_dense_serial_stdlib : std::false_type {};

template <typename T>
struct is_matrix_dense_serial_stdlib<T,
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
struct sparse_serial_eigen_same_storage : std::false_type{};

template <typename T1, typename T2>
struct sparse_serial_eigen_same_storage<
  T1, T2, typename
  std::enable_if<
	    (T1::isRowMajor && T2::isRowMajor) ||
	    (T1::isColMajor && T2::isColMajor)
	    >::type
  > : std::true_type{};
  
//----------------------------------------------------------------------


template <typename T, typename enable = void>
struct is_matrix_sparse_distributed_epetra
  : std::false_type {};

template<typename T>
struct is_matrix_sparse_distributed_epetra<T,
    typename std::enable_if<
      std::is_same<T, Epetra_CrsMatrix>::value
      >::type >
  : std::true_type{};


  
/////////////////////
} // namespace meta
/////////////////////

  
#define STATIC_ASSERT_IS_MATRIX_DENSE_SERIAL_EIGEN(TYPE) \
  static_assert( core::meta::is_matrix_dense_serial_eigen<TYPE>::value,	\
		 "THIS_IS_NOT_A_MATRIX_DENSE_SERIAL_EIGEN")
#define STATIC_ASSERT_IS_NOT_MATRIX_DENSE_SERIAL_EIGEN(TYPE) \
  static_assert( !core::meta::is_matrix_dense_serial_eigen<TYPE>::value, \
		 "THIS_IS_A_MATRIX_DENSE_SERIAL_EIGEN")

  
#define STATIC_ASSERT_IS_MATRIX_SPARSE_SERIAL_EIGEN(TYPE) \
  static_assert( core::meta::is_matrix_sparse_serial_eigen<TYPE>::value,	\
		 "THIS_IS_NOT_A_MATRIX_SPARSE_SERIAL_EIGEN")
#define STATIC_ASSERT_IS_NOT_MATRIX_SPARSE_SERIAL_EIGEN(TYPE) \
  static_assert( !core::meta::is_matrix_sparse_serial_eigen<TYPE>::value, \
		 "THIS_IS_A_MATRIX_SPARSE_SERIAL_EIGEN")

  
#define STATIC_ASSERT_IS_MATRIX_DENSE_SERIAL_STDLIB(TYPE)	      \
  static_assert( core::meta::is_matrix_dense_serial_stdlib<TYPE>::value, \
		 "THIS_IS_NOT_A_MATRIX_DENSE_SERIAL_STDLIB")
#define STATIC_ASSERT_IS_NOT_MATRIX_DENSE_SERIAL_STDLIB(TYPE) \
  static_assert( !core::meta::is_matrix_dense_serial_stdlib<TYPE>::value, \
		 "THIS_IS_A_MATRIX_DENSE_SERIAL_STDLIB")

  
#define STATIC_ASSERT_IS_MATRIX_SPARSE_DISTRIBUTED_EPETRA(TYPE)	      \
  static_assert( core::meta::is_matrix_sparse_distributed_epetra<TYPE>::value, \
		 "THIS_IS_NOT_A_MATRIX_SPARSE_DIST_EPETRA")
#define STATIC_ASSERT_IS_NOT_MATRIX_SPARSE_DISTRIBUTED_EPETRA(TYPE) \
  static_assert( !core::meta::is_matrix_sparse_distributed_epetra<TYPE>::value, \
		 "THIS_IS_A_MATRIX_SPARSE_DIST_EPETRA")
  
} // namespace core
#endif
