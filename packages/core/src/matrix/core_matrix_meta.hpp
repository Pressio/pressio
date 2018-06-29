
#ifndef CORE_MATRIX_META_HPP_
#define CORE_MATRIX_META_HPP_

#include "../meta/core_meta_basic.hpp"
#include "../vector/core_vector_meta.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"

namespace core{
namespace meta {


template <typename T, typename enable = void>
struct is_matrixDenseSerialEigen : std::false_type {};

/* a type is a dense eigen matrix if
   (a) it is not a eigen vector, 
   (b) matches a eigen matrix sized at compile time
   (c) its rows and cols are fully dynamic
*/
template<typename T>
struct is_matrixDenseSerialEigen<T, typename
			      std::enable_if< !is_vectorEigen<T>::value &&
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
struct is_matrixSparseSerialEigen : std::false_type {};

/* a type is a sparse eigen matrix if
   (a) it is not a eigen vector, 
   (b) matches a eigen sparse matrix sized at compile time
   (c) it is not a dense matrix
*/
template<typename T>
struct is_matrixSparseSerialEigen<T, typename
			      std::enable_if< !is_vectorEigen<T>::value &&
					      !is_matrixDenseSerialEigen<T>::value &&
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
struct is_matrixDenseSerialStdlib : std::false_type {};

template <typename T>
struct is_matrixDenseSerialStdlib<T,
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
struct sparseSerialEigenSameStorage : std::false_type{};

template <typename T1, typename T2>
struct sparseSerialEigenSameStorage<T1, T2,
				    typename
				    std::enable_if<
				      (T1::isRowMajor && T2::isRowMajor) ||
				      (T1::isColMajor && T2::isColMajor)
				      >::type
				    > : std::true_type{};
  
  
  
} // namespace meta

 
#define STATIC_ASSERT_IS_MATRIX_DENSE_SERIAL_EIGEN(TYPE) \
  static_assert( core::meta::is_matrixDenseSerialEigen<TYPE>::value,	\
		 "THIS_IS_NOT_A_MATRIX_DENSE_SERIAL_EIGEN")
#define STATIC_ASSERT_IS_NOT_MATRIX_DENSE_SERIAL_EIGEN(TYPE) \
  static_assert( !core::meta::is_matrixDenseSerialEigen<TYPE>::value, \
		 "THIS_IS_A_MATRIX_DENSE_SERIAL_EIGEN")

  
#define STATIC_ASSERT_IS_MATRIX_SPARSE_SERIAL_EIGEN(TYPE) \
  static_assert( core::meta::is_matrixSparseSerialEigen<TYPE>::value,	\
		 "THIS_IS_NOT_A_MATRIX_SPARSE_SERIAL_EIGEN")
#define STATIC_ASSERT_IS_NOT_MATRIX_SPARSE_SERIAL_EIGEN(TYPE) \
  static_assert( !core::meta::is_matrixSparseSerialEigen<TYPE>::value, \
		 "THIS_IS_A_MATRIX_SPARSE_SERIAL_EIGEN")

  
#define STATIC_ASSERT_IS_MATRIX_DENSE_SERIAL_STDLIB(TYPE)	      \
  static_assert( core::meta::is_matrixDenseSerialStdlib<TYPE>::value, \
		 "THIS_IS_NOT_A_MATRIX_DENSE_SERIAL_STDLIB")
#define STATIC_ASSERT_IS_NOT_MATRIX_DENSE_SERIAL_STDLIB(TYPE) \
  static_assert( !core::meta::is_matrixDenseSerialStdlib<TYPE>::value, \
		 "THIS_IS_A_MATRIX_DENSE_SERIAL_STDLIB")



} // namespace core
#endif
