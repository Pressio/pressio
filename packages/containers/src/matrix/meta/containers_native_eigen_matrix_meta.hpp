
#ifndef CONTAINERS_NATIVE_EIGEN_MATRIX_META_HPP_
#define CONTAINERS_NATIVE_EIGEN_MATRIX_META_HPP_

#include "../../meta/containers_meta_basic.hpp"
#include "../../vector/meta/containers_native_eigen_vector_meta.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_dense_static_matrix_eigen : std::false_type {};

/* T is a dense STATIC eigen matrix if
 * T is not an eigen vector
 * rows and cols are not = Eigen:Dynamic
 */
template<typename T>
struct is_dense_static_matrix_eigen<
  T,
  typename
  std::enable_if<
    !is_vector_eigen<T>::value and
    std::is_same<
      T,
      Eigen::Matrix<typename T::Scalar,
		    T::RowsAtCompileTime,
		    T::ColsAtCompileTime
		    >
      >::value and
    T::RowsAtCompileTime != Eigen::Dynamic and
    T::ColsAtCompileTime != Eigen::Dynamic
    >::type
  > : std::true_type{};
//----------------------------------------------------------------------


template <typename T, typename enable = void>
struct is_dense_dynamic_matrix_eigen : std::false_type {};

/* T is a dense DYNAMIC eigen matrix if
 * is not an eigen vector
 * is not a static dense matrix
 * both # rows and cols are dynamic
*/
template<typename T>
struct is_dense_dynamic_matrix_eigen<
  T,
  typename
  std::enable_if<
    !is_vector_eigen<T>::value and
    !is_dense_static_matrix_eigen<T>::value and
    std::is_same<
      T, Eigen::Matrix<typename T::Scalar,
		       Eigen::Dynamic, Eigen::Dynamic
		       >
      >::value
    >::type
  > : std::true_type{};


/* T is a dense DYNAMIC eigen matrix if
 * is not an eigen vector
 * is not a static dense matrix
 * # of rows is static but cols are dynamic
*/
template<typename T>
struct is_dense_dynamic_matrix_eigen<
  T,
  typename
  std::enable_if<
    !is_vector_eigen<T>::value and
    !is_dense_static_matrix_eigen<T>::value and
    std::is_same<
      T, Eigen::Matrix<typename T::Scalar,
		       T::RowsAtCompileTime,
		       Eigen::Dynamic>
      >::value and
    T::RowsAtCompileTime != Eigen::Dynamic
    >::type
  > : std::true_type{};


/* T is a dense DYNAMIC eigen matrix if
 * is not an eigen vector
 * is not a static dense matrix
 * # of rows is dynamic but cols are static
*/
template<typename T>
struct is_dense_dynamic_matrix_eigen<
  T,
  typename
  std::enable_if<
    !is_vector_eigen<T>::value and
    !is_dense_static_matrix_eigen<T>::value and
    std::is_same<
      T, Eigen::Matrix<typename T::Scalar,
		       Eigen::Dynamic,
		       T::ColsAtCompileTime>
      >::value and
    T::ColsAtCompileTime != Eigen::Dynamic
    >::type
  > : std::true_type{};

//----------------------------------------------------------------------


template <typename T, typename enable = void>
struct is_dense_matrix_eigen : std::false_type {};

template<typename T>
struct is_dense_matrix_eigen<
  T, typename
  std::enable_if<
       is_dense_static_matrix_eigen<T>::value or
       is_dense_dynamic_matrix_eigen<T>::value
    >::type
  > : std::true_type{};


//----------------------------------------------------------------------

template <typename T, typename enable = void>
struct is_sparse_matrix_eigen : std::false_type {};

/*
 * T is an eigen sparse matrix if is
 * not an eigen vector
 * is same type as sparse matrix
*/
template<typename T>
struct is_sparse_matrix_eigen<T, typename
  std::enable_if< !is_vector_eigen<T>::value and
		  /*!is_dense_matrix_eigen<T>::value &&*/
		  std::is_same<T,
			       Eigen::SparseMatrix<
			       typename T::Scalar,
				 T::Options,
			       typename T::StorageIndex>
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

}}}//end namespace pressio::containers::meta
#endif
