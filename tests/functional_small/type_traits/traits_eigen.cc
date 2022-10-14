#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"
#include "traits_shared.hpp"

namespace pressio { namespace traits { namespace test {

template <typename T, int rank, bool is_dynamic>
void test_eigen_container_traits()
{
  test_container_traits<
    T,
    rank,
    typename T::Scalar
  >();

  // negative
  test_is_not_teuchos_container<T>();
  test_is_not_epetra_container<T>();
  test_is_not_tpetra_container<T>();
  test_is_not_kokkos_container<T>();
}

//*******************************
// Eigen vector
//*******************************

/*
  Verify values of Eigen vector traits and relatad predicates
*/
template <
  typename T,
  bool is_dynamic,
  bool is_row_vector,
  typename traits = pressio::Traits<T>
>
void test_eigen_vector_type_traits()
{
  // traits
  test_eigen_container_traits<T, 1, is_dynamic>();
  // vector predicates
  static_assert(pressio::is_vector_eigen<T>::value, "");
  static_assert(pressio::is_dynamic_vector_eigen<T>::value == is_dynamic, "");
  static_assert(pressio::is_static_vector_eigen<T>::value == !is_dynamic, "");
  static_assert(pressio::is_dynamic_row_vector_eigen<T>::value == (is_dynamic && is_row_vector), "");
  static_assert(pressio::is_static_row_vector_eigen<T>::value  == (!is_dynamic && is_row_vector), "");
  static_assert(pressio::is_dynamic_column_vector_eigen<T>::value == (is_dynamic && !is_row_vector), "");
  static_assert(pressio::is_static_column_vector_eigen<T>::value == (!is_dynamic && !is_row_vector), "");
}

#define TEST_EIGEN_VECTOR(Type, is_dynamic, is_row_vector) \
TEST(type_traits, Type) { \
  test_eigen_vector_type_traits<Type, is_dynamic, is_row_vector>(); \
}

using eigen_vector_dynamic_row = Eigen::Matrix<double, 1, Eigen::Dynamic>;
using eigen_vector_dynamic_col = Eigen::Matrix<float, Eigen::Dynamic, 1>;
using eigen_vector_static_row  = Eigen::Matrix<double, 1, 32>;
using eigen_vector_static_col  = Eigen::Matrix<float, 32, 1>;

//                Type(Name)                dyamic row_vector
TEST_EIGEN_VECTOR(eigen_vector_dynamic_row, true,  true)
TEST_EIGEN_VECTOR(eigen_vector_dynamic_col, true,  false)
TEST_EIGEN_VECTOR(eigen_vector_static_row,  false, true)
TEST_EIGEN_VECTOR(eigen_vector_static_col,  false, false)

//*******************************
// Eigen dense matrix
//*******************************

/*
  Verify values of Eigen matrix traits and relatad predicates
*/
template <
  typename T,
  bool is_dynamic,
  typename traits = pressio::Traits<T>
>
void test_eigen_matrix_type_traits()
{
  // traits
  test_eigen_container_traits<T, 2, is_dynamic>();

  constexpr bool row_major = T::IsRowMajor == 1;
  // static_assert(traits::is_row_major == row_major);
  // static_assert(traits::is_col_major == !row_major);

  // matrix predicates
  static_assert(pressio::is_dense_matrix_eigen<T>::value, "");
  static_assert(!pressio::is_sparse_matrix_eigen<T>::value, "");
  static_assert(pressio::is_static_dense_matrix_eigen<T>::value == !is_dynamic, "");
  static_assert(pressio::is_dynamic_dense_matrix_eigen<T>::value == is_dynamic, "");
  static_assert(pressio::is_dense_row_major_matrix_eigen<T>::value == row_major, "");
}

#define TEST_EIGEN_MATRIX(Type, is_dynamic) \
TEST(type_traits, Type) { \
  test_eigen_matrix_type_traits<Type, is_dynamic>(); \
}

using eigen_dense_matrix_dynamic_rowmajor = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using eigen_dense_matrix_static_colmajor  = Eigen::Matrix<float, 32, 32, Eigen::ColMajor>;

//                Type(Name)                           dyamic
TEST_EIGEN_MATRIX(eigen_dense_matrix_dynamic_rowmajor, true)
TEST_EIGEN_MATRIX(eigen_dense_matrix_static_colmajor,  false)

//*******************************
// Eigen sparse matrix
//*******************************

/*
  Verify values of Eigen matrix traits and relatad predicates
*/
template <
  typename T,
  typename traits = pressio::Traits<T>
>
void test_eigen_sparse_matrix_type_traits()
{
  // traits
  test_eigen_container_traits<T, 2, true>();

  // sparse matrix predicates
  static_assert(pressio::is_sparse_matrix_eigen<T>::value, "");
  static_assert(!pressio::is_dense_matrix_eigen<T>::value, "");
  // static_assert(pressio::sparse_sharedmem_eigen_same_storage<T, T>::value); // NOT USED
}

#define TEST_EIGEN_SPARSE_MATRIX(Type) \
TEST(type_traits, Type) { \
  test_eigen_sparse_matrix_type_traits<Type>(); \
}

using eigen_sparse_matrix_dynamic_rowmajor = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using eigen_sparse_matrix_static_colmajor  = Eigen::SparseMatrix<float, Eigen::ColMajor>;

TEST_EIGEN_SPARSE_MATRIX(eigen_sparse_matrix_dynamic_rowmajor)
TEST_EIGEN_SPARSE_MATRIX(eigen_sparse_matrix_static_colmajor)

}}} // pressio::traits::test
