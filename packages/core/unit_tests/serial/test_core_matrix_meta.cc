
#include <gtest/gtest.h>
#include "matrix/core_matrix_meta.hpp"
#include "vector/core_vector_meta.hpp"

template<typename scalar, int rows, int cols>
struct typesEig{
  using sc_t = scalar;
  constexpr static int nr = rows;
  constexpr static int nc = cols;
};


template <typename T>
class MatrixEigenMetaTest : public ::testing::Test {
 public:
  using native_t = Eigen::Matrix<typename T::sc_t,T::nr, T::nc>;

  //this is not even necessary, but lets put it for consistency with other tests
  void check(){
   STATIC_ASSERT_IS_MATRIX_DENSE_SERIAL_EIGEN(native_t);
   STATIC_ASSERT_IS_NOT_MATRIX_DENSE_SERIAL_STDLIB(native_t);
   STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(native_t);
   STATIC_ASSERT_IS_NOT_VECTOR_STDLIB(native_t);
   STATIC_ASSERT_IS_NOT_VECTOR_EPETRA(native_t);
  }
};

typedef ::testing::Types< typesEig<int,-1,-1>, typesEig<int,2,2>,
			  typesEig<int,-1,4>, typesEig<int,4,-1>,
			  typesEig<int,11,11>, typesEig<int,11,11>,
			  //
			  typesEig<double,-1,-1>, typesEig<double,2,2>,
			  typesEig<double,-1,4>, typesEig<double,4,-1>,
			  typesEig<double,11,11>, typesEig<double,11,11>,
			  //
			  typesEig<float,-1,-1>, typesEig<float,2,2>,
			  typesEig<float,-1,4>, typesEig<float,4,-1>,
			  typesEig<float,11,11>, typesEig<float,11,11>,
			  //
			  typesEig<std::complex<double>,-1,-1>,
			  typesEig<std::complex<double>,2,2>,
			  typesEig<std::complex<double>,-1,4>,
			  typesEig<std::complex<double>,4,-1>,
			  typesEig<std::complex<double>,11,11>,
			  typesEig<std::complex<double>,11,11>,
			  //
			  typesEig<unsigned int,-1,-1>, typesEig<unsigned int,2,2>,
			  typesEig<unsigned int,-1,4>, typesEig<unsigned int,4,-1>,
			  typesEig<unsigned int,11,11>, typesEig<unsigned int,11,11>
			  > MyTypes;
TYPED_TEST_CASE(MatrixEigenMetaTest, MyTypes);

TYPED_TEST(MatrixEigenMetaTest, metaEigen)
{
  //this runs all types, no need to put anything else
  this->check();
}



template <typename T>
class MatrixStdlibTest{
 public:
  using native_t = std::vector<std::vector<T>>; 
  STATIC_ASSERT_IS_NOT_MATRIX_DENSE_SERIAL_EIGEN(native_t);
  STATIC_ASSERT_IS_MATRIX_DENSE_SERIAL_STDLIB(native_t);
  STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(native_t);
  STATIC_ASSERT_IS_NOT_VECTOR_STDLIB(native_t);
  STATIC_ASSERT_IS_NOT_VECTOR_EPETRA(native_t);
};

TEST(core_matrix_meta, metaStdlib)
{
  MatrixStdlibTest<double> a;
  MatrixStdlibTest<int> a1;
  MatrixStdlibTest<long int> a2;
  MatrixStdlibTest<unsigned> a3;
  MatrixStdlibTest<std::complex<double>> a4; 
}

TEST(core_matrix_meta, metaSparseEigenSerial)
{
  using native_t = Eigen::SparseMatrix<double>;
  STATIC_ASSERT_IS_MATRIX_SPARSE_SERIAL_EIGEN(native_t); 
  STATIC_ASSERT_IS_NOT_MATRIX_DENSE_SERIAL_EIGEN(native_t); 
  STATIC_ASSERT_IS_NOT_MATRIX_DENSE_SERIAL_STDLIB(native_t); 

  STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(native_t);
  STATIC_ASSERT_IS_NOT_VECTOR_STDLIB(native_t);
  STATIC_ASSERT_IS_NOT_VECTOR_EPETRA(native_t);  
}

