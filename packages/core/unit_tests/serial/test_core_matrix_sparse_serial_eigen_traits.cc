#include <gtest/gtest.h>
#include "matrix/core_matrix_dense_serial_eigen.hpp"
#include "matrix/core_matrix_dense_serial_stdlib.hpp"
#include "matrix/core_matrix_sparse_serial_eigen.hpp"
#include "matrix/core_matrix_meta.hpp"


template<typename scalar, int storage, typename index>
struct typesEig{
  using sc_t = scalar;
  using mat_t = Eigen::SparseMatrix<sc_t,storage,index>;
  using storeindex_t = index;
};

template<typename T>
struct core_matrix_sparse_serial_eigen_traitsTest
  : public ::testing::Test{
public:
  using native_t = typename T::mat_t;
  STATIC_ASSERT_IS_NOT_MATRIX_DENSE_SERIAL_EIGEN(native_t);
  STATIC_ASSERT_IS_MATRIX_SPARSE_SERIAL_EIGEN(native_t);
  using myvec_t = core::matrix<native_t>;
  using vecTrait = core::details::traits<myvec_t>;

  //need this void method otherwise the static assertion for type dont work
  void check(){
    ::testing::StaticAssertTypeEq<typename vecTrait::scalar_t,
				  typename T::sc_t>();
    ::testing::StaticAssertTypeEq<typename vecTrait::ordinal_t,
				  typename T::storeindex_t>();
    ::testing::StaticAssertTypeEq<typename
				  vecTrait::wrapped_t,native_t>();
    ::testing::StaticAssertTypeEq<typename
				  vecTrait::derived_t,myvec_t>();

    EXPECT_EQ(vecTrait::isMatrix, 1);
    EXPECT_EQ(vecTrait::isEigen, 1);
    EXPECT_EQ(vecTrait::isDense, 0);
    EXPECT_EQ(vecTrait::isSparse, 1);
    EXPECT_EQ(vecTrait::isVector, 0);
    EXPECT_EQ(vecTrait::isSerial, 1);
    EXPECT_EQ(vecTrait::isDistributed, 0);
    EXPECT_EQ(vecTrait::isStdlib, 0);
    EXPECT_EQ(vecTrait::isStatic, 0);
  }

  // virtual void SetUp(){}
  // virtual void TearDown(){}
};

typedef ::testing::Types< typesEig<int,Eigen::RowMajor,int>,
			  typesEig<int,Eigen::ColMajor,short>,
			  typesEig<int,Eigen::RowMajor,long long>,
			  //
			  typesEig<double,Eigen::RowMajor,int>,
			  typesEig<double,Eigen::ColMajor,short>,
			  typesEig<double,Eigen::RowMajor,long long>,
			  //
			  typesEig<float,Eigen::RowMajor,int>,
			  typesEig<float,Eigen::ColMajor,short>,
			  typesEig<float,Eigen::RowMajor,long long>,
			  //
			  typesEig<unsigned int,Eigen::RowMajor,int>,
			  typesEig<unsigned int,Eigen::ColMajor,short>,
			  typesEig<unsigned int,Eigen::RowMajor,long long>,
			  //
			  typesEig<std::complex<double>,Eigen::RowMajor,int>,
			  typesEig<std::complex<double>,Eigen::ColMajor,short>,
			  typesEig<std::complex<double>,Eigen::RowMajor,long long>
			  > MyTypes;
TYPED_TEST_CASE(core_matrix_sparse_serial_eigen_traitsTest, MyTypes);

TYPED_TEST(core_matrix_sparse_serial_eigen_traitsTest, traits)
{
  //this runs all types, no need to put anything
  this->check();
}
