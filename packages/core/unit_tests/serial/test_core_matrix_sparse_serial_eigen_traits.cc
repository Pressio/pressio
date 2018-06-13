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
  static constexpr int storemode = storage;
};

template<typename T>
struct core_matrix_sparse_serial_eigen_traitsTest
  : public ::testing::Test{
public:
  using native_t = typename T::mat_t;
  STATIC_ASSERT_IS_NOT_MATRIX_DENSE_SERIAL_EIGEN(native_t);
  STATIC_ASSERT_IS_MATRIX_SPARSE_SERIAL_EIGEN(native_t);
  using my_t = core::matrix<native_t>;
  using myTrait = core::details::traits<my_t>;

  //need this void method otherwise the static assertion for type dont work
  void check(){
    ::testing::StaticAssertTypeEq<typename myTrait::scalar_t,
				  typename T::sc_t>();
    ::testing::StaticAssertTypeEq<typename myTrait::ordinal_t,
				  typename T::storeindex_t>();
    ::testing::StaticAssertTypeEq<typename
				  myTrait::wrapped_t,native_t>();
    ::testing::StaticAssertTypeEq<typename
				  myTrait::derived_t,my_t>();

    if (T::storemode == Eigen::RowMajor)
      EXPECT_TRUE(myTrait::isRowMajor == Eigen::StorageOptions::RowMajor);
    else
       EXPECT_TRUE(myTrait::isRowMajor == Eigen::ColMajor);
    
    ASSERT_TRUE(myTrait::isColMajor == !myTrait::isRowMajor);
    ASSERT_TRUE(myTrait::isEigen == 1);
    ASSERT_TRUE(myTrait::isDense == 0);
    ASSERT_TRUE(myTrait::isSparse == 1);
    ASSERT_TRUE(myTrait::isVector == 0);
    ASSERT_TRUE(myTrait::isSerial == 1);
    ASSERT_TRUE(myTrait::isDistributed == 0);
    ASSERT_TRUE(myTrait::isStdlib == 0);
    ASSERT_TRUE(myTrait::isStatic == 0);
  }
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

  using type1 = core::matrix<Eigen::SparseMatrix<double,Eigen::RowMajor,int>>;  
  using type2 = core::matrix<Eigen::SparseMatrix<double,Eigen::ColMajor,int>>;
  static_assert( core::meta::sparseSerialEigenSameStorage<type1,type2>::value == false, "ohh" );
}
