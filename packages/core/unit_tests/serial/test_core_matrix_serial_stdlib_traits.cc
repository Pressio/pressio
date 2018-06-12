#include <gtest/gtest.h>
#include "matrix/core_matrix_dense_serial_eigen.hpp"
#include "matrix/core_matrix_dense_serial_stdlib.hpp"
#include "matrix/core_matrix_meta.hpp"

template<typename T>
struct core_matrix_dense_serial_stdlib_traitsTest
  : public ::testing::Test{
public:
  using native_t = std::vector<std::vector<T>>;
  STATIC_ASSERT_IS_MATRIX_DENSE_SERIAL_EIGEN(native_t);
  using myvec_t = core::matrix<native_t>;
  using vecTrait = core::details::traits<myvec_t>;

  //need to have this void method otherwise the static assertion for type cannot be used
  void check(){
    ::testing::StaticAssertTypeEq<typename
				  vecTrait::scalar_t, typename T::sc_t>();
    ::testing::StaticAssertTypeEq<typename
				  vecTrait::ordinal_t,int>();
    ::testing::StaticAssertTypeEq<typename
				  vecTrait::wrapped_t,native_t>();
    ::testing::StaticAssertTypeEq<typename
				  vecTrait::derived_t,myvec_t>();

    EXPECT_EQ(vecTrait::isMatrix, 1);
    EXPECT_EQ(vecTrait::isEigen, 0);
    EXPECT_EQ(vecTrait::isDense, 1);
    EXPECT_EQ(vecTrait::isVector, 0);
    EXPECT_EQ(vecTrait::isSerial, 1);
    EXPECT_EQ(vecTrait::isDistributed, 0);
    EXPECT_EQ(vecTrait::isStatic, 0);
  }

  // virtual void SetUp(){}
  // virtual void TearDown(){}
};

typedef ::testing::Types<double, int, unsigned int,
			 float, long long, std::comlex<double>
			 > MyTypes;
TYPED_TEST_CASE(core_matrix_dense_serial_stdlib_traitsTest, MyTypes);

TYPED_TEST(core_matrix_dense_serial_stdlib_traitsTest, traits)
{
  //this runs all types, no need to put anything
  this->check();
}
