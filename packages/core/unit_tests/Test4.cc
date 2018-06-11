
#include <gtest/gtest.h>
#include "vector/core_vector_serial_eigen.hpp"
#include "vector/core_vector_serial_stdlib.hpp"
#include "vector/core_vector_serial_userdefined.hpp"
#include "vector/core_vector_meta.hpp"
#include "core_static_assert_definitions.hpp"


template <typename T>
struct StdVecChecker
{
public:
  using stdV_t = std::vector<T>;
  STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(stdV_t);
  STATIC_ASSERT_IS_STDLIB_VECTOR(stdV_t);

  using myvec_t = core::vector<stdV_t>;
  STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(myvec_t);
  STATIC_ASSERT_IS_NOT_STDLIB_VECTOR(myvec_t);

  using vecTrait = core::details::traits<myvec_t>;
  void check(){   
    ::testing::StaticAssertTypeEq<typename
				  vecTrait::scalar_t,T>();
    ::testing::StaticAssertTypeEq<typename
				  vecTrait::ordinal_t,
				  core::defaultTypes::local_ordinal_t>();
    ::testing::StaticAssertTypeEq<typename
				  vecTrait::wrapped_t,stdV_t>();
    ::testing::StaticAssertTypeEq<typename
				  vecTrait::derived_t,myvec_t>();
    EXPECT_EQ(vecTrait::isVector, 1);
    EXPECT_EQ(vecTrait::isEigen, 0);
    EXPECT_EQ(vecTrait::isSerial, 1);
    EXPECT_EQ(vecTrait::isSTDVector, 1);
    EXPECT_EQ(vecTrait::isDistributed, 0);
  }
};

struct StdVecFixture : public ::testing::Test{
public:
  // virtual void SetUp(){}
  // virtual void TearDown(){}
  StdVecChecker<int> a;
  StdVecChecker<double> b;
  StdVecChecker<float> c;
  StdVecChecker<std::complex<double>> d;
};

TEST_F(StdVecFixture, StdVectorTraits)
{
  a.check();
  b.check();
  c.check();
  d.check(); 
 
  // check that a matrix is not a vector
  using stdmat_t = std::vector<std::vector<double>>;
  STATIC_ASSERT_IS_NOT_STDLIB_VECTOR(stdmat_t);
  using stdmat_t2 = std::vector<std::vector<int>>;
  STATIC_ASSERT_IS_NOT_STDLIB_VECTOR(stdmat_t2);
}
