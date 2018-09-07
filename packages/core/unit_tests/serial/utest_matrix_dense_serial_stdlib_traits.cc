

#include <gtest/gtest.h>
#include "CORE_MATRIX"

template<typename T>
struct core_matrix_dense_serial_stdlib_traitsTest
  : public ::testing::Test{
public:
  using native_t = std::vector<std::vector<T>>;
  STATIC_ASSERT_IS_NOT_MATRIX_DENSE_SHAREDMEM_EIGEN(native_t);
  STATIC_ASSERT_IS_MATRIX_DENSE_SHAREDMEM_STDLIB(native_t);
  using myvec_t = core::Matrix<native_t>;
  using vecTrait = core::details::traits<myvec_t>;

  //need to have this void method otherwise the static assertion for type cannot be used
  void check(){
    ::testing::StaticAssertTypeEq<typename
				  vecTrait::scalar_t,T>();
    ::testing::StaticAssertTypeEq<typename
				  vecTrait::ordinal_t,int>();
    ::testing::StaticAssertTypeEq<typename
				  vecTrait::wrapped_t,native_t>();
    ::testing::StaticAssertTypeEq<typename
				  vecTrait::derived_t,myvec_t>();

    ASSERT_TRUE(vecTrait::isMatrix == 1);
    ASSERT_TRUE(vecTrait::isEigen == 0);
    ASSERT_TRUE(vecTrait::isDense == 1);
    ASSERT_TRUE(vecTrait::isVector == 0);
    ASSERT_TRUE(vecTrait::isSharedMem == 1);
    ASSERT_TRUE(vecTrait::isStatic == 0);
  }

  // virtual void SetUp(){}
  // virtual void TearDown(){}
};

typedef ::testing::Types<double, int, unsigned int,
			 float, long long, std::complex<double>
			 > MyTypes;
TYPED_TEST_CASE(core_matrix_dense_serial_stdlib_traitsTest, MyTypes);

TYPED_TEST(core_matrix_dense_serial_stdlib_traitsTest, traits)
{
  //this runs all types, no need to put anything
  this->check();
}
