

#include <gtest/gtest.h>
#include "CONTAINERS_MATRIX"

template<typename T>
struct containers_matrix_dense_serial_stdlib_traitsTest
  : public ::testing::Test{
public:
  using native_t = std::vector<std::vector<T>>;
  STATIC_ASSERT_IS_NOT_MATRIX_DENSE_SHAREDMEM_EIGEN(native_t);
  STATIC_ASSERT_IS_MATRIX_DENSE_SHAREDMEM_STDLIB(native_t);
  using myvec_t = containers::Matrix<native_t>;
  using vecTrait = containers::details::traits<myvec_t>;

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

    ASSERT_TRUE(vecTrait::is_matrix == 1);
    ASSERT_TRUE(vecTrait::wrapped_matrix_identifier 
      == containers::details::WrappedMatrixIdentifier::CppStdLib);

    ASSERT_TRUE(vecTrait::is_dense == 1);
    ASSERT_TRUE(vecTrait::is_vector == 0);
    ASSERT_TRUE(vecTrait::is_shared_mem == 1);
    ASSERT_TRUE(vecTrait::is_static == 0);
  }

  // virtual void SetUp(){}
  // virtual void TearDown(){}
};

typedef ::testing::Types<double, int, unsigned int,
			 float, long long, std::complex<double>
			 > MyTypes;
TYPED_TEST_CASE(containers_matrix_dense_serial_stdlib_traitsTest, MyTypes);

TYPED_TEST(containers_matrix_dense_serial_stdlib_traitsTest, traits)
{
  //this runs all types, no need to put anything
  this->check();
}
