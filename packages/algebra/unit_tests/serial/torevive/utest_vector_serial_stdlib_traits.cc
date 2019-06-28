
#include <gtest/gtest.h>
#include "ALGEBRA_VECTOR"

template <typename T>
struct StdVecChecker
{
public:
  using stdV_t = std::vector<T>;
  STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(stdV_t);
  STATIC_ASSERT_IS_VECTOR_STDLIB(stdV_t);

  using myvec_t = algebra::Vector<stdV_t>;
  STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(myvec_t);
  STATIC_ASSERT_IS_NOT_VECTOR_STDLIB(myvec_t);

  using vecTrait = algebra::details::traits<myvec_t>;
  void check(){   
    ::testing::StaticAssertTypeEq<typename
				  vecTrait::scalar_t,T>();
    ::testing::StaticAssertTypeEq<typename
				  vecTrait::ordinal_t,
				  algebra::default_types::local_ordinal_t>();
    ::testing::StaticAssertTypeEq<typename
				  vecTrait::wrapped_t,stdV_t>();
    ::testing::StaticAssertTypeEq<typename
				  vecTrait::derived_t,myvec_t>();
    ASSERT_TRUE(vecTrait::is_vector == 1);
    ASSERT_TRUE(vecTrait::wrapped_vector_identifier 
      == algebra::details::WrappedVectorIdentifier::CppStdLib);
    ASSERT_TRUE(vecTrait::is_shared_mem == 1);
  }
};

struct algebra_vector_serial_stdlib_traits_Fixture : public ::testing::Test{
public:
  // virtual void SetUp(){}
  // virtual void TearDown(){}
  StdVecChecker<int> a;
  StdVecChecker<double> b;
  StdVecChecker<float> c;
  StdVecChecker<std::complex<double>> d;
};

TEST_F(algebra_vector_serial_stdlib_traits_Fixture, traits)
{
  a.check();
  b.check();
  c.check();
  d.check(); 
 
  // check that a matrix is not a vector
  using stdmat_t = std::vector<std::vector<double>>;
  STATIC_ASSERT_IS_NOT_VECTOR_STDLIB(stdmat_t);
  using stdmat_t2 = std::vector<std::vector<int>>;
  STATIC_ASSERT_IS_NOT_VECTOR_STDLIB(stdmat_t2);
}
