
#include <gtest/gtest.h>
#include "CORE_MATRIX"

//Eigen::MatrixNt
// N can be any one of 2, 3, 4, or X (meaning Dynamic).
// t = i (int), f (float), d (double),
// cf ( complex<float>), or cd (complex<double>)

template<typename scalar, int rows, int cols>
struct typesEig{
  using sc_t = scalar;
  constexpr static int nr = rows;
  constexpr static int nc = cols;
};

template<typename T>
struct core_matrix_dense_serial_eigen_traitsTest
  : public ::testing::Test{
public:
  using native_t = Eigen::Matrix<typename T::sc_t,T::nr, T::nc>;
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

    ASSERT_TRUE(vecTrait::isMatrix == 1);
    ASSERT_TRUE(vecTrait::isEigen == 1);
    ASSERT_TRUE(vecTrait::isDense == 1);
    ASSERT_TRUE(vecTrait::isVector == 0);
    ASSERT_TRUE(vecTrait::isSerial == 1);
    ASSERT_TRUE(vecTrait::isDistributed == 0);
    ASSERT_TRUE(vecTrait::isStatic == (T::nr!=-1 && T::nc!=-1) );
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
TYPED_TEST_CASE(core_matrix_dense_serial_eigen_traitsTest, MyTypes);

TYPED_TEST(core_matrix_dense_serial_eigen_traitsTest, traits)
{
  //this runs all types, no need to put anything
  this->check();
}
