#include <gtest/gtest.h>
#include "CORE_VECTOR"

//Eigen::MatrixNt
// N can be any one of 2, 3, 4, or X (meaning Dynamic).
// t = i (int), f (float), d (double),
// cf ( complex<float>), or cd (complex<double>)

struct core_vector_serial_eigen_traits_Fixture
  : public ::testing::Test{
public:

  template <typename T, int row, int col>
  struct EigenVecChecker
  {
  public:
    using eigV_t = Eigen::Matrix<T, row, col>;
    STATIC_ASSERT_IS_VECTOR_EIGEN(eigV_t);
    STATIC_ASSERT_IS_NOT_VECTOR_STDLIB(eigV_t);

    using myvec_t = core::Vector<eigV_t>;
    STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(myvec_t);
    STATIC_ASSERT_IS_NOT_VECTOR_STDLIB(myvec_t);
    STATIC_ASSERT_IS_CORE_VECTOR_WRAPPER(myvec_t);

    using vecTrait = core::details::traits<myvec_t>;
    void check(){   

      ::testing::StaticAssertTypeEq<typename
				    vecTrait::scalar_t,T>();
      ::testing::StaticAssertTypeEq<typename
				    vecTrait::ordinal_t,int>();
      ::testing::StaticAssertTypeEq<typename
				    vecTrait::wrapped_t,eigV_t>();
      ::testing::StaticAssertTypeEq<typename
				    vecTrait::derived_t,myvec_t>();
      ASSERT_TRUE(vecTrait::isVector == 1);
      ASSERT_TRUE(vecTrait::isEigen == 1);
      ASSERT_TRUE(vecTrait::isSerial == 1);
      ASSERT_TRUE(vecTrait::isSTDVector == 0);
      ASSERT_TRUE(vecTrait::isDistributed == 0);
      if (row == 1 && col != Eigen::Dynamic) // row vector
	ASSERT_TRUE(vecTrait::isStatic == 1);
      else if (col == 1 && row != Eigen::Dynamic) // col vector
	ASSERT_TRUE(vecTrait::isStatic == 1);
      else
	ASSERT_TRUE(vecTrait::isStatic == 0);
    }
  };
  
  static constexpr int dyn = Eigen::Dynamic;

  // row vectors
  EigenVecChecker<int,1,dyn> a;
  EigenVecChecker<double,1,dyn> b;
  EigenVecChecker<float,1,dyn> c;
  EigenVecChecker<std::complex<double>,1,dyn> d;
  EigenVecChecker<std::complex<int>,1,dyn> e;

  EigenVecChecker<int,1,5> a1;
  EigenVecChecker<double,1,5> b1;
  EigenVecChecker<float,1,5> c1;
  EigenVecChecker<std::complex<double>,1,5> d1;
  EigenVecChecker<std::complex<int>,1,5> e1;
  
  // column vectors
  EigenVecChecker<int,dyn,1> f;
  EigenVecChecker<double,dyn,1> g;
  EigenVecChecker<float,dyn,1> h;
  EigenVecChecker<std::complex<double>,dyn,1> k;
  EigenVecChecker<std::complex<int>,dyn,1> p;

  EigenVecChecker<int,4,1> f1;
  EigenVecChecker<double,4,1> g1;
  EigenVecChecker<float,4,1> h1;
  EigenVecChecker<std::complex<double>,4,1> k1;
  EigenVecChecker<std::complex<int>,4,1> p1;
};

TEST_F(core_vector_serial_eigen_traits_Fixture, traits)
{
  a.check(); b.check(); c.check(); d.check(); e.check();
  a1.check(); b1.check(); c1.check(); d1.check(); e1.check();
  f.check(); g.check(); h.check(); k.check(); p.check();
  f1.check(); g1.check(); h1.check(); k1.check(); p1.check();
  
  // check that a matrix from eigen is not a vector
  using eigmat_t = Eigen::Matrix<double, dyn, dyn>;
  STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(eigmat_t);
  using eigmat_t2 = Eigen::Matrix<double, 4, dyn>;
  STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(eigmat_t2);
  using eigmat_t3 = Eigen::Matrix<double, dyn, 5>;
  STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(eigmat_t3);
}
