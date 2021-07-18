
#include <gtest/gtest.h>
#include "pressio_ops.hpp"

TEST(ops_eigen, clone_vector)
{
  using T = Eigen::VectorXd;
  T a(6);
  for (int i=0; i<6; ++i){
   a(i)= (double) i;
  }

  auto b = pressio::ops::clone(a);
  ASSERT_EQ(b.size(), 6);
  ASSERT_FALSE( b.data()==a.data());
  for (int i=0; i<6; ++i){
    ASSERT_DOUBLE_EQ(b(i),a(i));
  }
}

TEST(ops_eigen, clone_matrix)
{
  using T = Eigen::MatrixXd;
  T a(6,8);

  int c=0;  
  for (int i=0; i<6; ++i){
    for (int j=0; j<8; ++j){
     a(i)= (double) ++c;
   }
  }

  auto b = pressio::ops::clone(a);
  ASSERT_EQ(b.rows(), 6);
  ASSERT_EQ(b.cols(), 8);
  ASSERT_FALSE( b.data()==a.data());

  for (int i=0; i<6; ++i){
    for (int j=0; j<8; ++j){
      ASSERT_DOUBLE_EQ(b(i,j),a(i,j));
     a(i)= (double) ++c;
   }
  }
}

TEST(ops_eigen, extent_vector)
{
  using vec_t = Eigen::VectorXd;
  vec_t x(6);
  ASSERT_TRUE(pressio::ops::extent(x,0)== 6);
}

TEST(ops_eigen, extent_matrix)
{
  using T = Eigen::MatrixXd;
  T x(6,8);
  ASSERT_TRUE(pressio::ops::extent(x,0) == 6);
  ASSERT_TRUE(pressio::ops::extent(x,1) == 8);
}

TEST(ops_eigen, abs)
{
  using vec_t = Eigen::VectorXd;
  vec_t x(6);
  for (int i=0; i<6; ++i){
    x(i) = -(double) i;
  }

  vec_t y(6);
  pressio::ops::abs(y,x);

  Eigen::VectorXd g(6);
  g(0) = 0.;
  g(1) = 1.;
  g(2) = 2.;
  g(3) = 3.;
  g(4) = 4.;
  g(5) = 5.;
  ASSERT_EQ( y.size(), 6 );
  for (int i=0; i<6; ++i){
    EXPECT_DOUBLE_EQ(y(i), g(i));
  }
}

TEST(ops_eigen, scale_vector)
{
  using T = Eigen::VectorXd;
  T a(6);
  a.setConstant(1.);

  pressio::ops::scale(a, 3.);
  for (int i=0; i<6; ++i){
    ASSERT_DOUBLE_EQ(a(i),3.);
  }
}

TEST(ops_eigen, scale_matrix)
{
  using T = Eigen::MatrixXd;
  T a(6,8);
  a.setConstant(1.);

  pressio::ops::scale(a, 3.);
  for (int i=0; i<6; ++i){
   for (int j=0; j<8; ++j){
    ASSERT_DOUBLE_EQ(a(i,j),3.);
   }
  }
}

TEST(ops_eigen, setzero_vector)
{
  using T = Eigen::VectorXd;
  T a(6);
  a.setConstant(11.);

  pressio::ops::set_zero(a);
  for (int i=0; i<6; ++i){
    ASSERT_DOUBLE_EQ(a(i),0.);
  }
}

TEST(ops_eigen, setzero_matrix)
{
  using T = Eigen::MatrixXd;
  T a(6,6);
  a.setConstant(11.);

  pressio::ops::set_zero(a);
  for (int i=0; i<6; ++i){
   for (int j=0; j<6; ++j){
    ASSERT_DOUBLE_EQ(a(i,j),0.);
   }
  }
}

TEST(ops_eigen, fill_vector)
{
  using T = Eigen::VectorXd;
  T a(6);
  a.setConstant(11.);

  pressio::ops::fill(a, 44.88);
  for (int i=0; i<6; ++i){
    ASSERT_DOUBLE_EQ(a(i),44.88);
  }
}

TEST(ops_eigen, fill_matrix)
{
  using T = Eigen::MatrixXd;
  T a(6,6);
  a.setConstant(11.);

  pressio::ops::fill(a, 55.);
  for (int i=0; i<6; ++i){
   for (int j=0; j<6; ++j){
    ASSERT_DOUBLE_EQ(a(i,j), 55.);
   }
  }
}

TEST(ops_eigen, resize_vector)
{
  using T = Eigen::VectorXd;
  T a(6);

  pressio::ops::resize(a,3);
  ASSERT_EQ(a.size(), 3);
}

TEST(ops_eigen, resize_matrix)
{
  using T = Eigen::MatrixXd;
  T a(6,6);
  pressio::ops::resize(a,3,4);
  ASSERT_EQ(a.rows(), 3);
  ASSERT_EQ(a.cols(), 4);
}

TEST(ops_eigen, deep_copy_vector)
{
  using T = Eigen::VectorXd;
  T a(6);
  a.setConstant(44.);

  T b(6);
  pressio::ops::deep_copy(b,a);
  for (int i=0; i<6; ++i){
    ASSERT_DOUBLE_EQ(a(i),44.);
  }
}

TEST(ops_eigen, deep_copy_matrix)
{
  using T = Eigen::MatrixXd;
  T a(6,5);
  a.setConstant(44.);

  T b(6,5);
  pressio::ops::deep_copy(b,a);
  for (int i=0; i<6; ++i){
    ASSERT_DOUBLE_EQ(a(i),44.);
  }
}

TEST(ops_eigen, norm1_vector)
{
  using T = Eigen::VectorXd;
  T a(5);
  for (int i=0; i<5; ++i){
   a(i)= (double) i;
  }
  ASSERT_DOUBLE_EQ(pressio::ops::norm1(a), a.lpNorm<1>());
}

TEST(ops_eigen, norm2_vector)
{
  using T = Eigen::VectorXd;
  T a(5);
  for (int i=0; i<5; ++i){
   a(i)= (double) i;
  }
 ASSERT_DOUBLE_EQ(pressio::ops::norm2(a), a.norm());
 ASSERT_DOUBLE_EQ(pressio::ops::norm2(a), a.lpNorm<2>());
}

TEST(ops_eigen, dot_vector)
{
  using T = Eigen::VectorXd;
  T a(5);
  a.setConstant(1.);
  T b(5);
  b.setConstant(2.);

  ASSERT_DOUBLE_EQ(pressio::ops::dot(a,b), 10.);
}


TEST(ops_eigen, vectorPow)
{
  using vec_t = Eigen::VectorXd;
  vec_t x(6);
  for (int i=0; i<6; ++i) x(i) = (double) i;

  ::pressio::ops::pow(x, 2.);
  Eigen::VectorXd g(6);
  g(0) = 0.;
  g(1) = 1.;
  g(2) = 4.;
  g(3) = 9.;
  g(4) = 16.;
  g(5) = 25.;
  ASSERT_EQ( x.size(), 6 );
  for (int i=0; i<6; ++i)
    EXPECT_DOUBLE_EQ(x(i), g(i));
}

TEST(ops_eigen, vectorAbsPowPos)
{
  using vec_t = Eigen::VectorXd;
  vec_t y(6);
  vec_t x(6);
  for (int i=0; i<6; ++i){
    x(i) = (double) i; x(i)*=-1.;
  }

  ::pressio::ops::abs_pow(y, x, 3.);
  Eigen::VectorXd g(6);
  g(0) = 0.;
  g(1) = 1.;
  g(2) = 8.;
  g(3) = 27.;
  g(4) = 64.;
  g(5) = 125.;
  ASSERT_EQ( y.size(), 6 );
  for (int i=0; i<6; ++i)
    EXPECT_DOUBLE_EQ(y(i), g(i));
}

TEST(ops_eigen, vectorAbsPowNeg)
{
  using vec_t = Eigen::VectorXd;
  vec_t y(6);
  vec_t x(6);
  for (int i=0; i<6; ++i){
    x(i) = (double) i; 
    x(i)*=-1.;
  }

  ::pressio::ops::abs_pow(y, x, -3., 0.00001);
  // std::cout << y << std::endl;

  Eigen::VectorXd g(6);
  g(0) = 1./0.00001; // because we guard against diving by zero
  g(1) = 1.;
  g(2) = 1./8.;
  g(3) = 1./27.;
  g(4) = 1./64.;
  g(5) = 1./125.;
  ASSERT_EQ( y.size(), 6 );
  for (int i=0; i<6; ++i)
    EXPECT_DOUBLE_EQ(y(i), g(i));
}


TEST(ops_eigen, VectorUpdate1)
{
  using V_t = Eigen::Matrix<double,3,1>;
  V_t v; v << 1.,2.,3.;
  V_t a; a << 1.,2.,3.;

  pressio::ops::update(v, 1., a, 1.);
  EXPECT_DOUBLE_EQ( v(0), 2.0);
  EXPECT_DOUBLE_EQ( v(1), 4.0);
  EXPECT_DOUBLE_EQ( v(2), 6.0);

  pressio::ops::update(v, a, 1.);
  EXPECT_DOUBLE_EQ( v(0), 1.0);
  EXPECT_DOUBLE_EQ( v(1), 2.0);
  EXPECT_DOUBLE_EQ( v(2), 3.0);
}

TEST(ops_eigen, VectorUpdate2)
{
  using V_t = Eigen::Matrix<double,3,1>;
  V_t v; v << 1.,2.,3.;
  V_t a; a << 1.,2.,3.;
  V_t b; b << 1.,2.,3.;

  pressio::ops::update(v, 1., a, 1., b, 1.);
  EXPECT_DOUBLE_EQ( v(0), 3.0);
  EXPECT_DOUBLE_EQ( v(1), 6.0);
  EXPECT_DOUBLE_EQ( v(2), 9.0);

  pressio::ops::update(v, a, 1., b, 1.);
  EXPECT_DOUBLE_EQ( v(0), 2.0);
  EXPECT_DOUBLE_EQ( v(1), 4.0);
  EXPECT_DOUBLE_EQ( v(2), 6.0);
}

TEST(ops_eigen, VectorUpdate3)
{
  using V_t = Eigen::Matrix<double,3,1>;
  V_t v; v << 1.,2.,3.;
  V_t a; a << 1.,2.,3.;
  V_t b; b << 1.,2.,3.;
  V_t c; c << 1.,2.,3.;

  pressio::ops::update(v, 1., a, 1., b, 1., c, 1.);
  EXPECT_DOUBLE_EQ( v(0), 4.0);
  EXPECT_DOUBLE_EQ( v(1), 8.0);
  EXPECT_DOUBLE_EQ( v(2), 12.0);

  pressio::ops::update(v, a, 1., b, 1., c, 1.);
  EXPECT_DOUBLE_EQ( v(0), 3.0);
  EXPECT_DOUBLE_EQ( v(1), 6.0);
  EXPECT_DOUBLE_EQ( v(2), 9.0);
}

TEST(ops_eigen, VectorUpdate4)
{
  using V_t = Eigen::Matrix<double,3,1>;
  V_t v; v << 1.,2.,3.;
  V_t a; a << 1.,2.,3.;
  V_t b; b << 1.,2.,3.;
  V_t c; c << 1.,2.,3.;
  V_t d; d << 1.,2.,3.;

  pressio::ops::update(v, 1., a, 1., b, 1., c, 1., d, 1.);
  EXPECT_DOUBLE_EQ( v(0), 5.0);
  EXPECT_DOUBLE_EQ( v(1), 10.0);
  EXPECT_DOUBLE_EQ( v(2), 15.0);

  pressio::ops::update(v, a, 1., b, 1., c, 1., d, 1.);
  EXPECT_DOUBLE_EQ( v(0), 4.0);
  EXPECT_DOUBLE_EQ( v(1), 8.0);
  EXPECT_DOUBLE_EQ( v(2), 12.0);
}

TEST(ops_eigen, elementwiseMultiply)
{
  using V_t = Eigen::Matrix<double,3,1>;
  V_t y; y << 1.,2.,3.;
  V_t x; x << 2.,3.,4.;
  V_t z; z << 3.,4.,5.;

  pressio::ops::elementwise_multiply(1., x, z, 1., y);
  EXPECT_DOUBLE_EQ( y(0), 7.0);
  EXPECT_DOUBLE_EQ( y(1), 14.0);
  EXPECT_DOUBLE_EQ( y(2), 23.0);
}

TEST(ops_eigen, VectorDenseMatrixProd)
{
  using V_t = Eigen::VectorXd;
  V_t a(3); a << 4.,2.,6;

  using M_t = Eigen::MatrixXd;
  M_t M(3,3);
  M << 1,0,2,2,1,3,0,0,1;

  V_t myR(3);
  constexpr auto beta  = ::pressio::utils::constants<double>::zero();
  constexpr auto alpha = ::pressio::utils::constants<double>::one();
  pressio::ops::product(pressio::nontranspose(), alpha, M, a, beta, myR);
  EXPECT_DOUBLE_EQ( myR(0), 16.0);
  EXPECT_DOUBLE_EQ( myR(1), 28.0);
  EXPECT_DOUBLE_EQ( myR(2), 6.0);
}
