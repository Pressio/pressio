
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"


TEST(containers_ops_dot,
     eigen_dynamic_vector_dot_same){

  using nat_t = Eigen::VectorXd;

  pressio::containers::Vector<nat_t> a(6);
  a(0) = 1.; a(1) = 1.; a(2) = 1.;
  a(3) = 1.; a(4) = 2.; a(5) = 1.;

  pressio::containers::Vector<nat_t> b(6);
  b(0) = 1.; b(1) = 1.; b(2) = 1.;
  b(3) = 1.; b(4) = 2.; b(5) = 1.;

  auto res = pressio::containers::ops::dot(a,b);
  EXPECT_DOUBLE_EQ( res, 9.);

  double res2 = {};
  pressio::containers::ops::dot(a,b,res2);
  EXPECT_DOUBLE_EQ( res2, 9.);
}



TEST(containers_ops_dot,
     eigen_static_vector_dot_eigen_static_vector){

  using nat_t = Eigen::Matrix<double, 6, 1>;

  pressio::containers::Vector<nat_t> a;
  a(0) = 1.; a(1) = 1.; a(2) = 1.;
  a(3) = 1.; a(4) = 2.; a(5) = 1.;

  pressio::containers::Vector<nat_t> b;
  b(0) = 1.; b(1) = 1.; b(2) = 1.;
  b(3) = 1.; b(4) = 2.; b(5) = 1.;

  auto res = pressio::containers::ops::dot(a,b);
  EXPECT_DOUBLE_EQ( res, 9.);

  double res2 = {};
  pressio::containers::ops::dot(a,b,res2);
  EXPECT_DOUBLE_EQ( res2, 9.);
}
