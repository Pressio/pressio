
#include <gtest/gtest.h>
#include "CORE_ALL"


TEST(core_ops_dot,
     eigen_dynamic_vector_dot_same){

  using nat_t = Eigen::VectorXd;

  rompp::core::Vector<nat_t> a(6);
  a(0) = 1.; a(1) = 1.; a(2) = 1.;
  a(3) = 1.; a(4) = 2.; a(5) = 1.;

  rompp::core::Vector<nat_t> b(6);
  b(0) = 1.; b(1) = 1.; b(2) = 1.;
  b(3) = 1.; b(4) = 2.; b(5) = 1.;

  auto res = rompp::core::ops::dot(a,b);
  EXPECT_DOUBLE_EQ( res, 9.);

  double res2 = {};
  rompp::core::ops::dot(a,b,res2);
  EXPECT_DOUBLE_EQ( res2, 9.);
}



TEST(core_ops_dot,
     eigen_static_vector_dot_eigen_static_vector){

  using nat_t = Eigen::Matrix<double, 6, 1>;

  rompp::core::Vector<nat_t> a;
  a(0) = 1.; a(1) = 1.; a(2) = 1.;
  a(3) = 1.; a(4) = 2.; a(5) = 1.;

  rompp::core::Vector<nat_t> b;
  b(0) = 1.; b(1) = 1.; b(2) = 1.;
  b(3) = 1.; b(4) = 2.; b(5) = 1.;

  auto res = rompp::core::ops::dot(a,b);
  EXPECT_DOUBLE_EQ( res, 9.);

  double res2 = {};
  rompp::core::ops::dot(a,b,res2);
  EXPECT_DOUBLE_EQ( res2, 9.);
}
