
#include <gtest/gtest.h>
#include "CONTAINERS_VECTOR"

TEST(containers_vector_dot_product, eigenVector)
{
  using natV_t = Eigen::Matrix<double,4,1>;
  STATIC_ASSERT_IS_VECTOR_EIGEN(natV_t);  
  natV_t a; a << 3.,4.,1.,2.;
  natV_t b; b << 1.,1.,1.,1.;
  using my_t = containers::Vector<natV_t>;
  my_t va(a);
  my_t vb(b);

  // double res; containers::vec_ops::dot(va,vb,res);
  // EXPECT_DOUBLE_EQ( res, 10.0);
}//end TEST


// TEST(containers_vector_dot_product_DeathTest, eigenVector) {
//   using natV1_t = Eigen::Matrix<double,Eigen::Dynamic,1>;
//   natV1_t a; a << 1.,3.,4.,7;
//   containers::Vector<natV1_t> va(a);

//   using natV1_t2 = Eigen::Matrix<double,Eigen::Dynamic,1>;
//   natV1_t2 b; a << 1.,3.,4.,7., 3., 10.;
//   containers::Vector<natV1_t2> vb(b);

//   double res;
//   ASSERT_DEATH( containers::dotProduct(va,vb,res), "Assertion failed:*");// (vecA.size() == vecB.size())*");
// }
