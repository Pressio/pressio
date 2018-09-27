
#include <gtest/gtest.h>
#include "CORE_VECTOR"

TEST(core_vector_blaze_class, constructor)
{
  using namespace rompp;

  using nvec_t = blaze::DynamicVector<double>;
  using myvec_t = core::Vector<nvec_t>;
  using vecTrait = core::details::traits<myvec_t>;
  ASSERT_TRUE(vecTrait::wrapped_vector_identifier 
    == core::details::WrappedVectorIdentifier::BlazeDynamic);

  STATIC_ASSERT_IS_VECTOR_BLAZE(nvec_t);
  STATIC_ASSERT_IS_NOT_VECTOR_BLAZE(myvec_t);
    
  myvec_t m_v2(5);
  ASSERT_TRUE( m_v2.size() == 5 );

  // create an native vector
  nvec_t e_v1(4);
  myvec_t m_v3(e_v1);
  ASSERT_TRUE( m_v3.size() == 4 );
}


// TEST(core_vector_serial_eigen_class, queryWrappedData)
// {
// 	using namespace rompp;

//   using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
//   using myvec_t = core::Vector<eigvec_t>;

//   myvec_t m_v1(4);
//   ::testing::StaticAssertTypeEq<decltype(m_v1.data()),
// 				eigvec_t * >(); 
//   const myvec_t m_v2(4);
//   ::testing::StaticAssertTypeEq< decltype(m_v2.data()),
// 				 const eigvec_t * >();
// }


TEST(core_vector_serial_eigen_class, sizeResize)
{
  using namespace rompp;
  using nvec_t = blaze::DynamicVector<double>;
  using myvec_t = core::Vector<nvec_t>; 

  myvec_t m_v1(11);
  ASSERT_TRUE( m_v1.size() == 11 );
  myvec_t m_v2(22);
  ASSERT_TRUE( m_v2.size() == 22 );
  m_v2.resize(33);
  ASSERT_TRUE( m_v2.size() == 33 );  
}


TEST(core_vector_serial_eigen_class, subscriptOperator)
{
  using namespace rompp;
  using nvec_t = blaze::DynamicVector<double>;
  using myvec_t = core::Vector<nvec_t>; 

  // create an native vector
  nvec_t e_v1{34.0, 22.5, 11.5, 75.};

  myvec_t m_v3(e_v1);
  ASSERT_TRUE( m_v3.size() == 4 );
  EXPECT_DOUBLE_EQ( m_v3[0], 34.0);
  EXPECT_DOUBLE_EQ( m_v3[1], 22.5);
  EXPECT_DOUBLE_EQ( m_v3[2], 11.5);
  EXPECT_DOUBLE_EQ( m_v3[3], 75.0);

  m_v3[0] = 56.;
  m_v3[3] = 44.;
  EXPECT_DOUBLE_EQ( m_v3[0], 56.0);
  EXPECT_DOUBLE_EQ( m_v3[3], 44.0);  
}


TEST(core_vector_serial_eigen_class, additionOperator)
{
  using namespace rompp;
  using nvec_t = blaze::DynamicVector<double>;
  using myvec_t = core::Vector<nvec_t>; 

  
  myvec_t v1(4);
  v1[0] = 3.; v1[1] = 2.;
  v1[2] = 4.; v1[3] = 5.;
  myvec_t v2(4);
  v2[0] = 1.; v2[1] = 1.;
  v2[2] = 1.; v2[3] = 1.;

  myvec_t c(v1 + v2);
  EXPECT_DOUBLE_EQ( c[0], 4.0);
  EXPECT_DOUBLE_EQ( c[1], 3.);
  EXPECT_DOUBLE_EQ( c[2], 5.);
  EXPECT_DOUBLE_EQ( c[3], 6.0);

  myvec_t c2(v1 + v2 + v2);
  EXPECT_DOUBLE_EQ( c2[0], 5.0);
  EXPECT_DOUBLE_EQ( c2[1], 4.);
  EXPECT_DOUBLE_EQ( c2[2], 6.);
  EXPECT_DOUBLE_EQ( c2[3], 7.0);
  
}
//   // myvec_t v4 = (v1+v2)*2.;
//   // v4 = 3.*(v1+v2);
//   // v4 = 3.*v2;
//   // v4 = v2*2.;
  
//   // myvec_t v4(4);
//   // v4 = v1 + v2 + v1;
//   // std::cout << *v4.data() << "\n";

//   // myvec_t v5(v1 + v2);
//   // std::cout << *v5.data() << "\n";

//   myvec_t v6 = v1;
//   std::cout << *v6.data() << "\n";
  
//   // v4 = v2 - v1;
//   // v4.data()->Print(std::cout);
  
//   // v4 = v2 + v1 + v3;
//   // v4.data()->Print(std::cout);
  
//   // v4 = v2 + v1 + v3;
//   // v4.data()->Print(std::cout);

//   // v4 = v2*2. + 1.*v1;
//   // v4.data()->Print(std::cout);

//   // v4 = v2*2. + 1.*v1 + v3;
  
//   // myvec_t res = m_v1 + m_v1 - m_v2;  
//   // std::cout << *res.data() << "\n";
  
//   // EXPECT_DOUBLE_EQ(res[0], 4.);
//   // EXPECT_DOUBLE_EQ(res[1], 3.);
//   // EXPECT_DOUBLE_EQ(res[2], 5.);
//   // EXPECT_DOUBLE_EQ(res[3], 6.);
// }

// // TEST(core_vector_serial_eigen_class, substractOperator)
// // {
// // 	using namespace rompp;
// //   using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
// //   using myvec_t = core::Vector<eigvec_t>;

// //   myvec_t m_v1(4);
// //   m_v1[0] = 3.; m_v1[1] = 2.;
// //   m_v1[2] = 4.; m_v1[3] = 5.;
// //   myvec_t m_v2(4);
// //   m_v2[0] = 1.; m_v2[1] = 1.;
// //   m_v2[2] = 1.; m_v2[3] = 1.;
 
// //   myvec_t res = m_v1 - m_v2;
// //   EXPECT_DOUBLE_EQ(res[0], 2.);
// //   EXPECT_DOUBLE_EQ(res[1], 1.);
// //   EXPECT_DOUBLE_EQ(res[2], 3.);
// //   EXPECT_DOUBLE_EQ(res[3], 4.);
// // }

// // TEST(core_vector_serial_eigen_class, starOperator)
// // {
// // 	using namespace rompp;
// //   using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
// //   using myvec_t = core::Vector<eigvec_t>;

// //   myvec_t m_v1(4);
// //   m_v1[0] = 3.; m_v1[1] = 2.;
// //   m_v1[2] = 4.; m_v1[3] = 5.;
// //   myvec_t m_v2(4);
// //   m_v2[0] = 1.; m_v2[1] = 1.;
// //   m_v2[2] = 1.; m_v2[3] = 1.;

// //   myvec_t res = m_v1 * m_v2;
// //   EXPECT_DOUBLE_EQ(res[0], 3.);
// //   EXPECT_DOUBLE_EQ(res[1], 2.);
// //   EXPECT_DOUBLE_EQ(res[2], 4.);
// //   EXPECT_DOUBLE_EQ(res[3], 5.);
// // }

// TEST(core_vector_serial_eigen_class, CompoundAssignAddOperator)
// {
// 	using namespace rompp;
//   using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
//   using myvec_t = core::Vector<eigvec_t>;

//   myvec_t m_v1(4);
//   m_v1[0] = 3.; m_v1[1] = 2.;
//   m_v1[2] = 4.; m_v1[3] = 5.;
//   myvec_t m_v2(4);
//   m_v2[0] = 1.; m_v2[1] = 1.;
//   m_v2[2] = 1.; m_v2[3] = 1.;

//   m_v1 += m_v2;
//   EXPECT_DOUBLE_EQ(m_v1[0], 4.);
//   EXPECT_DOUBLE_EQ(m_v1[1], 3.);
//   EXPECT_DOUBLE_EQ(m_v1[2], 5.);
//   EXPECT_DOUBLE_EQ(m_v1[3], 6.);
// }


// TEST(core_vector_serial_eigen_class, CompoundAssignSubtractOperator)
// {
// 	using namespace rompp;
//   using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
//   using myvec_t = core::Vector<eigvec_t>;

//   myvec_t m_v1(4);
//   m_v1[0] = 3.; m_v1[1] = 2.;
//   m_v1[2] = 4.; m_v1[3] = 5.;
//   myvec_t m_v2(4);
//   m_v2[0] = 1.; m_v2[1] = 1.;
//   m_v2[2] = 1.; m_v2[3] = 1.;

//   m_v1 -= m_v2;
//   EXPECT_DOUBLE_EQ(m_v1[0], 2.);
//   EXPECT_DOUBLE_EQ(m_v1[1], 1.);
//   EXPECT_DOUBLE_EQ(m_v1[2], 3.);
//   EXPECT_DOUBLE_EQ(m_v1[3], 4.);
// }
