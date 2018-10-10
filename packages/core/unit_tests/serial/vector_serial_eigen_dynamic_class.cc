
#include <gtest/gtest.h>
#include "CORE_VECTOR"

using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using myvec_t = rompp::core::Vector<eigvec_t>;


TEST(core_vector_serial_eigen_dynamic_class,
     constructor){

  using vecTrait = rompp::core::details::traits<myvec_t>;
  ASSERT_TRUE(vecTrait::wrapped_vector_identifier 
	      == rompp::core::details::WrappedVectorIdentifier::Eigen);

  //construct by passing the size 
  myvec_t m_v2(5);
  ASSERT_FALSE( m_v2.empty() );
  ASSERT_TRUE( m_v2.size() == 5 );
  for (size_t i=0; i<5; i++)
    EXPECT_DOUBLE_EQ( m_v2[i], 0.);
  
  // pass native eigen vector
  eigvec_t e_v1(45);
  e_v1(2) = 2.2; e_v1(4) = 4.4;  
  myvec_t m_v1(e_v1);
  ASSERT_FALSE( m_v1.empty() );
  ASSERT_TRUE( m_v1.size() == 45 );
  EXPECT_DOUBLE_EQ( m_v1[2], 2.2);
  EXPECT_DOUBLE_EQ( m_v1[4], 4.4);
}


TEST(core_vector_serial_eigen_dynamic_class,
     queryWrappedData){

  myvec_t m_v1(4);
  ::testing::StaticAssertTypeEq<decltype(m_v1.data()),
				eigvec_t * >(); 
  const myvec_t m_v2(4);
  ::testing::StaticAssertTypeEq< decltype(m_v2.data()),
				 const eigvec_t * >();
}


TEST(core_vector_serial_eigen_dynamic_class,
     size){

  myvec_t m_v1(11);
  ASSERT_FALSE( m_v1.empty() );
  ASSERT_TRUE( m_v1.size() == 11 );
}


TEST(core_vector_serial_eigen_dynamic_class,
     resize){

  myvec_t m_v1(11);
  ASSERT_FALSE( m_v1.empty() );
  ASSERT_TRUE( m_v1.size() == 11 );
  m_v1.resize(22);
  ASSERT_FALSE( m_v1.empty() );
  ASSERT_TRUE( m_v1.size() == 22 );
}


TEST(core_vector_serial_eigen_dynamic_class,
     subscriptOperatorSquareBrack){

  myvec_t m_v3(4);
  ::testing::StaticAssertTypeEq<
    decltype(m_v3[1]), double & >();
  
  ASSERT_TRUE( m_v3.size() == 4 );
  m_v3[0] = 34.0;
  m_v3[1] = 22.5;
  m_v3[2] = 11.5;
  m_v3[3] = 75.0;  
  EXPECT_DOUBLE_EQ( m_v3[0], 34.0);
  EXPECT_DOUBLE_EQ( m_v3[1], 22.5);
  EXPECT_DOUBLE_EQ( m_v3[2], 11.5);
  EXPECT_DOUBLE_EQ( m_v3[3], 75.0);
  m_v3[0] = 56.;
  m_v3[3] = 44.;
  EXPECT_DOUBLE_EQ( m_v3[0], 56.0);
  EXPECT_DOUBLE_EQ( m_v3[3], 44.0);  

  const myvec_t m_v4(4);
  ::testing::StaticAssertTypeEq<
    decltype(m_v4[1]), const double & >();
}

TEST(core_vector_serial_eigen_dynamic_class,
     subscriptOperatorParenth){

  myvec_t m_v3(4);
  ASSERT_TRUE( m_v3.size() == 4 );
  ::testing::StaticAssertTypeEq<
    decltype(m_v3(1)), double & >();
  m_v3(0) = 34.0;
  m_v3(1) = 22.5;
  m_v3(2) = 11.5;
  m_v3(3) = 75.0;  
  EXPECT_DOUBLE_EQ( m_v3(0), 34.0);
  EXPECT_DOUBLE_EQ( m_v3(1), 22.5);
  EXPECT_DOUBLE_EQ( m_v3(2), 11.5);
  EXPECT_DOUBLE_EQ( m_v3(3), 75.0);

  m_v3(0) = 56.;
  m_v3(3) = 44.;
  EXPECT_DOUBLE_EQ( m_v3(0), 56.0);
  EXPECT_DOUBLE_EQ( m_v3(3), 44.0);

  const myvec_t m_v4(4);
  ::testing::StaticAssertTypeEq<
    decltype(m_v4(1)), const double & >();  
}

TEST(core_vector_serial_eigen_dynamic_class,
     matchingLayout){

  myvec_t a(4);
  ASSERT_TRUE( a.size() == 4 );
  myvec_t b(6);
  a.matchLayoutWith(b);
  ASSERT_TRUE( a.size() == 6 );
  ASSERT_FALSE( a.size() == 4 );
}

TEST(core_vector_serial_eigen_dynamic_class,
     setToScalar){

  myvec_t a(4);
  ASSERT_TRUE( a.size() == 4 );
  a.putScalar(1.12);
  EXPECT_DOUBLE_EQ( a(0), 1.12);
  EXPECT_DOUBLE_EQ( a(1), 1.12);
  EXPECT_DOUBLE_EQ( a(2), 1.12);
  EXPECT_DOUBLE_EQ( a(3), 1.12);
}

TEST(core_vector_serial_eigen_dynamic_class,
     assignScalar){

  myvec_t a(4);
  a = 1.12;
  EXPECT_DOUBLE_EQ( a(0), 1.12);
  EXPECT_DOUBLE_EQ( a(1), 1.12);
  EXPECT_DOUBLE_EQ( a(2), 1.12);
  EXPECT_DOUBLE_EQ( a(3), 1.12);
}



TEST(core_vector_serial_eigen_dynamic_class,
     settingZero){

  myvec_t a(4);
  a.setZero();
  EXPECT_DOUBLE_EQ( a(0), 0.0);
  EXPECT_DOUBLE_EQ( a(1), 0.0);
  EXPECT_DOUBLE_EQ( a(2), 0.0);
  EXPECT_DOUBLE_EQ( a(3), 0.0);
}


TEST(core_vector_serial_eigen_dynamic_class,
     scaleByFactor){

  myvec_t a(4);
  a = 4.;
  a.scale(2.);
  EXPECT_DOUBLE_EQ( a(0), 8.0);
  EXPECT_DOUBLE_EQ( a(1), 8.0);
  EXPECT_DOUBLE_EQ( a(2), 8.0);
  EXPECT_DOUBLE_EQ( a(3), 8.0);
}





// TEST(core_vector_serial_eigen_class, additionOperator)
// {
// 	using namespace rompp;
//   using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
//   using myvec_t = core::Vector<eigvec_t>;
  
//   myvec_t v1(4);
//   v1[0] = 3.; v1[1] = 2.;
//   v1[2] = 4.; v1[3] = 5.;
//   myvec_t v2(4);
//   v2[0] = 1.; v2[1] = 1.;
//   v2[2] = 1.; v2[3] = 1.;

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
