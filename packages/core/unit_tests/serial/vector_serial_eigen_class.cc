
#include <gtest/gtest.h>
#include "CORE_VECTOR"

TEST(core_vector_serial_eigen_class, EigenVectorConstructor)
{
	using namespace rompp;

  // using eigvec_t1 = Eigen::Matrix<double, 4, 1>;
  // using myvec_t1 = core::Vector<eigvec_t1>;
  // myvec_t1 m333(5);
  using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  using myvec_t = core::Vector<eigvec_t>;
  using vecTrait = core::details::traits<myvec_t>;
  ASSERT_TRUE(vecTrait::wrapped_vector_identifier 
    == core::details::WrappedVectorIdentifier::Eigen);
 
  myvec_t m_v2(5);
  ASSERT_FALSE( m_v2.empty() );
  ASSERT_TRUE( m_v2.size() == 5 );

  // create an eigen-type vector
  Eigen::Vector4d e_v1;
  myvec_t m_v3(e_v1);
  ASSERT_FALSE( m_v3.empty() );
  ASSERT_TRUE( m_v3.size() == 4 );
}


TEST(core_vector_serial_eigen_class, queryWrappedData)
{
	using namespace rompp;

  using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  using myvec_t = core::Vector<eigvec_t>;

  myvec_t m_v1(4);
  ::testing::StaticAssertTypeEq<decltype(m_v1.data()),
				eigvec_t * >(); 
  const myvec_t m_v2(4);
  ::testing::StaticAssertTypeEq< decltype(m_v2.data()),
				 const eigvec_t * >();
}


TEST(core_vector_serial_eigen_class, sizeResize)
{
	using namespace rompp;
  using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  using myvec_t = core::Vector<eigvec_t>;
 
  // myvec_t m_v1;
  // ASSERT_TRUE( m_v1.empty() );
  // ASSERT_TRUE( m_v1.size() == 0 );
  myvec_t m_v1(11);
  ASSERT_FALSE( m_v1.empty() );
  ASSERT_TRUE( m_v1.size() == 11 );

  myvec_t m_v2(22);
  ASSERT_FALSE( m_v2.empty() );
  ASSERT_TRUE( m_v2.size() == 22 );
}


TEST(core_vector_serial_eigen_class, subscriptOperator)
{
  using namespace rompp;
  using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  using myvec_t = core::Vector<eigvec_t>;

  // create an eigen-type vector
  Eigen::Vector4d e_v1;
  e_v1 << 34.0, 22.5, 11.5, 75.;

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
  using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  using myvec_t = core::Vector<eigvec_t>;
  
  myvec_t m_v1(4);
  m_v1[0] = 3.; m_v1[1] = 2.;
  m_v1[2] = 4.; m_v1[3] = 5.;
  myvec_t m_v2(4);
  m_v2[0] = 1.; m_v2[1] = 1.;
  m_v2[2] = 1.; m_v2[3] = 1.;

  // v4 = (v1+v2)*2.;
  // v4.data()->Print(std::cout);

  // v4 = 3.*(v1+v2);
  // v4.data()->Print(std::cout);
  
  // v4 = 3.*v2;
  // v4.data()->Print(std::cout);

  // v4 = v2*2.;
  // v4.data()->Print(std::cout);
  
  myvec_t v4(4);
  v4 = m_v1 + m_v2 + m_v1;
  std::cout << *v4.data() << "\n";

  myvec_t v5(m_v1 + m_v2);
  std::cout << *v5.data() << "\n";

  v5 = m_v1 - m_v2;
  std::cout << *v5.data() << "\n";
  
  // v4 = v2 - v1;
  // v4.data()->Print(std::cout);
  
  // v4 = v2 + v1 + v3;
  // v4.data()->Print(std::cout);
  
  // v4 = v2 + v1 + v3;
  // v4.data()->Print(std::cout);

  // v4 = v2*2. + 1.*v1;
  // v4.data()->Print(std::cout);

  // v4 = v2*2. + 1.*v1 + v3;
  
  // myvec_t res = m_v1 + m_v1 - m_v2;  
  // std::cout << *res.data() << "\n";
  
  // EXPECT_DOUBLE_EQ(res[0], 4.);
  // EXPECT_DOUBLE_EQ(res[1], 3.);
  // EXPECT_DOUBLE_EQ(res[2], 5.);
  // EXPECT_DOUBLE_EQ(res[3], 6.);
}

// TEST(core_vector_serial_eigen_class, substractOperator)
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
 
//   myvec_t res = m_v1 - m_v2;
//   EXPECT_DOUBLE_EQ(res[0], 2.);
//   EXPECT_DOUBLE_EQ(res[1], 1.);
//   EXPECT_DOUBLE_EQ(res[2], 3.);
//   EXPECT_DOUBLE_EQ(res[3], 4.);
// }

// TEST(core_vector_serial_eigen_class, starOperator)
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

//   myvec_t res = m_v1 * m_v2;
//   EXPECT_DOUBLE_EQ(res[0], 3.);
//   EXPECT_DOUBLE_EQ(res[1], 2.);
//   EXPECT_DOUBLE_EQ(res[2], 4.);
//   EXPECT_DOUBLE_EQ(res[3], 5.);
// }

TEST(core_vector_serial_eigen_class, CompoundAssignAddOperator)
{
	using namespace rompp;
  using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  using myvec_t = core::Vector<eigvec_t>;

  myvec_t m_v1(4);
  m_v1[0] = 3.; m_v1[1] = 2.;
  m_v1[2] = 4.; m_v1[3] = 5.;
  myvec_t m_v2(4);
  m_v2[0] = 1.; m_v2[1] = 1.;
  m_v2[2] = 1.; m_v2[3] = 1.;

  m_v1 += m_v2;
  EXPECT_DOUBLE_EQ(m_v1[0], 4.);
  EXPECT_DOUBLE_EQ(m_v1[1], 3.);
  EXPECT_DOUBLE_EQ(m_v1[2], 5.);
  EXPECT_DOUBLE_EQ(m_v1[3], 6.);
}


TEST(core_vector_serial_eigen_class, CompoundAssignSubtractOperator)
{
	using namespace rompp;
  using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  using myvec_t = core::Vector<eigvec_t>;

  myvec_t m_v1(4);
  m_v1[0] = 3.; m_v1[1] = 2.;
  m_v1[2] = 4.; m_v1[3] = 5.;
  myvec_t m_v2(4);
  m_v2[0] = 1.; m_v2[1] = 1.;
  m_v2[2] = 1.; m_v2[3] = 1.;

  m_v1 -= m_v2;
  EXPECT_DOUBLE_EQ(m_v1[0], 2.);
  EXPECT_DOUBLE_EQ(m_v1[1], 1.);
  EXPECT_DOUBLE_EQ(m_v1[2], 3.);
  EXPECT_DOUBLE_EQ(m_v1[3], 4.);
}
