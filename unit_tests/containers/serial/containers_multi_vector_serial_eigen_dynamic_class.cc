
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"

using eigdmat_t = Eigen::MatrixXd;
using myMV_t = pressio::containers::MultiVector<eigdmat_t>;


TEST(containers_multi_vector_serial_eigen_dynamic_class,
     constructor){

  using MVTrait = pressio::containers::details::traits<myMV_t>;
  ASSERT_TRUE(MVTrait::wrapped_multi_vector_identifier
  == pressio::containers::details::WrappedMultiVectorIdentifier::Eigen);

  ASSERT_TRUE(
  pressio::containers::meta::is_multi_vector_wrapper_eigen<myMV_t>::value);
  
  //construct by passing the size 
  myMV_t A(6,3);

  // pass native eigen data
  eigdmat_t eA(45,32);
  myMV_t AW(eA);
}


TEST(containers_multi_vector_serial_eigen_dynamic_class,
     constructorAndCheckVals){

  //construct by passing the sizes 
  myMV_t A(6,3);
  ASSERT_TRUE( A.extent(0) == 6 );
  ASSERT_TRUE( A.numVectors() == 3 );
  for (size_t i=0; i<6; i++)
    for (size_t j=0; j<3; j++)
      EXPECT_DOUBLE_EQ( A(i,j), 0.);
  
  // pass native eigen vector
  eigdmat_t eA(45,12);
  eA(2,2) = 2.2;
  eA(4,11) = 4.4;

  myMV_t B(eA);
  ASSERT_TRUE( B.extent(0) == 45 );
  ASSERT_TRUE( B.numVectors() == 12 );
  ASSERT_FALSE( B.extent(0) == 4 );
  ASSERT_FALSE( B.numVectors() == 1 );
  EXPECT_DOUBLE_EQ( B(2,2), 2.2);
  EXPECT_DOUBLE_EQ( B(4,11), 4.4);
}


TEST(containers_multi_vector_serial_eigen_dynamic_class,
     copyConstructor){

  //construct by passing the size 
  myMV_t a(3,2);
  a(0,0)=1.1;
  a(0,1)=1.2;
  a(2,1)=1.3;
  
  myMV_t b(a);
  EXPECT_DOUBLE_EQ( b(0,0), 1.1);
  EXPECT_DOUBLE_EQ( b(0,1), 1.2);
  EXPECT_DOUBLE_EQ( b(2,1), 1.3);
}
