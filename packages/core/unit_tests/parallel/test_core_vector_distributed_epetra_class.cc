
#include <gtest/gtest.h>
#include "vector/core_vector_meta.hpp"
#include "vector/core_vector_distributed_epetra.hpp"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Version.h"
#include "Epetra_MpiComm.h"



struct core_vector_distributed_epetraFix
  : public ::testing::Test{
public:

  int rank_;
  Epetra_MpiComm * Comm_;
  int MyPID_;
  int NumProc_;
  const int localSize_ = 5;
  int numGlobalEntries_;
  Epetra_Map * contigMap_;
  Epetra_Vector * x_;
  
  virtual void SetUp()
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    Comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
    MyPID_ = Comm_->MyPID();
    NumProc_ = Comm_->NumProc();
    numGlobalEntries_ = Comm_->NumProc() * localSize_;
    contigMap_ = new Epetra_Map(numGlobalEntries_, 0, *Comm_);
    x_ = new Epetra_Vector(*contigMap_);
  }

  int getNumProc() const{
    return NumProc_;
  }
  int numGlobalEntries() const{
    return numGlobalEntries_;
  }
  int numLocalEntries() const{
    return localSize_;
  }
  const Epetra_Map * getMap() const{
    return contigMap_;
  }    
  const Epetra_Vector * getVector() const{
    return x_;
  }
  void fillScalar(double value){
    x_->PutScalar(value);
  }
    
  virtual void TearDown(){
    delete Comm_;
    delete contigMap_;
    delete x_;
  }
};

using core_vector_distributed_epetra_DeathTest = core_vector_distributed_epetraFix;


TEST_F(core_vector_distributed_epetraFix, EpetraVectorConstructor)
{
  using myvec_t = core::vector<Epetra_Vector>;
  ASSERT_TRUE( core::details::traits<myvec_t>::isEpetra==1 );

  myvec_t a( *getMap() );
  ASSERT_EQ( a.globalSize(), numGlobalEntries() );
  ASSERT_EQ( a.localSize(), numLocalEntries() );
  myvec_t a2( *getVector() );
  ASSERT_EQ( a2.globalSize(), numGlobalEntries() );
  ASSERT_EQ( a2.localSize(), numLocalEntries() );
}


TEST_F(core_vector_distributed_epetraFix, EpetraVectorQueryWrappedData)
{    
  using myvec_t = core::vector<Epetra_Vector>;
  myvec_t v1( *getMap() );
  ::testing::StaticAssertTypeEq<decltype(v1.data()),
  				Epetra_Vector * >(); 
  const myvec_t v2( *getVector() );
  ::testing::StaticAssertTypeEq< decltype(v2.data()),
  				 const Epetra_Vector * >();
}

TEST_F(core_vector_distributed_epetraFix, EpetraVectorSubscriptOperator)
{
  using myvec_t = core::vector<Epetra_Vector>;

  fillScalar(11.2);
  myvec_t v1( *getMap() );
  for (size_t i=0; i<v1.localSize(); i++){
    v1[i] = 11.2;
  }
  for (size_t i=0; i<v1.localSize(); i++){
    EXPECT_DOUBLE_EQ( v1[i], (*getVector())[i] );
  }
  v1[3] = 56.;
  EXPECT_DOUBLE_EQ( v1[3], 56.0);
}

TEST_F(core_vector_distributed_epetra_DeathTest, EpetraVectorSubscriptOperator)
{
  using myvec_t = core::vector<Epetra_Vector>;
  myvec_t v1( *getMap() );
  int localSize = numLocalEntries();
  ASSERT_DEATH(v1[localSize+1]=4.0, "Assertion failed:");
}




// TEST(core_vector_serial_eigen, EigenVectorAdditionOperator)
// {
//   using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
//   using myvec_t = core::vector<eigvec_t>;
//   using vecTrait = core::details::traits<myvec_t>;
//   ASSERT_EQ(vecTrait::isEigen, 1);

//   myvec_t m_v1(4);
//   m_v1[0] = 3.; m_v1[1] = 2.;
//   m_v1[2] = 4.; m_v1[3] = 5.;
//   myvec_t m_v2(4);
//   m_v2[0] = 1.; m_v2[1] = 1.;
//   m_v2[2] = 1.; m_v2[3] = 1.;

//   myvec_t res = m_v1 + m_v2;
//   EXPECT_DOUBLE_EQ(res[0], 4.);
//   EXPECT_DOUBLE_EQ(res[1], 3.);
//   EXPECT_DOUBLE_EQ(res[2], 5.);
//   EXPECT_DOUBLE_EQ(res[3], 6.);
// }

// TEST(core_vector_serial_eigen, EigenVectorSubstractOperator)
// {
//   using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
//   using myvec_t = core::vector<eigvec_t>;
//   using vecTrait = core::details::traits<myvec_t>;
//   ASSERT_EQ(vecTrait::isEigen, 1);

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

// TEST(core_vector_serial_eigen, EigenVectorStarOperator)
// {
//   using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
//   using myvec_t = core::vector<eigvec_t>;
//   using vecTrait = core::details::traits<myvec_t>;
//   ASSERT_EQ(vecTrait::isEigen, 1);

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

// TEST(core_vector_serial_eigen, EigenVectorCompoundAssignAddOperator)
// {
//   using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
//   using myvec_t = core::vector<eigvec_t>;
//   using vecTrait = core::details::traits<myvec_t>;
//   ASSERT_EQ(vecTrait::isEigen, 1);

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


// TEST(core_vector_serial_eigen, EigenVectorCompoundAssignSubtractOperator)
// {
//   using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
//   using myvec_t = core::vector<eigvec_t>;
//   using vecTrait = core::details::traits<myvec_t>;
//   ASSERT_EQ(vecTrait::isEigen, 1);

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
