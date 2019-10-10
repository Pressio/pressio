
#include "epetra_only_fixtures.hpp"
#include "Epetra_MpiComm.h"
#include "CONTAINERS_ALL"


TEST_F(epetraSparseMatR7MultiVectorR9C4Fixture,
      CrsMatrixWrapperConstructor){
  
  sm_->fillingIsCompleted();
  using mymat_w_t = pressio::containers::Matrix<Epetra_CrsMatrix>;
  mymat_w_t Bw(*sm_);
}

TEST_F(epetraSparseMatR7MultiVectorR9C4Fixture,
      rescaling){

  fillCrsMatrix();
  sm_->fillingIsCompleted();
  sm_->data()->Print(std::cout);

  
  sm_->addToDiagonal(-2.);
  Epetra_Vector dd(*dataMapSM_);
  sm_->data()->ExtractDiagonalCopy(dd);
  dd.Print(std::cout);

  if(rank_==0){
    EXPECT_DOUBLE_EQ( dd[0], -1.);
    EXPECT_DOUBLE_EQ( dd[1], 0.);
    EXPECT_DOUBLE_EQ( dd[2], -1.);
  }
  if(rank_==1){
    EXPECT_DOUBLE_EQ( dd[0], 0.);
    EXPECT_DOUBLE_EQ( dd[1], 0.);
  }
  if(rank_==2){
    EXPECT_DOUBLE_EQ( dd[0], -1.);
    EXPECT_DOUBLE_EQ( dd[1], 0.);
  }
}

// WIP 
