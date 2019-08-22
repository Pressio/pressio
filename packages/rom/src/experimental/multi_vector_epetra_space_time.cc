
#include <gtest/gtest.h>
#include "Epetra_MpiComm.h"
#include "CONTAINERS_ALL"
#include "../../src/experimental/codes.hpp"

using namespace pressio;

struct MultiVectorEpetraSpaceTimeR9C10Fixture
  : public ::testing::Test{

public:
  int rank_;
  // Epetra_MpiComm * Comm_;
  std::shared_ptr<Epetra_MpiComm> Comm_;
  int MyPID_;
  int NumProc_;
  const int localSize_ = 3;
  const int numVectors_ = 10;
  int numGlobalEntries_;
  // Epetra_Map * contigMap_;  
  std::shared_ptr<Epetra_Map> contigMap_;
  std::shared_ptr<Epetra_MultiVector> A_;
  // containers::Matrix<Epetra_CrsMatrix> * A_;

  virtual void SetUp()
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    // Comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
    Comm_ = std::make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    MyPID_ = Comm_->MyPID();
    NumProc_ = Comm_->NumProc();
    numGlobalEntries_ = Comm_->NumProc() * localSize_;

    // contigMap_ = new Epetra_Map(numGlobalEntries_, 0, *Comm_);
    contigMap_ = std::make_shared<Epetra_Map>(numGlobalEntries_, 0, *Comm_);
    A_ = std::make_shared<Epetra_MultiVector>(*contigMap_, numVectors_);

    // A_ = new containers::Matrix<Epetra_CrsMatrix>(*contigMap_, 3);
  }
  
  virtual void TearDown(){
    // delete Comm_;
    // delete contigMap_;
    // delete A_;
  }
};


TEST_F(MultiVectorEpetraSpaceTimeR9C10Fixture, Constructor)
{
  // int myN = contigMap_->NumMyElements();

  // // get my global indices
  // std::vector<int> mygid(myN);
  // contigMap_->MyGlobalElements( mygid.data() );
  
  std::cout << "print A_ " << std::endl;
  A_->Print(std::cout);

  
  // A_.data()->Print(std::cout);

  // std::array<double,3> vals;
  // std::array<int,3> colind;
  
  // for (auto const & it : mygid){
  //   if(it == 0){
  //     vals = {1., 5., 2.};
  //     colind = {0, 3, 4};
  //     A_->insertGlobalValues(it, 3, vals.data(), colind.data());
  //   }
  //   if(it == 3){
  //     vals = {11., 55., 22.};
  //     colind = {6, 8, 11};
  //     A_->insertGlobalValues(it, 3, vals.data(), colind.data());
  //   }
  //   if(it == 11){
  //     vals = {111., 555.};
  //     colind = {3, 12};
  //     A_->insertGlobalValues(it, 2, vals.data(), colind.data());
  //   }
  // }
  // A_->fillingIsCompleted();
  // //A_->data()->Print(std::cout);

  // std::cout << A_->data()->NumGlobalRows() << " " << 
  //              A_->data()->NumGlobalCols() << std::endl;
  
}

