
#ifndef CORE_FIXTURES_TPETRA_ONLY_FIXTURES_HPP_
#define CORE_FIXTURES_TPETRA_ONLY_FIXTURES_HPP_

#include <gtest/gtest.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_CommHelpers.hpp>
#include "CORE_ALL"

struct tpetraVectorGlobSize15Fixture
  : public ::testing::Test{

public:
  using ST = double;
  using LO = int;
  using GO = int;
  using NT = Tpetra::Vector<>::node_type;
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<LO,GO,NT>;
  using vec_t = Tpetra::Vector<ST,LO,GO,NT>;

  int rank_;
  int numProc_;
  const int localSize_ = 5;
  int numGlobalEntries_;
  Teuchos::RCP<const tcomm> comm_;
  Teuchos::RCP<const map_t> contigMap_;
  vec_t * x_;
  
  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    rank_ = comm_->getRank();
    numProc_ = comm_->getSize();
    assert(numProc_==3);
    numGlobalEntries_ = numProc_ * localSize_;
    contigMap_ = Teuchos::rcp(new map_t(numGlobalEntries_, 0,
    					comm_,
    					Tpetra::GloballyDistributed));
    x_ = new vec_t(contigMap_);
  }
  
  virtual void TearDown(){
    delete x_;
  }
};
//-----------------------------------------------------------


struct tpetraMultiVectorGlobSize15Fixture
  : public ::testing::Test{

public:
  using ST = double;
  using LO = int;
  using GO = int;
  using NT = Tpetra::Vector<>::node_type;
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<LO,GO,NT>;
  using mvec_t = Tpetra::MultiVector<ST,LO,GO,NT>;

  int rank_;
  int numProc_;
  const int localSize_ = 5;
  const int numVecs_ = 4;
  int numGlobalEntries_;
  Teuchos::RCP<const tcomm> comm_;
  Teuchos::RCP<const map_t> contigMap_;
  mvec_t * x_;
  
  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    rank_ = comm_->getRank();
    numProc_ = comm_->getSize();
    assert(numProc_==3);
    numGlobalEntries_ = numProc_ * localSize_;
    contigMap_ = Teuchos::rcp(new map_t(numGlobalEntries_, 0,
    					comm_,
    					Tpetra::GloballyDistributed));
    x_ = new mvec_t(contigMap_, numVecs_);
  }
  
  virtual void TearDown(){
    delete x_;
  }
};
//-----------------------------------------------------------


struct tpetraMultiVectorGlobSize9Fixture
  : public ::testing::Test{

public:
  using ST = double;
  using LO = int;
  using GO = int;
  using NT = Tpetra::Vector<>::node_type;
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<LO,GO,NT>;
  using mvec_t = Tpetra::MultiVector<ST,LO,GO,NT>;

  int rank_;
  int numProc_;
  const int localSize_ = 3;
  const int numVecs_ = 4;
  int numGlobalEntries_;
  Teuchos::RCP<const tcomm> comm_;
  Teuchos::RCP<const map_t> contigMap_;
  mvec_t * x_;
  
  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    rank_ = comm_->getRank();
    numProc_ = comm_->getSize();
    assert(numProc_==3);
    numGlobalEntries_ = numProc_ * localSize_;
    contigMap_ = Teuchos::rcp(new map_t(numGlobalEntries_, 0,
    					comm_,
    					Tpetra::GloballyDistributed));
    x_ = new mvec_t(contigMap_, numVecs_);
  }
  
  virtual void TearDown(){
    delete x_;
  }
};
//-----------------------------------------------------------



struct tpetraMultiVectorR9C4VecS9Fixture 
  : public ::testing::Test{

public:
  using ST = double;
  using LO = int;
  using GO = int;
  using NT = Tpetra::Vector<>::node_type;
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<LO,GO,NT>;
  using mvec_t = Tpetra::MultiVector<ST,LO,GO,NT>;
  using vec_t = Tpetra::Vector<ST,LO,GO,NT>;

  int rank_;
  int numProc_;
  const int localSize_ = 3;
  const int numVecs_ = 4;
  int numGlobalEntries_;
  Teuchos::RCP<const tcomm> comm_;
  Teuchos::RCP<const map_t> contigMap_;
  mvec_t * mv_;
  vec_t * x_;
  
  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    rank_ = comm_->getRank();
    numProc_ = comm_->getSize();
    assert(numProc_==3);
    numGlobalEntries_ = numProc_ * localSize_;
    contigMap_ = Teuchos::rcp(new map_t(numGlobalEntries_, 0,
		comm_,Tpetra::GloballyDistributed));
    mv_ = new mvec_t(contigMap_, numVecs_);
    x_ = new vec_t(contigMap_);
  }
  
  virtual void TearDown(){
    delete x_;
    delete mv_;
  }
};
//-----------------------------------------------------------


#endif /* CORE_FIXTURES_EPETRA_ONLY_FIXTURES_HPP_ */

