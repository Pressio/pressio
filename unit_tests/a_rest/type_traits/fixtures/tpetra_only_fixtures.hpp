
#ifndef CONTAINERS_FIXTURES_TPETRA_ONLY_FIXTURES_HPP_
#define CONTAINERS_FIXTURES_TPETRA_ONLY_FIXTURES_HPP_

#include <gtest/gtest.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_CommHelpers.hpp>
#include "pressio_containers.hpp"
#include <Tpetra_Map_decl.hpp>

/* the tpetra data structures below are
 * left without templates such that it picks
 * up the default. So if kokkos is built with
 * CUDA, then tpetra vector, matrices operate on CUDA.
 * IF kokkos is built for openmp, then on-node parallelim
 * in tpetra operations is handled via openmp
 */

struct tpetraVectorGlobSize15Fixture
  : public ::testing::Test{

public:
  using NT = Tpetra::Vector<>::node_type;
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;
  using vec_t = Tpetra::Vector<>;

  using ST = typename vec_t::scalar_type;
  using LO = typename vec_t::local_ordinal_type;
  using GO = typename vec_t::global_ordinal_type;

  int rank_;
  int numProc_;
  const int localSize_ = 5;
  int numGlobalEntries_;
  Teuchos::RCP<const tcomm> comm_;
  Teuchos::RCP<const map_t> contigMap_;
  std::shared_ptr<vec_t> x_;

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    rank_ = comm_->getRank();
    numProc_ = comm_->getSize();
    EXPECT_EQ(numProc_,3);
    numGlobalEntries_ = numProc_ * localSize_;
    contigMap_ = Teuchos::rcp(new map_t(numGlobalEntries_, 0,
    					comm_));
    x_ = std::make_shared<vec_t>(contigMap_);
  }

  virtual void TearDown(){}
};
//-----------------------------------------------------------


struct tpetraMultiVectorGlobSize15Fixture
  : public ::testing::Test{

public:
  using NT = Tpetra::Vector<>::node_type;
  using tcomm = Teuchos::Comm<int>;
  using mvec_t = Tpetra::MultiVector<>;
  using map_t = Tpetra::Map<>;
  using vec_t = Tpetra::Vector<>;
  using ST = typename vec_t::scalar_type;
  using LO = typename vec_t::local_ordinal_type;
  using GO = typename vec_t::global_ordinal_type;

  int rank_;
  int numProc_;
  const int localSize_ = 5;
  const int numVecs_ = 4;
  int numGlobalEntries_;
  Teuchos::RCP<const tcomm> comm_;
  Teuchos::RCP<const map_t> contigMap_;
  std::shared_ptr<mvec_t> x_;

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    rank_ = comm_->getRank();
    numProc_ = comm_->getSize();
    EXPECT_EQ(numProc_,3);
    numGlobalEntries_ = numProc_ * localSize_;
    contigMap_ = Teuchos::rcp(new map_t(numGlobalEntries_, 0,
    					comm_));
    x_ = std::make_shared<mvec_t>(contigMap_, numVecs_);
  }

  virtual void TearDown(){}
};
//-----------------------------------------------------------


struct tpetraMultiVectorGlobSize9Fixture
  : public ::testing::Test{

public:
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;
  using mvec_t = Tpetra::MultiVector<>;
  using ST = typename mvec_t::scalar_type;
  using LO = typename mvec_t::local_ordinal_type;
  using GO = typename mvec_t::global_ordinal_type;

  int rank_;
  int numProc_;
  const int localSize_ = 3;
  const int numVecs_ = 4;
  int numGlobalEntries_;
  Teuchos::RCP<const tcomm> comm_;
  Teuchos::RCP<const map_t> contigMap_;
  std::shared_ptr<mvec_t> x_;

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    rank_ = comm_->getRank();
    numProc_ = comm_->getSize();
    EXPECT_EQ(numProc_,3);
    numGlobalEntries_ = numProc_ * localSize_;
    contigMap_ = Teuchos::rcp(new map_t(numGlobalEntries_, 0, comm_));
    x_ = std::make_shared<mvec_t>(contigMap_, numVecs_);
  }

  virtual void TearDown(){}
};
//-----------------------------------------------------------


struct tpetraMultiVectorR9C4VecS9Fixture
  : public ::testing::Test{

public:
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;
  using mvec_t = Tpetra::MultiVector<>;
  using vec_t = Tpetra::Vector<>;
  using ST = typename mvec_t::scalar_type;
  using LO = typename mvec_t::local_ordinal_type;
  using GO = typename mvec_t::global_ordinal_type;

  int rank_;
  int numProc_;
  const int localSize_ = 3;
  const int numVecs_ = 4;
  int numGlobalEntries_;
  Teuchos::RCP<const tcomm> comm_;
  Teuchos::RCP<const map_t> contigMap_;
  std::shared_ptr<mvec_t> mv_;
  std::shared_ptr<vec_t> x_;

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    rank_ = comm_->getRank();
    numProc_ = comm_->getSize();
    EXPECT_EQ(numProc_,3);
    numGlobalEntries_ = numProc_ * localSize_;
    contigMap_ = Teuchos::rcp(new map_t(numGlobalEntries_,0,comm_));
    mv_ = std::make_shared<mvec_t>(contigMap_, numVecs_);
    x_ = std::make_shared<vec_t>(contigMap_);
  }

  virtual void TearDown(){}
};
//-----------------------------------------------------------


struct tpetraSparseMatR7MultiVectorR7C4Fixture
  : public ::testing::Test{

public:
  using tcomm = Teuchos::Comm<int>;
  using mat_t = Tpetra::CrsMatrix<>;
  using mvec_t = Tpetra::MultiVector<>;
  using ST = typename mvec_t::scalar_type;
  using LO = typename mvec_t::local_ordinal_type;
  using GO = typename mvec_t::global_ordinal_type;
  using map_t = Tpetra::Map<>;

  int rank_;
  int numProc_;
  const int localSize_ = 3;
  int numGlobalEntries_;
  Teuchos::RCP<const tcomm> comm_;
  Teuchos::RCP<const map_t> contigMap_;

  // for crs matrix
  std::shared_ptr<mat_t> A_;

  // for multivector
  const int numVecs_ = 4;
  std::shared_ptr<mvec_t> mv_;

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    rank_ = comm_->getRank();
    numProc_ = comm_->getSize();
    EXPECT_EQ(numProc_,3);
    numGlobalEntries_ = numProc_ * localSize_;
    contigMap_ = Teuchos::rcp(new map_t(numGlobalEntries_,0,comm_));

    A_ = std::make_shared<mat_t>(contigMap_, 4);
    mv_ = std::make_shared<mvec_t>(contigMap_, numVecs_);
  }

  void fillCrsMatrix(){

    // 1 0 2 0 0 3 0
    // 1 0 1 0 1 1 0
    // 0 0 1 2 1 0 0
    //--------------
    // 1 0 0 0 2 1 3
    // 0 0 1 1 0 0 0
    //--------------
    // 0 0 3 3 0 1 0
    // 0 0 0 3 0 1 0

    using tarr_dt = Teuchos::ArrayView<ST>;
    using tarr_it = Teuchos::ArrayView<GO>;
    std::array<ST,4> va;
    std::array<GO,4> ci;
    if (rank_==0){
      //row 0
      va = {{1., 2., 3.}};  ci = {{0, 2, 5}};
      A_->insertGlobalValues( 0, tarr_it(ci.data(), 3), tarr_dt(va.data(), 3) );
      //row 1
      va = {{1.,1.,1.,1.}};  ci = {{0,2,4,5}};
      A_->insertGlobalValues(1, tarr_it(ci.data(), 4), tarr_dt(va.data(), 4) );
      //row 2
      va = {{1.,2.,1.}};  ci = {{2,3,4}};
      A_->insertGlobalValues(2, tarr_it(ci.data(), 3), tarr_dt(va.data(), 3) );
    }

    if (rank_==1){
      //row 3
      va = {{1.,2.,1,3.}};  ci = {{0,4,5,6}};
      A_->insertGlobalValues(3, tarr_it(ci.data(), 4), tarr_dt(va.data(), 4) );
      //row 4
      va = {{1.,1.}};  ci = {{2,3}};
      A_->insertGlobalValues(4, tarr_it(ci.data(), 3), tarr_dt(va.data(), 3) );
    }

    if (rank_==2){
      //row 5
      va = {{3.,3.,1}};  ci = {{2,3,5}};
      A_->insertGlobalValues(5, tarr_it(ci.data(), 3), tarr_dt(va.data(), 3) );
      //row 6
      va = {{3,1.}};  ci = {{3,5}};
      A_->insertGlobalValues(6, tarr_it(ci.data(), 3), tarr_dt(va.data(), 3) );
    }
  }//end fill

  virtual void TearDown(){}
};

#endif
