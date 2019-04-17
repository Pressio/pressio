
#ifndef CORE_FIXTURES_BLOCK_TPETRA_FIXTURES_HPP_
#define CORE_FIXTURES_BLOCK_TPETRA_FIXTURES_HPP_

#include <gtest/gtest.h>
#include "CORE_ALL"
#include <Tpetra_Experimental_BlockVector_decl.hpp>
#include <Tpetra_Experimental_BlockMultiVector_decl.hpp>
#include <Tpetra_Experimental_BlockVector_def.hpp>
#include <Tpetra_Experimental_BlockVector.hpp>
#include <Tpetra_Map.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_Map_decl.hpp>



/* the tpetra data structures below are
 * left without templates such that it picks
 * up the default. So if kokkos is built with
 * CUDA, then tpetra vector, matrices operate on CUDA.
 * IF kokkos is built for openmp, then on-node parallelim
 * in tpetra operations is handled via openmp
 */

struct tpetraBlockVectorGlobSize15BlockSize5Fixture
  : public ::testing::Test{

public:
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;
  using vec_t = Tpetra::Experimental::BlockVector<>;

  using ST = typename vec_t::scalar_type;
  using LO = typename vec_t::local_ordinal_type;
  using GO = typename vec_t::global_ordinal_type;
  using NT = typename vec_t::node_type;

  int rank_;
  int numProc_;
  const int blockSize_ = 4;
  int numGlobalEntries_;
  Teuchos::RCP<const tcomm> comm_;
  Teuchos::RCP<const map_t> contigMap_;
  std::shared_ptr<vec_t> x_;

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    rank_ = comm_->getRank();
    numProc_ = comm_->getSize();

    numGlobalEntries_ = 15;
    contigMap_ = Teuchos::rcp(new map_t(numGlobalEntries_, 0, comm_));
    x_ = std::make_shared<vec_t>(*contigMap_, blockSize_);
  }

  virtual void TearDown(){}
};

#endif
