
#ifndef CONTAINERS_FIXTURES_BLOCK_TPETRA_FIXTURES_HPP_
#define CONTAINERS_FIXTURES_BLOCK_TPETRA_FIXTURES_HPP_

#include <gtest/gtest.h>
#include <Tpetra_BlockVector.hpp>
#include <Tpetra_BlockMultiVector.hpp>
#include <Tpetra_Map.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_Map_decl.hpp>

#include "pressio/ops.hpp"

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
  using vec_t = Tpetra::BlockVector<>;

  using ST = typename vec_t::scalar_type;
  using LO = typename vec_t::local_ordinal_type;
  using GO = typename vec_t::global_ordinal_type;
  using NT = typename vec_t::node_type;

  int rank_;
  int numProc_;
  const int localSize_ = 5;
  const int blockSize_ = 5;
  int numGlobalEntries_;
  Teuchos::RCP<const tcomm> comm_;
  Teuchos::RCP<const map_t> contigMap_;
  std::shared_ptr<vec_t> myVector_;

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    rank_ = comm_->getRank();
    numProc_ = comm_->getSize();
    numGlobalEntries_ = localSize_*numProc_;
    contigMap_ = Teuchos::rcp(new map_t(numGlobalEntries_, 0, comm_));
    myVector_ = std::make_shared<vec_t>(*contigMap_, blockSize_);
  }

  virtual void TearDown(){}
};

struct tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture
  : public ::testing::Test
{

public:
  using tcomm  = Teuchos::Comm<int>;
  using map_t  = Tpetra::Map<>;
  using vec_t  = Tpetra::BlockVector<>;
  using mvec_t = Tpetra::BlockMultiVector<>;

  using ST = typename mvec_t::scalar_type;
  using LO = typename mvec_t::local_ordinal_type;
  using GO = typename mvec_t::global_ordinal_type;
  using NT = typename mvec_t::node_type;

  int rank_;
  int numProc_;
  const int localSize_ = 5;
  const int blockSize_ = 4;
  const int numVecs_ = 3;
  int numGlobalEntries_;
  Teuchos::RCP<const tcomm> comm_;
  Teuchos::RCP<const map_t> contigMap_;

  // level2
  std::shared_ptr<mvec_t> myMv_;
  std::shared_ptr<vec_t> x_tpetra;
  std::shared_ptr<vec_t> y_tpetra;

  // level3
  decltype(::pressio::ops::clone(*myMv_)) A;
  std::shared_ptr<mvec_t> B;
  static std::array<double, 4> ac;
  static std::array<double, 3> bc;
  Eigen::MatrixXd C_eigen;
  Kokkos::View<double**, Kokkos::LayoutLeft> C_kokkos;

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    rank_ = comm_->getRank();
    numProc_ = comm_->getSize();

    numGlobalEntries_ = localSize_*numProc_;
    contigMap_ = Teuchos::rcp(new map_t(numGlobalEntries_, 0, comm_));

    // initialize data for computations
    myMv_ = std::make_shared<mvec_t>(*contigMap_, blockSize_, numVecs_);
    auto myMv_h = myMv_->getMultiVectorView().getLocalViewHost(Tpetra::Access::ReadWrite);
    for (int i = 0; i < localSize_ * blockSize_; ++i){
      for (int j = 0; j < numVecs_; ++j){
        // generate rank-unique int values
        myMv_h(i, j) = (double)((rank_ * localSize_ * blockSize_ + i) * numVecs_ + j + 1.);
        printf("\t@ [%d](%d, %d) -> %g\n", rank_, i, j, myMv_h(i, j));
      }
    }
    x_tpetra = std::make_shared<vec_t>(*contigMap_, blockSize_);
    auto x_tpetra_h = x_tpetra->getVectorView().getLocalViewHost(Tpetra::Access::ReadWrite);
    for (int j = 0; j < localSize_; ++j) {
      // generate rank-unique int values
      x_tpetra_h(j, 0) = rank_ * localSize_ + (double)(j + 1.);
    }
    y_tpetra = std::make_shared<vec_t>(*contigMap_, blockSize_);

    // level3 data
    A = pressio::ops::clone(*myMv_);
    for (decltype(A.getNumVectors()) i=0; i<A.getNumVectors(); ++i) {
      A.getMultiVectorView().getVectorNonConst(i)->putScalar(ac[i]);
    }
    B = std::make_shared<mvec_t>(*contigMap_, blockSize_, bc.size());
    for (int i=0; i<(int)bc.size(); ++i) {
      B->getMultiVectorView().getVectorNonConst(i)->putScalar(bc[i]);
    }
    C_eigen = Eigen::MatrixXd(A.getNumVectors(), A.getNumVectors());
    C_eigen.setConstant(0.);
    C_kokkos = Kokkos::View<double**, Kokkos::LayoutLeft>("C", A.getNumVectors(), A.getNumVectors());
    Kokkos::deep_copy(C_kokkos, 0.);
  }

  virtual void TearDown(){}

  auto gold_a(size_t i, size_t j) {
      return ac[i]*A.getMultiVectorView().getGlobalLength()*1.5*ac[j];
  }

  auto gold_ab(size_t i, size_t j) {
      return ac[i]*A.getMultiVectorView().getGlobalLength()*1.5*bc[j];
  }
};

std::array<double, 4> tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture::ac{1.,2.,3.,4.};
std::array<double, 3> tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture::bc{1.2, 2.2, 3.2};

#endif
