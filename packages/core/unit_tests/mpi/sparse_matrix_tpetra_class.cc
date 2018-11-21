
#include "tpetra_only_fixtures.hpp"
#include "CORE_ALL"

TEST_F(tpetraSparseMatR7MultiVectorR7C4Fixture,
      CrsMatrixWrapperConstructor){

  using namespace rompp;

  using mat_t = typename tpetraSparseMatR7MultiVectorR7C4Fixture::mat_t;
  STATIC_ASSERT_IS_MATRIX_SPARSE_DISTRIBUTED_TPETRA(mat_t);
  using mymat_w_t = core::Matrix<mat_t>;
  static_assert(
   core::meta::is_tpetra_sparse_matrix_wrapper<mymat_w_t>::value, " ");

  fillCrsMatrix();
  A_->fillComplete();
  // A_->print(std::cout);
  // //  A_->fillingIsCompleted();
  mymat_w_t Aw(contigMap_, 4);
}

// WIP
