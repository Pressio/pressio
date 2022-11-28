
#include "tpetra_only_fixtures.hpp"
#include "pressio/ops.hpp"

TEST_F(tpetraMultiVectorGlobSize15Fixture, mv_prod_kokkos_vector)
{
    auto myMv_h = myMv_->getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
        myMv_h(i,j) = (double)j;
     }
    }

    Kokkos::View<double*> a("a", numVecs_);
    KokkosBlas::fill(a, 1.);

    vec_t y(contigMap_);
    y.putScalar(0.);
    pressio::ops::product(::pressio::nontranspose{}, 1., *myMv_, a, 1., y);

    auto y_h = y.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
    EXPECT_DOUBLE_EQ(y_h(0,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(1,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(2,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(3,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(4,0), 6.);

    // simulate beta=0 with uninitialized y containing NaN
    const auto nan = std::nan("0");
    y.putScalar(nan);
    pressio::ops::product(::pressio::nontranspose{}, 1., *myMv_, a, 0., y);

    y_h = y.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
    EXPECT_DOUBLE_EQ(y_h(0,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(1,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(2,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(3,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(4,0), 6.);

    // simulate alpha=0 with NaN in input matrix
    myMv_h(0,0) = nan;
    pressio::ops::product(::pressio::nontranspose{}, 0., *myMv_, a, 1., y);
    y_h = y.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
    EXPECT_DOUBLE_EQ(y_h(0,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(1,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(2,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(3,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(4,0), 6.);

    // alpha == beta == 0
    pressio::ops::product(::pressio::nontranspose{}, 0., *myMv_, a, 0., y);
    y_h = y.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
    EXPECT_DOUBLE_EQ(y_h(0,0), 0.);
    EXPECT_DOUBLE_EQ(y_h(1,0), 0.);
    EXPECT_DOUBLE_EQ(y_h(2,0), 0.);
    EXPECT_DOUBLE_EQ(y_h(3,0), 0.);
    EXPECT_DOUBLE_EQ(y_h(4,0), 0.);
}

TEST_F(tpetraMultiVectorGlobSize15Fixture, mv_T_vector_storein_kokkos_vector)
{
    auto myMv_h = myMv_->getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
        myMv_h(i,j) = (double)i;
     }
    }

    vec_t y(contigMap_);
    y.putScalar(2.);

    Kokkos::View<double*> a("a", numVecs_);
    KokkosBlas::fill(a, 1.);
    pressio::ops::product(::pressio::transpose{}, 1., *myMv_, y, 1., a);

    auto a_h = Kokkos::create_mirror_view(a);
    Kokkos::deep_copy(a_h, a);
    const auto a1_ref = numProc_ * 2. * 10. + 1.;
    EXPECT_DOUBLE_EQ(a_h(0), a1_ref);
    EXPECT_DOUBLE_EQ(a_h(1), a1_ref);
    EXPECT_DOUBLE_EQ(a_h(2), a1_ref);
    EXPECT_DOUBLE_EQ(a_h(3), a1_ref);

    // simulate beta=0 with uninitialized y containing NaN
    const auto nan = std::nan("0");
    KokkosBlas::fill(a, nan);
    pressio::ops::product(::pressio::transpose{}, 1., *myMv_, y, 0., a);

    a_h = Kokkos::create_mirror_view(a);
    Kokkos::deep_copy(a_h, a);
    const auto a0_ref = numProc_ * 2. * 10. + 0.;
    EXPECT_DOUBLE_EQ(a_h(0), a0_ref);
    EXPECT_DOUBLE_EQ(a_h(1), a0_ref);
    EXPECT_DOUBLE_EQ(a_h(2), a0_ref);
    EXPECT_DOUBLE_EQ(a_h(3), a0_ref);

    // simulate alpha=0 with NaN in input matrix
    myMv_h(0,0) = nan;
    pressio::ops::product(::pressio::transpose{}, 0., *myMv_, y, 1., a);
    Kokkos::deep_copy(a_h, a);
    EXPECT_DOUBLE_EQ(a_h(0), a0_ref);
    EXPECT_DOUBLE_EQ(a_h(1), a0_ref);
    EXPECT_DOUBLE_EQ(a_h(2), a0_ref);
    EXPECT_DOUBLE_EQ(a_h(3), a0_ref);

    // alpha == beta == 0
    pressio::ops::product(::pressio::transpose{}, 0., *myMv_, y, 0., a);
    Kokkos::deep_copy(a_h, a);
    EXPECT_DOUBLE_EQ(a_h(0), 0.);
    EXPECT_DOUBLE_EQ(a_h(1), 0.);
    EXPECT_DOUBLE_EQ(a_h(2), 0.);
    EXPECT_DOUBLE_EQ(a_h(3), 0.);
}

#ifdef PRESSIO_ENABLE_TPL_EIGEN

template <typename XType>
void test_with_eigen_x(const tpetraMultiVectorGlobSize15Fixture &test, XType x) {
    auto A = *test.myMv_;
    auto A_h = A.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i = 0; i < test.localSize_; ++i) {
      for (int j = 0; j < test.numVecs_; ++j) {
        A_h(i, j) = (double)j;
      }
    }

    tpetraMultiVectorGlobSize15Fixture::vec_t y(test.contigMap_);
    y.putScalar(0.);
    pressio::ops::product(::pressio::nontranspose{}, 1., A, x, 1., y);

    auto y_h = y.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
    EXPECT_DOUBLE_EQ(y_h(0,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(1,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(2,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(3,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(4,0), 6.);

    // simulate beta=0 with uninitialized y containing NaN
    const auto nan = std::nan("0");
    y.putScalar(nan);
    pressio::ops::product(::pressio::nontranspose{}, 1., A, x, 0., y);

    y_h = y.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
    EXPECT_DOUBLE_EQ(y_h(0,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(1,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(2,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(3,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(4,0), 6.);

    // simulate alpha=0 with NaN in input matrix
    A_h(0, 0) = nan;
    pressio::ops::product(::pressio::nontranspose{}, 0., A, x, 1., y);
    y_h = y.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
    EXPECT_DOUBLE_EQ(y_h(0,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(1,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(2,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(3,0), 6.);
    EXPECT_DOUBLE_EQ(y_h(4,0), 6.);

    // alpha == beta == 0
    pressio::ops::product(::pressio::nontranspose{}, 0., A, x, 0., y);
    y_h = y.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
    EXPECT_DOUBLE_EQ(y_h(0,0), 0.);
    EXPECT_DOUBLE_EQ(y_h(1,0), 0.);
    EXPECT_DOUBLE_EQ(y_h(2,0), 0.);
    EXPECT_DOUBLE_EQ(y_h(3,0), 0.);
    EXPECT_DOUBLE_EQ(y_h(4,0), 0.);
}

TEST_F(tpetraMultiVectorGlobSize15Fixture, mv_prod_eigen_vector)
{
    Eigen::VectorXd x0(numVecs_);
    x0.setConstant(1.);
    test_with_eigen_x(*this, x0);
}

TEST_F(tpetraMultiVectorGlobSize15Fixture, mv_prod_eigen_span)
{
    Eigen::VectorXd x0(numVecs_ + 3);
    x0.setConstant(1.);
    test_with_eigen_x(*this, pressio::span(x0, 2, numVecs_));
}

TEST_F(tpetraMultiVectorGlobSize15Fixture, mv_prod_eigen_diag)
{
    Eigen::MatrixXd M0(numVecs_, numVecs_);
    for (int i = 0; i < numVecs_; ++i) {
        M0(i, i) = 1.;
    }
    test_with_eigen_x(*this, pressio::diag(M0));
}

template <typename YType>
void test_with_eigen_y(const tpetraMultiVectorGlobSize15Fixture &test, YType y) {
    auto A = *test.myMv_;
    auto A_h = A.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i = 0; i < test.localSize_; ++i) {
      for (int j = 0; j < test.numVecs_; ++j) {
        A_h(i, j) = (double)i;
      }
    }

    tpetraMultiVectorGlobSize15Fixture::vec_t x(test.contigMap_);
    x.putScalar(2.);

    pressio::ops::product(::pressio::transpose{}, 1., A, x, 1., y);

    const auto a1_ref = test.numProc_ * 2. * 10. + 1.;
    EXPECT_DOUBLE_EQ(y(0), a1_ref);
    EXPECT_DOUBLE_EQ(y(1), a1_ref);
    EXPECT_DOUBLE_EQ(y(2), a1_ref);
    EXPECT_DOUBLE_EQ(y(3), a1_ref);

    // simulate beta=0 with uninitialized y containing NaN
    const auto nan = std::nan("0");
    ::pressio::ops::fill(y, nan);
    pressio::ops::product(::pressio::transpose{}, 1., A, x, 0., y);

    const auto a0_ref = test.numProc_ * 2. * 10. + 0.;
    EXPECT_DOUBLE_EQ(y(0), a0_ref);
    EXPECT_DOUBLE_EQ(y(1), a0_ref);
    EXPECT_DOUBLE_EQ(y(2), a0_ref);
    EXPECT_DOUBLE_EQ(y(3), a0_ref);

    // simulate alpha=0 with NaN in input matrix
    A_h(0, 0) = nan;
    pressio::ops::product(::pressio::transpose{}, 0., A, x, 1., y);
    EXPECT_DOUBLE_EQ(y(0), a0_ref);
    EXPECT_DOUBLE_EQ(y(1), a0_ref);
    EXPECT_DOUBLE_EQ(y(2), a0_ref);
    EXPECT_DOUBLE_EQ(y(3), a0_ref);

    // alpha == beta == 0
    pressio::ops::product(::pressio::transpose{}, 0., A, x, 0., y);
    EXPECT_DOUBLE_EQ(y(0), 0.);
    EXPECT_DOUBLE_EQ(y(1), 0.);
    EXPECT_DOUBLE_EQ(y(2), 0.);
    EXPECT_DOUBLE_EQ(y(3), 0.);
}

TEST_F(tpetraMultiVectorGlobSize15Fixture, mv_T_vector_storein_eigen_vector)
{
    Eigen::VectorXd y0(numVecs_);
    y0.setConstant(1.);
    test_with_eigen_y(*this, y0);
}

TEST_F(tpetraMultiVectorGlobSize15Fixture, mv_T_vector_storein_eigen_span)
{
    Eigen::VectorXd y0(numVecs_ + 3);
    y0.setConstant(1.);
    test_with_eigen_y(*this, pressio::span(y0, 2, numVecs_));
}

TEST_F(tpetraMultiVectorGlobSize15Fixture, mv_T_vector_storein_eigen_diag)
{
    Eigen::MatrixXd M0(numVecs_, numVecs_);
    for (int i = 0; i < numVecs_; ++i) {
        M0(i, i) = 1.;
    }
    test_with_eigen_y(*this, pressio::diag(M0));
}
#endif
