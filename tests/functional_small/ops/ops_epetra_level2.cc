
#include "epetra_only_fixtures.hpp"
#include "pressio/ops.hpp"

TEST_F(epetraMultiVectorGlobSize15Fixture, mv_prod_teuchos_vector)
{
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
        (*myMv_)[j][i] = (double)j;
     }
    }

    Teuchos::SerialDenseVector<int, double> a(numVecs_);
    a = 1.;

    Epetra_Vector y(*contigMap_);
    y.PutScalar(0.);
    pressio::ops::product(::pressio::nontranspose{}, 1., *myMv_, a, 1., y);

    EXPECT_DOUBLE_EQ(y[0], 6.);
    EXPECT_DOUBLE_EQ(y[1], 6.);
    EXPECT_DOUBLE_EQ(y[2], 6.);
    EXPECT_DOUBLE_EQ(y[3], 6.);

    // simulate beta=0 with uninitialized y containing NaN
    const auto nan = std::nan("0");
    y.PutScalar(nan);
    pressio::ops::product(::pressio::nontranspose{}, 1., *myMv_, a, 0., y);

    EXPECT_DOUBLE_EQ(y[0], 6.);
    EXPECT_DOUBLE_EQ(y[1], 6.);
    EXPECT_DOUBLE_EQ(y[2], 6.);
    EXPECT_DOUBLE_EQ(y[3], 6.);

    // simulate alpha=0 with NaN in input matrix
    (*myMv_)[0][0] = nan;
    pressio::ops::product(::pressio::nontranspose{}, 0., *myMv_, a, 1., y);

    EXPECT_DOUBLE_EQ(y[0], 6.);
    EXPECT_DOUBLE_EQ(y[1], 6.);
    EXPECT_DOUBLE_EQ(y[2], 6.);
    EXPECT_DOUBLE_EQ(y[3], 6.);

    // alpha == beta == 0
    pressio::ops::product(::pressio::nontranspose{}, 0., *myMv_, a, 0., y);

    EXPECT_DOUBLE_EQ(y[0], 0.);
    EXPECT_DOUBLE_EQ(y[1], 0.);
    EXPECT_DOUBLE_EQ(y[2], 0.);
    EXPECT_DOUBLE_EQ(y[3], 0.);
}

TEST_F(epetraMultiVectorGlobSize15Fixture, mv_prod_storein_teuchos_vector)
{
    auto & myMv_h = *myMv_;
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
        myMv_h[j][i] = (double)i;
     }
    }

    Epetra_Vector y(*contigMap_);
    y.PutScalar(2.);

    Teuchos::SerialDenseVector<int, double> a(numVecs_);
    a = 1.;
    pressio::ops::product(::pressio::transpose{}, 1., *myMv_, y, 1., a);

    const auto a1_ref = numProc_ * 2. * 10. + 1.;
    EXPECT_DOUBLE_EQ(a(0), a1_ref);
    EXPECT_DOUBLE_EQ(a(1), a1_ref);
    EXPECT_DOUBLE_EQ(a(2), a1_ref);
    EXPECT_DOUBLE_EQ(a(3), a1_ref);

    // simulate beta=0 with uninitialized y containing NaN
    const auto nan = std::nan("0");
    a = nan;
    pressio::ops::product(::pressio::transpose{}, 1., *myMv_, y, 0., a);

    const auto a0_ref = numProc_ * 2. * 10. + 0.;
    EXPECT_DOUBLE_EQ(a(0), a0_ref);
    EXPECT_DOUBLE_EQ(a(1), a0_ref);
    EXPECT_DOUBLE_EQ(a(2), a0_ref);
    EXPECT_DOUBLE_EQ(a(3), a0_ref);

    // simulate alpha=0 with NaN in input matrix
    myMv_h[0][0] = nan;
    pressio::ops::product(::pressio::transpose{}, 0., *myMv_, y, 1., a);

    EXPECT_DOUBLE_EQ(a(0), a0_ref);
    EXPECT_DOUBLE_EQ(a(1), a0_ref);
    EXPECT_DOUBLE_EQ(a(2), a0_ref);
    EXPECT_DOUBLE_EQ(a(3), a0_ref);

    // alpha == beta == 0
    pressio::ops::product(::pressio::transpose{}, 0., *myMv_, y, 0., a);

    EXPECT_DOUBLE_EQ(a(0), 0.);
    EXPECT_DOUBLE_EQ(a(1), 0.);
    EXPECT_DOUBLE_EQ(a(2), 0.);
    EXPECT_DOUBLE_EQ(a(3), 0.);
}

#ifdef PRESSIO_ENABLE_TPL_EIGEN
TEST_F(epetraMultiVectorGlobSize15Fixture, mv_prod_eigen_vector)
{
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
        (*myMv_)[j][i] = (double)j;
     }
    }

    Eigen::VectorXd a(numVecs_);
    a.setConstant(1.);

    Epetra_Vector y(*contigMap_);
    y.PutScalar(0.);
    pressio::ops::product(::pressio::nontranspose{}, 1., *myMv_, a, 1., y);

    EXPECT_DOUBLE_EQ(y[0], 6.);
    EXPECT_DOUBLE_EQ(y[1], 6.);
    EXPECT_DOUBLE_EQ(y[2], 6.);
    EXPECT_DOUBLE_EQ(y[3], 6.);

    // simulate beta=0 with uninitialized y containing NaN
    const auto nan = std::nan("0");
    y.PutScalar(nan);
    pressio::ops::product(::pressio::nontranspose{}, 1., *myMv_, a, 0., y);

    EXPECT_DOUBLE_EQ(y[0], 6.);
    EXPECT_DOUBLE_EQ(y[1], 6.);
    EXPECT_DOUBLE_EQ(y[2], 6.);
    EXPECT_DOUBLE_EQ(y[3], 6.);

    // simulate alpha=0 with NaN in input matrix
    (*myMv_)[0][0] = nan;
    pressio::ops::product(::pressio::nontranspose{}, 0., *myMv_, a, 1., y);

    EXPECT_DOUBLE_EQ(y[0], 6.);
    EXPECT_DOUBLE_EQ(y[1], 6.);
    EXPECT_DOUBLE_EQ(y[2], 6.);
    EXPECT_DOUBLE_EQ(y[3], 6.);

    // alpha == beta == 0
    pressio::ops::product(::pressio::nontranspose{}, 0., *myMv_, a, 0., y);

    EXPECT_DOUBLE_EQ(y[0], 0.);
    EXPECT_DOUBLE_EQ(y[1], 0.);
    EXPECT_DOUBLE_EQ(y[2], 0.);
    EXPECT_DOUBLE_EQ(y[3], 0.);
}

TEST_F(epetraMultiVectorGlobSize15Fixture, mv_T_vector_storein_eigen_vector)
{
    auto & myMv_h = *myMv_;
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
        myMv_h[j][i] = (double)i;
     }
    }

    Epetra_Vector y(*contigMap_);
    y.PutScalar(2.);

    Eigen::VectorXd a(numVecs_);
    a.setConstant(1.);
    pressio::ops::product(::pressio::transpose{}, 1., *myMv_, y, 1., a);

    const auto a1_ref = numProc_ * 2. * 10. + 1.;
    EXPECT_DOUBLE_EQ(a(0), a1_ref);
    EXPECT_DOUBLE_EQ(a(1), a1_ref);
    EXPECT_DOUBLE_EQ(a(2), a1_ref);
    EXPECT_DOUBLE_EQ(a(3), a1_ref);

    // simulate beta=0 with uninitialized y containing NaN
    const auto nan = std::nan("0");
    a.setConstant(nan);
    pressio::ops::product(::pressio::transpose{}, 1., *myMv_, y, 0., a);

    const auto a0_ref = numProc_ * 2. * 10. + 0.;
    EXPECT_DOUBLE_EQ(a(0), a0_ref);
    EXPECT_DOUBLE_EQ(a(1), a0_ref);
    EXPECT_DOUBLE_EQ(a(2), a0_ref);
    EXPECT_DOUBLE_EQ(a(3), a0_ref);

    // simulate alpha=0 with NaN in input matrix
    myMv_h[0][0] = nan;
    pressio::ops::product(::pressio::transpose{}, 0., *myMv_, y, 1., a);

    EXPECT_DOUBLE_EQ(a(0), a0_ref);
    EXPECT_DOUBLE_EQ(a(1), a0_ref);
    EXPECT_DOUBLE_EQ(a(2), a0_ref);
    EXPECT_DOUBLE_EQ(a(3), a0_ref);

    // alpha == beta == 0
    pressio::ops::product(::pressio::transpose{}, 0., *myMv_, y, 0., a);

    EXPECT_DOUBLE_EQ(a(0), 0.);
    EXPECT_DOUBLE_EQ(a(1), 0.);
    EXPECT_DOUBLE_EQ(a(2), 0.);
    EXPECT_DOUBLE_EQ(a(3), 0.);
}
#endif
