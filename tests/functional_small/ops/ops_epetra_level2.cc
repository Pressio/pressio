
#include "epetra_only_fixtures.hpp"
#include "pressio/ops.hpp"

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

    EXPECT_DOUBLE_EQ(a(0), 61.);
    EXPECT_DOUBLE_EQ(a(1), 61.);
    EXPECT_DOUBLE_EQ(a(2), 61.);
    EXPECT_DOUBLE_EQ(a(3), 61.);
}
#endif
