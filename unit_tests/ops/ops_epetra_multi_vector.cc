
#include "epetra_only_fixtures.hpp"
#include "pressio_ops.hpp"

TEST_F(epetraMultiVectorGlobSize15Fixture, multi_vector_clone)
{
    auto a = pressio::ops::clone(*myMv_);
    ASSERT_TRUE(a.Values() != myMv_->Values());

    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
        EXPECT_DOUBLE_EQ(a[j][i], 0.0);
     }
    }

    myMv_->PutScalar(23.);
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
        EXPECT_DOUBLE_EQ(a[j][i], 0.0);
     }
    }
}

TEST_F(epetraMultiVectorGlobSize15Fixture, multi_vector_extent)
{
    ASSERT_TRUE(pressio::ops::extent(*myMv_,0) == 15);
    ASSERT_TRUE(pressio::ops::extent(*myMv_,1) == 4);
}

TEST_F(epetraMultiVectorGlobSize15Fixture, multi_vector_deep_copy)
{
    myMv_->PutScalar(-5.);
    auto a = pressio::ops::clone(*myMv_);
    pressio::ops::deep_copy(a, *myMv_);
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){        
        EXPECT_DOUBLE_EQ(a[j][i], -5.);
     }
    }
}

TEST_F(epetraMultiVectorGlobSize15Fixture, multi_vector_setzero)
{
    myMv_->PutScalar(23.);
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
        EXPECT_DOUBLE_EQ((*myMv_)[j][i], 23.);
     }
    }

    pressio::ops::set_zero(*myMv_);
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
        EXPECT_DOUBLE_EQ((*myMv_)[j][i], 0.);
     } 
    }
}

TEST_F(epetraMultiVectorGlobSize15Fixture, multi_vector_fill)
{
    pressio::ops::fill(*myMv_, 55.);
    auto & myMv_h = *myMv_;
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
        EXPECT_DOUBLE_EQ(myMv_h[j][i], 55.);
     }
    }
}

TEST_F(epetraMultiVectorGlobSize15Fixture, multi_vector_update1_a)
{
    auto v = pressio::ops::clone(*myMv_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myMv_);
    pressio::ops::fill(a, 2.);
    pressio::ops::update(v, a, 1.);
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(v[j][i], 2.);
     }
    }
}

TEST_F(epetraMultiVectorGlobSize15Fixture, multi_vector_update1_b)
{
    auto v = pressio::ops::clone(*myMv_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myMv_);
    pressio::ops::fill(a, 2.);
    pressio::ops::update(v, 1.0, a, 1.);
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(v[j][i], 3.);
     }
    }
}
