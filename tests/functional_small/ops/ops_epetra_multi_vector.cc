
#include "epetra_only_fixtures.hpp"
#include "pressio/ops.hpp"

TEST_F(epetraMultiVectorGlobSize15Fixture, multi_vector_clone)
{
    auto a = ::pressio::ops::clone(*myMv_);
    ASSERT_TRUE(a.Values() != myMv_->Values());

    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
        EXPECT_DOUBLE_EQ(a[j][i], (*myMv_)[j][i]);
     }
    }

    myMv_->PutScalar(82631.);
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
        EXPECT_NE(a[j][i], 82631.);
     }
    }
}

TEST_F(epetraMultiVectorGlobSize15Fixture, multi_vector_extent)
{
    ASSERT_TRUE(pressio::ops::extent(*myMv_,0) == numProc_ * localSize_);
    ASSERT_TRUE(pressio::ops::extent(*myMv_,1) == numVecs_);
    ASSERT_TRUE(pressio::ops::extent(*myMv_,2) == 1); // check extent over the rank
}

TEST_F(epetraMultiVectorGlobSize15Fixture, multi_vector_deep_copy)
{
    myMv_->PutScalar(-5.);
    auto a = ::pressio::ops::clone(*myMv_);
    ::pressio::ops::deep_copy(a, *myMv_);
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

    ::pressio::ops::set_zero(*myMv_);
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
        EXPECT_DOUBLE_EQ((*myMv_)[j][i], 0.);
     }
    }
}

TEST_F(epetraMultiVectorGlobSize15Fixture, multi_vector_scale)
{
    myMv_->PutScalar(2.);
    ::pressio::ops::scale(*myMv_, 3);
    for (int i = 0; i < localSize_; ++i){
     for (int j = 0; j < numVecs_; ++j){
        EXPECT_DOUBLE_EQ((*myMv_)[j][i], 6.);
     }
    }
}

TEST_F(epetraMultiVectorGlobSize15Fixture, multi_vector_fill)
{
    ::pressio::ops::fill(*myMv_, 55.);
    auto & myMv_h = *myMv_;
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
        EXPECT_DOUBLE_EQ(myMv_h[j][i], 55.);
     }
    }
}

TEST_F(epetraMultiVectorGlobSize15Fixture, multi_vector_update1_a)
{
    auto v = ::pressio::ops::clone(*myMv_);
    ::pressio::ops::fill(v, 1.);
    auto a = ::pressio::ops::clone(*myMv_);
    ::pressio::ops::fill(a, 2.);
    ::pressio::ops::update(v, 0., a, 1.);
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(v[j][i], 2.);
     }
    }
    ::pressio::ops::update(v, a, 1.);
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(v[j][i], 2.);
     }
    }
}

TEST_F(epetraMultiVectorGlobSize15Fixture, multi_vector_update1_b)
{
    auto v = ::pressio::ops::clone(*myMv_);
    ::pressio::ops::fill(v, 1.);
    auto a = ::pressio::ops::clone(*myMv_);
    ::pressio::ops::fill(a, 2.);
    ::pressio::ops::update(v, 1., a, 1.);
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(v[j][i], 3.);
     }
    }
}

TEST_F(epetraMultiVectorGlobSize15Fixture, multi_vector_update2_nan)
{
    auto v = ::pressio::ops::clone(*myMv_);
    auto a = ::pressio::ops::clone(*myMv_);

    // NaN injection through alpha=0
    const auto nan = std::nan("0");
    ::pressio::ops::fill(v, nan);
    ::pressio::ops::fill(a, 1.);
    ::pressio::ops::update(v, 0., a, 2.);
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(v[j][i], 2.);
     }
    }

    // NaN injection through beta=0
    ::pressio::ops::fill(a, nan);
    ::pressio::ops::update(v, -1., a, 0.);
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(v[j][i], -2.);
     }
    }

    // alpha=beta=0 corner case
    ::pressio::ops::fill(v, nan);
    ::pressio::ops::update(v, 0., a, 0.);
    for (int i=0; i<localSize_; ++i){
     for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(v[j][i], 0.);
     }
    }
}

TEST_F(epetraMultiVectorGlobSize15Fixture, multi_vector_min_max)
{
  ASSERT_DOUBLE_EQ(pressio::ops::min(*myMv_), 1.);
  ASSERT_DOUBLE_EQ(pressio::ops::max(*myMv_), numProc_ * localSize_ * numVecs_);
}
