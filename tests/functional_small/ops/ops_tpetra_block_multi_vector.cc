
#include "tpetra_block_only_fixtures.hpp"
#include "pressio/ops.hpp"

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture, multi_vector_clone)
{
  auto a = ::pressio::ops::clone(*myMv_);
  auto a_h = a.getMultiVectorView().getLocalViewDevice(Tpetra::Access::ReadOnly);
  auto myMv_h = myMv_->getMultiVectorView().getLocalViewDevice(Tpetra::Access::ReadOnly);
  ASSERT_NE(a_h.data(), myMv_h.data());

  for (int i=0; i<localSize_*blockSize_; ++i){
    for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(a_h(i,j), myMv_h(i, j));
    }
  }

  myMv_->putScalar(82347.);
  for (int i=0; i<localSize_*blockSize_; ++i){
    for (int j=0; j<numVecs_; ++j){
      EXPECT_NE(a_h(i,j), 82347.0);
    }
  }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture, multi_vector_extent)
{
  ASSERT_TRUE(pressio::ops::extent(*myMv_,0) == numProc_ * localSize_);
  ASSERT_TRUE(pressio::ops::extent(*myMv_,1) == numVecs_);
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture, multi_vector_setzero)
{
  myMv_->putScalar(23.);
  auto myMv_h = myMv_->getMultiVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(myMv_h(i,j), 23.);
    }
  }

  ::pressio::ops::set_zero(*myMv_);
  auto myMv_h2 = myMv_->getMultiVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(myMv_h2(i,j), 0.);
    }
  }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture, multi_vector_fill)
{
  ::pressio::ops::fill(*myMv_, 55.);
  auto myMv_h = myMv_->getMultiVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(myMv_h(i,j), 55.);
    }
  }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture, multi_vector_deep_copy)
{
  ::pressio::ops::fill(*myMv_, -5.);
  auto a = ::pressio::ops::clone(*myMv_);
  ::pressio::ops::deep_copy(a, *myMv_);

  auto a_h = a.getMultiVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(a_h(i,j), -5.);
    }
  }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture, multi_vector_update)
{
  auto v = ::pressio::ops::clone(*myMv_);
  ::pressio::ops::fill(v, 1.);
  auto a = ::pressio::ops::clone(*myMv_);
  ::pressio::ops::fill(a, 2.);
  ::pressio::ops::update(v, 0., a, 1.);
  auto v_h = v.getMultiVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(v_h(i,j), 2.);
    }
  }
  ::pressio::ops::update(v, a, 1.);
  for (int i=0; i<localSize_*blockSize_; ++i){
    for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(v_h(i,j), 2.);
    }
  }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture, multi_vector_update_nan)
{
    auto v = ::pressio::ops::clone(*myMv_);
    auto a = ::pressio::ops::clone(*myMv_);
    auto v_h = v.getMultiVectorView().getLocalViewHost(Tpetra::Access::ReadOnly);

    // NaN injection through alpha=0
    const auto nan = std::nan("0");
    ::pressio::ops::fill(v, nan);
    ::pressio::ops::fill(a, 1.);
    ::pressio::ops::update(v, 0., a, 2.);
    for (int i = 0; i < localSize_; ++i) {
      for (int j = 0; j < numVecs_; ++j) {
        EXPECT_DOUBLE_EQ(v_h(i, j), 2.);
      }
    }

    // NaN injection through beta=0
    ::pressio::ops::fill(a, nan);
    ::pressio::ops::update(v, -1., a, 0.);
    for (int i = 0; i < localSize_; ++i) {
      for (int j = 0; j < numVecs_; ++j) {
        EXPECT_DOUBLE_EQ(v_h(i, j), -2.);
      }
    }

    // alpha=beta=0 corner case
    ::pressio::ops::fill(v, nan);
    ::pressio::ops::update(v, 0., a, 0.);
    for (int i = 0; i < localSize_; ++i) {
      for (int j = 0; j < numVecs_; ++j) {
        EXPECT_DOUBLE_EQ(v_h(i, j), 0.);
      }
    }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture, multi_vector_min_max)
{
  ASSERT_DOUBLE_EQ(pressio::ops::min(*myMv_), 1.);
  ASSERT_DOUBLE_EQ(pressio::ops::max(*myMv_), numProc_ * localSize_ * numVecs_ * blockSize_);
}
