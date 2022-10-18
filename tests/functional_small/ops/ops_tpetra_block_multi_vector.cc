
#include "tpetra_block_only_fixtures.hpp"
#include "pressio/ops.hpp"

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture, multi_vector_clone)
{
  auto a = pressio::ops::clone(*myMv_);
  ASSERT_TRUE(a.getMultiVectorView().getLocalViewDevice(Tpetra::Access::ReadWriteStruct()).data() !=
	      myMv_->getMultiVectorView().getLocalViewDevice(Tpetra::Access::ReadWriteStruct()).data());

  auto a_h = a.getMultiVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(a_h(i,j), 0.0);
    }
  }

  myMv_->putScalar(23.);
  for (int i=0; i<localSize_*blockSize_; ++i){
    for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(a_h(i,j), 0.0);
    }
  }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture, multi_vector_extent)
{
  ASSERT_TRUE(pressio::ops::extent(*myMv_,0) == 15);
  ASSERT_TRUE(pressio::ops::extent(*myMv_,1) == 3);
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

  pressio::ops::set_zero(*myMv_);
  auto myMv_h2 = myMv_->getMultiVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(myMv_h2(i,j), 0.);
    }
  }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture, multi_vector_fill)
{
  pressio::ops::fill(*myMv_, 55.);
  auto myMv_h = myMv_->getMultiVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(myMv_h(i,j), 55.);
    }
  }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture, multi_vector_deep_copy)
{
  pressio::ops::fill(*myMv_, -5.);
  auto a = pressio::ops::clone(*myMv_);
  pressio::ops::deep_copy(a, *myMv_);

  auto a_h = a.getMultiVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(a_h(i,j), -5.);
    }
  }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture, multi_vector_update)
{
  auto v = pressio::ops::clone(*myMv_);
  pressio::ops::fill(v, 1.);
  auto a = pressio::ops::clone(*myMv_);
  pressio::ops::fill(a, 2.);
  pressio::ops::update(v, 0., a, 1.);
  auto v_h = v.getMultiVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    for (int j=0; j<numVecs_; ++j){
      EXPECT_DOUBLE_EQ(v_h(i,j), 2.);
    }
  }
}
