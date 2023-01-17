
#include "tpetra_block_only_fixtures.hpp"
#include "pressio/ops.hpp"

// convenient alias for nice test names
using ops_tpetra_block = tpetraBlockVectorGlobSize15BlockSize5Fixture;

TEST_F(ops_tpetra_block, vector_clone)
{
  auto a = pressio::ops::clone(*myVector_);
  ASSERT_TRUE(a.getVectorView().getLocalViewDevice(Tpetra::Access::ReadWriteStruct()).data()
    != myVector_->getVectorView().getLocalViewDevice(Tpetra::Access::ReadWriteStruct()).data());

  auto a_h = a.getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
      EXPECT_DOUBLE_EQ(a_h(i,0), 0.0);
  }

  myVector_->putScalar(23.);
  for (int i=0; i<localSize_*blockSize_; ++i){
      EXPECT_DOUBLE_EQ(a_h(i,0), 0.0);
  }
}

TEST_F(ops_tpetra_block, vector_extent)
{
    ASSERT_TRUE(pressio::ops::extent(*myVector_,0) == numProc_ * 5.);
}

TEST_F(ops_tpetra_block, vector_deep_copy)
{
  myVector_->putScalar(-5.);
  auto a = pressio::ops::clone(*myVector_);
  pressio::ops::deep_copy(a, *myVector_);

  auto a_h = a.getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    EXPECT_DOUBLE_EQ(a_h(i,0), -5.);
  }
}

TEST_F(ops_tpetra_block, vector_setzero)
{
  myVector_->putScalar(23.);
  auto x_h = myVector_->getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    EXPECT_DOUBLE_EQ(x_h(i,0), 23.);
  }

  pressio::ops::set_zero(*myVector_);
  auto x_h2 = myVector_->getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    EXPECT_DOUBLE_EQ(x_h2(i,0), 0.);
  }
}

TEST_F(ops_tpetra_block, vector_fill)
{
  pressio::ops::fill(*myVector_, 55.);
  auto x_h = myVector_->getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    EXPECT_DOUBLE_EQ(x_h(i,0), 55.);
  }
}

TEST_F(ops_tpetra_block, vector_abs)
{
  pressio::ops::fill(*myVector_, -5.);
  auto x_h = myVector_->getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    EXPECT_DOUBLE_EQ(x_h(i,0), -5.);
  }

  auto a = pressio::ops::clone(*myVector_);
  pressio::ops::abs(a, *myVector_);
  auto a_h = a.getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    EXPECT_DOUBLE_EQ(a_h(i,0), 5.);
    EXPECT_DOUBLE_EQ(x_h(i,0), -5.);
  }
}

TEST_F(ops_tpetra_block, vector_dot)
{
  auto a = pressio::ops::clone(*myVector_);
  auto b = pressio::ops::clone(*myVector_);
  pressio::ops::fill(a, 1.0);
  pressio::ops::fill(b, 1.0);

  auto res = ::pressio::ops::dot(a, b);
  EXPECT_DOUBLE_EQ(res, numProc_ * 5. * blockSize_);

  res = 0.0;
  ::pressio::ops::dot(a, b, res);
  EXPECT_DOUBLE_EQ(res, numProc_ * 5. * blockSize_);
}

TEST_F(ops_tpetra_block, vector_norm2)
{
  pressio::ops::fill(*myVector_, 1.0);
  auto mynorm = pressio::ops::norm2(*myVector_);
  EXPECT_NEAR(mynorm, std::sqrt(numProc_ * 25.0), 1e-15);
}

TEST_F(ops_tpetra_block, vector_norm1)
{
  pressio::ops::fill(*myVector_, 1.0);
  auto mynorm = pressio::ops::norm1(*myVector_);
  EXPECT_DOUBLE_EQ(mynorm, numProc_ * 25.);
}

TEST_F(ops_tpetra_block, vector_pow)
{
  pressio::ops::fill(*myVector_, 3.0);
  pressio::ops::pow(*myVector_, 2.);

  auto a_h = myVector_->getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    EXPECT_DOUBLE_EQ(a_h(i,0), 9.);
  }
}

TEST_F(ops_tpetra_block, vector_absPowPos)
{
  pressio::ops::fill(*myVector_, -3.0);
  auto y = ::pressio::ops::clone(*myVector_);
  pressio::ops::abs_pow(y, *myVector_, 3.);

  auto a_h = y.getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    EXPECT_DOUBLE_EQ(a_h(i,0), 27.);
  }
}


TEST_F(ops_tpetra_block, vector_absPowNeg)
{
  pressio::ops::fill(*myVector_, -3.0);
  auto y = ::pressio::ops::clone(*myVector_);
  pressio::ops::abs_pow(y, *myVector_, -3., 1.e-6);

  auto a_h = y.getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    EXPECT_DOUBLE_EQ(a_h(i,0), 1./27.);
  }
}

TEST_F(ops_tpetra_block, vector_update1_a)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    pressio::ops::update(v, 0., a, 1.);
    auto v_h = v.getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_*blockSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 2.);
    }
}

TEST_F(ops_tpetra_block, vector_update1_b)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    pressio::ops::update(v, 1., a, 1.);
    auto v_h = v.getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_*blockSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 3.);
    }
}

TEST_F(ops_tpetra_block, vector_update2_a)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    auto b = pressio::ops::clone(*myVector_);
    pressio::ops::fill(b, 3.);

    pressio::ops::update(v, 0., a, 1., b, 1.);
    auto v_h = v.getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_*blockSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 5.);
    }
}

TEST_F(ops_tpetra_block, vector_update2_b)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    auto b = pressio::ops::clone(*myVector_);
    pressio::ops::fill(b, 3.);

    pressio::ops::update(v, 1., a, 1., b, 1.);
    auto v_h = v.getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_*blockSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 6.);
    }
}

TEST_F(ops_tpetra_block, vector_update3_a)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    auto b = pressio::ops::clone(*myVector_);
    pressio::ops::fill(b, 3.);
    auto c = pressio::ops::clone(*myVector_);
    pressio::ops::fill(c, 4.);

    pressio::ops::update(v, 0., a, 1., b, 1., c, 1.);
    auto v_h = v.getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_*blockSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 9.);
    }
}

TEST_F(ops_tpetra_block, vector_update3_b)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    auto b = pressio::ops::clone(*myVector_);
    pressio::ops::fill(b, 3.);
    auto c = pressio::ops::clone(*myVector_);
    pressio::ops::fill(c, 4.);

    pressio::ops::update(v, 1., a, 1., b, 1., c, 1.);
    auto v_h = v.getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_*blockSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 10.);
    }
}

TEST_F(ops_tpetra_block, vector_update4_a)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    auto b = pressio::ops::clone(*myVector_);
    pressio::ops::fill(b, 3.);
    auto c = pressio::ops::clone(*myVector_);
    pressio::ops::fill(c, 4.);
    auto d = pressio::ops::clone(*myVector_);
    pressio::ops::fill(d, 5.);

    pressio::ops::update(v, 0., a, 1., b, 1., c, 1., d, 1.);
    auto v_h = v.getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_*blockSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 14.);
    }
}

TEST_F(ops_tpetra_block, vector_update4_b)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    auto b = pressio::ops::clone(*myVector_);
    pressio::ops::fill(b, 3.);
    auto c = pressio::ops::clone(*myVector_);
    pressio::ops::fill(c, 4.);
    auto d = pressio::ops::clone(*myVector_);
    pressio::ops::fill(d, 5.);

    pressio::ops::update(v, 1., a, 1., b, 1., c, 1., d, 1.);
    auto v_h = v.getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_*blockSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 15.);
    }
}

TEST_F(ops_tpetra_block, vector_update_nan1)
{
    auto v = pressio::ops::clone(*myVector_);
    auto v2_h = v.getVectorView().getLocalViewHost(Tpetra::Access::ReadOnly);
    auto v_h = Kokkos::subview(v2_h, Kokkos::ALL, 0);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 1.);
    auto nan = pressio::ops::clone(*myVector_);
    pressio::ops::fill(nan, std::nan("0"));

    // Note: this test covers just enough nan/non-nan combinations
    // to trigger and verify all execution paths in our update()
    // implementations, which include anti-NaN-injection variants
    pressio::ops::update(v, 1., nan, 0.);
    EXPECT_DOUBLE_EQ(v_h(0), 1.0);

    pressio::ops::update(v, 1., nan, 0., nan, 0.);
    EXPECT_DOUBLE_EQ(v_h(0), 1.0);
    pressio::ops::update(v, 1., a, 1., nan, 0.);
    EXPECT_DOUBLE_EQ(v_h(0), 2.);

    pressio::ops::update(v, 1., nan, 0., nan, 0., nan, 0.);
    EXPECT_DOUBLE_EQ(v_h(0), 2.0);
    pressio::ops::update(v, 1., a, 1., nan, 0., a, 1.);
    EXPECT_DOUBLE_EQ(v_h(0), 4.);
    pressio::ops::update(v, 1., a, 1., a, 1., nan, 0.);
    EXPECT_DOUBLE_EQ(v_h(0), 6.);

    pressio::ops::update(v, 1., nan, 0., nan, 0., nan, 0., nan, 0.);
    EXPECT_DOUBLE_EQ(v_h(0), 6.0);
    pressio::ops::update(v, 1., a, 1., nan, 0., a, 1., a, 1.);
    EXPECT_DOUBLE_EQ(v_h(0), 9.);
    pressio::ops::update(v, 1., a, 1., a, 1., nan, 0., a, 1.);
    EXPECT_DOUBLE_EQ(v_h(0), 12.);
    pressio::ops::update(v, 1., a, 1., a, 1., a, 1., nan, 0.);
    EXPECT_DOUBLE_EQ(v_h(0), 15.);
}

// injects NaN through the updated vector
TEST_F(ops_tpetra_block, vector_update_nan2)
{
  const auto nan = std::nan("0");
  auto v = pressio::ops::clone(*myVector_);
  auto v2_h = v.getVectorView().getLocalViewHost(Tpetra::Access::ReadOnly);
  auto v_h = Kokkos::subview(v2_h, Kokkos::ALL, 0);
  auto a = pressio::ops::clone(*myVector_);
  pressio::ops::fill(a, 1.);

  pressio::ops::fill(v, nan);
  pressio::ops::update(v, 0., a, 1.);
  EXPECT_DOUBLE_EQ(v_h(0), 1.0);

  pressio::ops::fill(v, nan);
  pressio::ops::update(v, 0., a, 0.);
  EXPECT_DOUBLE_EQ(v_h(0), 0.0);

  pressio::ops::fill(v, nan);
  pressio::ops::update(v, 0., a, 1., a, 1.);
  EXPECT_DOUBLE_EQ(v_h(0), 2.0);

  pressio::ops::fill(v, nan);
  pressio::ops::update(v, 0., a, 1., a, 1., a, 1.);
  EXPECT_DOUBLE_EQ(v_h(0), 3.0);

  pressio::ops::fill(v, nan);
  pressio::ops::update(v, 0., a, 1., a, 1., a, 1., a, 1.);
  EXPECT_DOUBLE_EQ(v_h(0), 4.0);
}

TEST_F(ops_tpetra_block, vector_elementwiseMultiply)
{
  // computing elementwise:  y = beta * y + alpha * x * z

  auto x = pressio::ops::clone(*myVector_);
  pressio::ops::fill(x, 3.);

  auto z = pressio::ops::clone(*myVector_);
  pressio::ops::fill(z, 2.);

  auto y = pressio::ops::clone(*myVector_);
  pressio::ops::fill(y, 2.);

  pressio::ops::elementwise_multiply(2., x,z, 3., z);
  auto z_h = z.getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    EXPECT_DOUBLE_EQ(z_h(i,0), 18.);
  }

  // test beta=0 with simulated NaN in uninitialized y
  pressio::ops::fill(x, 3.);
  pressio::ops::fill(z, 2.);
  pressio::ops::fill(y, std::nan("0"));
  pressio::ops::elementwise_multiply(2., x,z, 0., y);
  auto y_h = y.getVectorView().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    EXPECT_DOUBLE_EQ(y_h(i,0), 12);
  }

}
