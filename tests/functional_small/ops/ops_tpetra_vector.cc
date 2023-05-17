
#include "tpetra_only_fixtures.hpp"
#include "pressio/ops.hpp"

// convenient alias for nice test names
using ops_tpetra = tpetraVectorGlobSize15Fixture;

TEST_F(ops_tpetra, vector_clone)
{
    auto a = pressio::ops::clone(*myVector_);
    ASSERT_TRUE(a.getLocalViewDevice(Tpetra::Access::ReadOnlyStruct()).data() !=
		myVector_->getLocalViewDevice(Tpetra::Access::ReadOnlyStruct()).data());

    auto a_h = a.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(a_h(i,0), 0.0);
    }

    myVector_->putScalar(23.);
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(a_h(i,0), 0.0);
    }
}

TEST_F(ops_tpetra, vector_extent)
{
    ASSERT_TRUE( (std::size_t)pressio::ops::extent(*myVector_,0) == (std::size_t)numProc_ * 5);
    ASSERT_TRUE( (std::size_t)pressio::ops::extent(*myVector_,1) == (std::size_t)1); // check extent over the rank
}

TEST_F(ops_tpetra, vector_deep_copy)
{
    pressio::ops::fill(*myVector_, -5.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::deep_copy(a, *myVector_);

    auto a_h = a.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(a_h(i,0), -5.);
    }
}

TEST_F(ops_tpetra, vector_setzero)
{
    myVector_->putScalar(23.);
    auto x_h = myVector_->getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(x_h(i,0), 23.);
    }

    pressio::ops::set_zero(*myVector_);
    auto x_h2 = myVector_->getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(x_h2(i,0), 0.);
    }
}

TEST_F(ops_tpetra, vector_scale)
{
    myVector_->putScalar(2.);
    pressio::ops::scale(*myVector_, 3.);
    auto x_h2 = myVector_->getLocalViewHost(Tpetra::Access::ReadOnly);
    for (int i = 0; i < localSize_; ++i){
        EXPECT_DOUBLE_EQ(x_h2(i, 0), 6.);
    }

    // NaN injection
    myVector_->putScalar(std::nan("0"));
    pressio::ops::scale(*myVector_, 0.);
    x_h2 = myVector_->getLocalViewHost(Tpetra::Access::ReadOnly);
    for (int i = 0; i < localSize_; ++i){
        EXPECT_DOUBLE_EQ(x_h2(i, 0), 0.);
    }
}

TEST_F(ops_tpetra, vector_fill)
{
    pressio::ops::fill(*myVector_, 55.);
    auto x_h = myVector_->getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(x_h(i,0), 55.);
    }
}

TEST_F(ops_tpetra, vector_abs)
{
    pressio::ops::fill(*myVector_, -5.);
    auto x_h = myVector_->getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(x_h(i,0), -5.);
    }

    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::abs(a, *myVector_);
    auto a_h = a.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(a_h(i,0), 5.);
        EXPECT_DOUBLE_EQ(x_h(i,0), -5.);
    }
}

TEST_F(ops_tpetra, vector_dot)
{
  auto a = pressio::ops::clone(*myVector_);
  auto b = pressio::ops::clone(*myVector_);
  pressio::ops::fill(a, 1.0);
  pressio::ops::fill(b, 1.0);

  auto res = ::pressio::ops::dot(a, b);
  EXPECT_DOUBLE_EQ(res, numProc_ * 5.);

  res = 0.0;
  ::pressio::ops::dot(a, b, res);
  EXPECT_DOUBLE_EQ(res, numProc_ * 5.);
}

TEST_F(ops_tpetra, vector_min_max)
{
  auto a = pressio::ops::clone(*myVector_);
  auto a2_h = a.getLocalViewHost(Tpetra::Access::ReadWrite);
  auto a_h = Kokkos::subview(a2_h, Kokkos::ALL, 0);
  for (int i = 0; i < localSize_; ++i) {
    a_h(i) = 100.0 - (rank_ * localSize_ + i);
  }
  ASSERT_DOUBLE_EQ(pressio::ops::min(a), 100. - (numProc_ * localSize_ - 1.0));
  ASSERT_DOUBLE_EQ(pressio::ops::max(a), 100.);
}

TEST_F(ops_tpetra, vector_norm2)
{
  pressio::ops::fill(*myVector_, 1.0);
  auto mynorm = pressio::ops::norm2(*myVector_);
  EXPECT_NEAR(mynorm, std::sqrt(numProc_ * 5.), 1e-15);
}

TEST_F(ops_tpetra, vector_norm1)
{
  pressio::ops::fill(*myVector_, 1.0);
  auto mynorm = pressio::ops::norm1(*myVector_);
  EXPECT_DOUBLE_EQ(mynorm, numProc_ * 5.);
}

TEST_F(ops_tpetra, vector_pow)
{
  {
    auto shift = rank_*localSize_;

    // 1. create vector x and fill with data
    vec_t x(contigMap_);
    auto xh = x.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    std::vector<ST> x0(numGlobalEntries_);
    for (int i=0; i<numGlobalEntries_; ++i) x0[i]= (double) i;
    for (int i=0; i<localSize_; ++i) xh(i,0) = x0[shift+i];

    pressio::ops::pow(x, 2.);

    // 5. check correctness
    Eigen::VectorXd g(numGlobalEntries_);
    for (int i=0; i<numGlobalEntries_; ++i){
      g(i) = i * i;
    }
    ASSERT_EQ( pressio::ops::extent(x,0), numProc_ * localSize_ );
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(xh(i,0), g(shift+i));
    }
  }
}

TEST_F(ops_tpetra, vector_absPowPos)
{
  {
    auto shift = rank_*localSize_;

    // 1. create vector x and fill with data
    vec_t y(contigMap_);
    auto yh = y.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    vec_t x(contigMap_);
    auto xh = x.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    std::vector<ST> x0(numGlobalEntries_);
    for (int i=0; i<numGlobalEntries_; ++i){
      x0[i] = (double) i;
      x0[i] = -x0[i];
    }
    for (int i=0; i<localSize_; ++i) xh(i,0) = x0[shift+i];

    // 2. compute
    ::pressio::ops::abs_pow(y, x, 3.);

    // 5. check correctness
    Eigen::VectorXd g(numGlobalEntries_);
    for (int i=0; i<numGlobalEntries_; ++i){
      g(i) = i * i * i;
    }
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(yh(i,0), g(shift+i));
    }
  }
}

TEST_F(ops_tpetra, vector_absPowNeg)
{
  {
    auto shift = rank_*localSize_;

    // 1. create vector x and fill with data
    vec_t y(contigMap_);
    auto yh = y.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    vec_t x(contigMap_);
    auto xh = x.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    std::vector<ST> x0(numGlobalEntries_);
    for (int i=0; i<numGlobalEntries_; ++i){
      x0[i] = (double) i;
      x0[i] = -x0[i];
    }
    for (int i=0; i<localSize_; ++i) xh(i,0) = x0[shift+i];

    // 2. compute
    pressio::ops::abs_pow(y, x, -3., 0.001);

    // 5. check correctness
    Eigen::VectorXd g(numGlobalEntries_);
    g(0) = 1./0.001; //because we guard against div by 0
    for (int i=1; i<numGlobalEntries_; ++i){
      g(i) = 1. / (i * i * i);
    }
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(yh(i,0), g(shift+i));
    }
  }
}


TEST_F(ops_tpetra, vector_update1_a)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    pressio::ops::update(v, 0., a, 1.);
    auto v_h = v.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 2.);
    }
}

TEST_F(ops_tpetra, vector_update1_b)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    pressio::ops::update(v, 1., a, 1.);
    auto v_h = v.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 3.);
    }
}

TEST_F(ops_tpetra, vector_update2_a)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    auto b = pressio::ops::clone(*myVector_);
    pressio::ops::fill(b, 3.);

    pressio::ops::update(v, 0., a, 1., b, 1.);
    auto v_h = v.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 5.);
    }
}

TEST_F(ops_tpetra, vector_update2_b)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    auto b = pressio::ops::clone(*myVector_);
    pressio::ops::fill(b, 3.);

    pressio::ops::update(v, 1., a, 1., b, 1.);
    auto v_h = v.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 6.);
    }
}

TEST_F(ops_tpetra, vector_update3_a)
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
    auto v_h = v.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 9.);
    }
}

TEST_F(ops_tpetra, vector_update3_b)
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
    auto v_h = v.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 10.);
    }
}

TEST_F(ops_tpetra, vector_update4_a)
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
    auto v_h = v.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 14.);
    }
}

TEST_F(ops_tpetra, vector_update4_b)
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
    auto v_h = v.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 15.);
    }
}

TEST_F(ops_tpetra, vector_update_nan1)
{
    auto v = pressio::ops::clone(*myVector_);
    auto v2_h = v.getLocalViewHost(Tpetra::Access::ReadOnly);
    auto v_h = Kokkos::subview(v2_h, Kokkos::ALL, 0);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 1.);
    auto vecOfNans = pressio::ops::clone(*myVector_);
    pressio::ops::fill(vecOfNans, std::nan("0"));

    // Note: this test covers just enough nan/non-nan combinations
    // to trigger and verify all execution paths in our update()
    // implementations, which include anti-NaN-injection variants
    pressio::ops::update(v, 1., vecOfNans, 0.);
    EXPECT_DOUBLE_EQ(v_h(0), 1.0);

    pressio::ops::update(v, 1., vecOfNans, 0., vecOfNans, 0.);
    EXPECT_DOUBLE_EQ(v_h(0), 1.0);
    pressio::ops::update(v, 1., a, 1., vecOfNans, 0.);
    EXPECT_DOUBLE_EQ(v_h(0), 2.);

    pressio::ops::update(v, 1., vecOfNans, 0., vecOfNans, 0., vecOfNans, 0.);
    EXPECT_DOUBLE_EQ(v_h(0), 2.0);
    pressio::ops::update(v, 1., a, 1., vecOfNans, 0., a, 1.);
    EXPECT_DOUBLE_EQ(v_h(0), 4.);
    pressio::ops::update(v, 1., a, 1., a, 1., vecOfNans, 0.);
    EXPECT_DOUBLE_EQ(v_h(0), 6.);

    pressio::ops::update(v, 1., vecOfNans, 0., vecOfNans, 0., vecOfNans, 0., vecOfNans, 0.);
    EXPECT_DOUBLE_EQ(v_h(0), 6.0);
    pressio::ops::update(v, 1., a, 1., vecOfNans, 0., a, 1., a, 1.);
    EXPECT_DOUBLE_EQ(v_h(0), 9.);
    pressio::ops::update(v, 1., a, 1., a, 1., vecOfNans, 0., a, 1.);
    EXPECT_DOUBLE_EQ(v_h(0), 12.);
    pressio::ops::update(v, 1., a, 1., a, 1., a, 1., vecOfNans, 0.);
    EXPECT_DOUBLE_EQ(v_h(0), 15.);
}

// injects NaN through the updated vector
TEST_F(ops_tpetra, vector_update_nan2)
{
  const auto nan = std::nan("0");
  auto v = pressio::ops::clone(*myVector_);
  auto v2_h = v.getLocalViewHost(Tpetra::Access::ReadOnly);
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

TEST_F(ops_tpetra, vector_elementwiseMultiply)
{
  auto shift = rank_*localSize_;

  // c = beta*c + alpha*x*b;

  // 1. create vector x and fill with data
  vec_t x(contigMap_);
  auto xh = x.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  std::vector<ST> x0(numGlobalEntries_);
  for (int i=0; i<numGlobalEntries_; ++i) x0[i]= (double) i;
  for (int i=0; i<localSize_; ++i) xh(i,0) = x0[shift+i];

  // 2. create vector b and fill with data
  vec_t b(contigMap_);
  auto bh = b.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  std::vector<ST> b0(numGlobalEntries_);
  for (int i=0; i<numGlobalEntries_; ++i) b0[i]= (double) i + 2.;
  for (int i=0; i<localSize_; ++i) bh(i,0) = b0[shift+i];

  // 3. compute ewise multiply
  const auto alpha1 = 1.;
  vec_t c(contigMap_);
  const double offset = -5.0;
  pressio::ops::fill(c, offset);
  pressio::ops::elementwise_multiply(alpha1, x, b, 1., c);
  auto ch = c.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  vec_t c0(contigMap_);
  pressio::ops::fill(c0, std::nan("0")); // test beta=0 with simulated NaN in uninitialized c
  pressio::ops::elementwise_multiply(alpha1, x, b, 0., c0);
  auto ch0 = c0.getLocalViewHost(Tpetra::Access::ReadWriteStruct());

  // 5. check correctness
  for (int i=0; i<localSize_; ++i){
    const auto gold = alpha1 * xh(i, 0) * bh(i, 0);
    EXPECT_DOUBLE_EQ(ch0(i,0), gold);
    EXPECT_DOUBLE_EQ(ch(i,0), gold + offset);
  }
}
