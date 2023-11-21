
#include "tpetra_only_fixtures.hpp"
#include "pressio/ops.hpp"

using ops_tpetra = tpetraMultiVectorGlobSize15Fixture;

TEST_F(ops_tpetra, column_expr_extent)
{
  auto e = pressio::column(*myMv_, 0);
  ASSERT_TRUE( pressio::ops::extent(e, 0) == std::size_t(numGlobalEntries_) );
}

TEST_F(ops_tpetra, column_expr_setzero)
{
  auto e = pressio::column(*myMv_, 0);

  // before setting to zero
  auto v = e.native();
  auto v_h = v.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
  std::vector<ST> gold(localSize_);
  if (rank_ == 0){ gold = {1,5,9,13,17}; }
  else if (rank_==1){ gold = {21,25,29,33,37}; }
  else{ gold = {41,45,49,53,57}; }
  for (std::size_t i=0; i<e.extentLocal(0); ++i){
    ASSERT_TRUE(v_h(i,0) == gold[i]);
  }

  pressio::ops::set_zero(e);

  std::vector<ST> goldAfter(localSize_, 0);
  for (std::size_t i=0; i<e.extentLocal(0); ++i){
    ASSERT_TRUE(v_h(i,0) == goldAfter[i]);
  }
}

TEST_F(ops_tpetra, column_expr_scale)
{
  auto e = pressio::column(*myMv_, 0);

  auto v = e.native();
  auto v_h = v.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());

  // before
  {
    std::vector<ST> gold(localSize_);
    if (rank_ == 0){ gold = {1,5,9,13,17}; }
    else if (rank_==1){ gold = {21,25,29,33,37}; }
    else{ gold = {41,45,49,53,57}; }
    for (std::size_t i=0; i<e.extentLocal(0); ++i){
      ASSERT_TRUE(v_h(i,0) == gold[i]);
    }
  }

  pressio::ops::scale(e, 2.);
  {
    std::vector<ST> gold(localSize_);
    if (rank_ == 0){ gold = {1,5,9,13,17}; }
    else if (rank_==1){ gold = {21,25,29,33,37}; }
    else{ gold = {41,45,49,53,57}; }
    for (std::size_t i=0; i<e.extentLocal(0); ++i){
      ASSERT_TRUE(v_h(i,0) == gold[i]*2.);
    }
  }
}

TEST_F(ops_tpetra, column_expr_fill)
{
  auto e = pressio::column(*myMv_, 0);

  auto v = e.native();
  auto v_h = v.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());

  // before
  {
    std::vector<ST> gold(localSize_);
    if (rank_ == 0){ gold = {1,5,9,13,17}; }
    else if (rank_==1){ gold = {21,25,29,33,37}; }
    else{ gold = {41,45,49,53,57}; }
    for (std::size_t i=0; i<e.extentLocal(0); ++i){
      ASSERT_TRUE(v_h(i,0) == gold[i]);
    }
  }

  pressio::ops::fill(e, 33.);
  {
    std::vector<ST> gold(localSize_);
    for (std::size_t i=0; i<e.extentLocal(0); ++i){
      ASSERT_TRUE(v_h(i,0) == 33.);
    }
  }
}

TEST_F(ops_tpetra, column_expr_abs)
{
  auto e = pressio::column(*myMv_, 0);

  auto v = e.native();
  auto v_h = v.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());

  // before
  pressio::ops::scale(e, -1);
  {
    std::vector<ST> gold(localSize_);
    if      (rank_ == 0){ gold = {1,5,9,13,17};    }
    else if (rank_==1)  { gold = {21,25,29,33,37}; }
    else                { gold = {41,45,49,53,57}; }
    for (std::size_t i=0; i<e.extentLocal(0); ++i){
      ASSERT_TRUE(v_h(i,0) == -gold[i]);
    }
  }

  // after
  vec_t vec(contigMap_);
  pressio::ops::abs(vec, e);
  {
    auto v2_h = vec.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
    std::vector<ST> gold(localSize_);
    if      (rank_ == 0){ gold = {1,5,9,13,17};    }
    else if (rank_==1)  { gold = {21,25,29,33,37}; }
    else                { gold = {41,45,49,53,57}; }
    for (std::size_t i=0; i<e.extentLocal(0); ++i){
      ASSERT_TRUE(v2_h(i,0) == gold[i]);
    }
  }
}

TEST_F(ops_tpetra, column_expr_dot_with_vector)
{
  auto a = pressio::column(*myMv_, 0);
  pressio::ops::fill(a, 1.0);
  vec_t b(contigMap_);
  pressio::ops::fill(b, 1.);

  {
    auto res = ::pressio::ops::dot(a, b);
    EXPECT_DOUBLE_EQ(res, numProc_ * 5.);
    res = 0.0;
    ::pressio::ops::dot(a, b, res);
    EXPECT_DOUBLE_EQ(res, numProc_ * 5.);
  }

  {
    auto res = ::pressio::ops::dot(b, a);
    EXPECT_DOUBLE_EQ(res, numProc_ * 5.);
    res = 0.0;
    ::pressio::ops::dot(a, b, res);
    EXPECT_DOUBLE_EQ(res, numProc_ * 5.);
  }
}

TEST_F(ops_tpetra, column_expr_dot_with_expr)
{
  auto a = pressio::column(*myMv_, 0);
  auto b = pressio::column(*myMv_, 1);
  pressio::ops::fill(a, 1.0);
  pressio::ops::fill(b, 1.0);
  auto res = ::pressio::ops::dot(a, b);
  EXPECT_DOUBLE_EQ(res, numProc_ * 5.);
  res = 0.0;
  ::pressio::ops::dot(a, b, res);
  EXPECT_DOUBLE_EQ(res, numProc_ * 5.);
}

TEST_F(ops_tpetra, column_expr_min_max)
{
  auto e = pressio::column(*myMv_, 0);
  ASSERT_DOUBLE_EQ(pressio::ops::min(e), 1.);
  ASSERT_DOUBLE_EQ(pressio::ops::max(e), 57.);
}

TEST_F(ops_tpetra, column_expr_norm2)
{
  auto e = pressio::column(*myMv_, 0);
  pressio::ops::fill(e, 1.0);
  auto mynorm = pressio::ops::norm2(e);
  EXPECT_NEAR(mynorm, std::sqrt(numProc_ * 5.), 1e-15);
}

TEST_F(ops_tpetra, column_expr_norm1)
{
  auto e = pressio::column(*myMv_, 0);
  pressio::ops::fill(e, 1.0);
  auto mynorm = pressio::ops::norm1(e);
  EXPECT_DOUBLE_EQ(mynorm, numProc_ * 5.);
}

TEST_F(ops_tpetra, column_expr_pow)
{
  {
    auto shift = rank_*localSize_;

    // 1. create mvector x and fill with data
    mvec_t x(contigMap_, 3);
    auto xh = x.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    std::vector<ST> x0(numGlobalEntries_);
    for (int i=0; i<numGlobalEntries_; ++i){ x0[i] = (double) i; }

    for (int i=0; i<localSize_; ++i){
      for (int j=0; j<3; ++j){
	xh(i,j) = x0[shift+i];
      }
    }

    auto e = pressio::column(x, 0);
    pressio::ops::pow(e, 2.);

    // 5. check correctness
    Eigen::VectorXd g(numGlobalEntries_);
    for (int i=0; i<numGlobalEntries_; ++i){ g(i) = i * i; }
    ASSERT_EQ( pressio::ops::extent(x,0), numProc_ * localSize_ );
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(xh(i,0), g(shift+i));
    }
  }
}

TEST_F(ops_tpetra, column_expr_absPowPos)
{
  {
    auto shift = rank_*localSize_;

    // 1. create and fill with data
    mvec_t y(contigMap_, 3);
    auto yh = y.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    mvec_t x(contigMap_, 3);
    auto xh = x.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    std::vector<ST> x0(numGlobalEntries_);
    for (int i=0; i<numGlobalEntries_; ++i){
      x0[i] = (double) i;
      x0[i] = -x0[i];
    }
    for (int i=0; i<localSize_; ++i){
      for (int j=0; j<3; ++j){
	xh(i,j) = x0[shift+i];
      }
    }

    // 2. compute
    auto e_x = pressio::column(x, 0);
    auto e_y = pressio::column(y, 0);
    ::pressio::ops::abs_pow(e_y, e_x, 3.);

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

TEST_F(ops_tpetra, column_expr_absPowNeg)
{
  {
    auto shift = rank_*localSize_;

    // 1. create and fill with data
    mvec_t y(contigMap_, 3);
    auto yh = y.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    mvec_t x(contigMap_, 3);
    auto xh = x.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    std::vector<ST> x0(numGlobalEntries_);
    for (int i=0; i<numGlobalEntries_; ++i){
      x0[i] = (double) i;
      x0[i] = -x0[i];
    }
    for (int i=0; i<localSize_; ++i) {
      for (int j=0; j<3; ++j){
	xh(i,j) = x0[shift+i];
      }
    }

    // 2. compute
    auto e_x = pressio::column(x, 0);
    auto e_y = pressio::column(y, 0);
    pressio::ops::abs_pow(e_y, e_x, -3., 0.001);

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


TEST_F(ops_tpetra, column_expr_update1_a)
{
    auto a = pressio::column(*myMv_, 0);
    pressio::ops::fill(a, 1.);
    auto b = pressio::column(*myMv_, 1);
    pressio::ops::fill(b, 2.);
    pressio::ops::update(a, 0., b, 1.);
    auto a_h = a.native().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(a_h(i,0), 2.);
    }
}

TEST_F(ops_tpetra, column_expr_update1_b)
{
    auto v = pressio::column(*myMv_, 0);
    pressio::ops::fill(v, 1.);
    auto a = pressio::column(*myMv_, 1);
    pressio::ops::fill(a, 2.);
    pressio::ops::update(v, 1., a, 1.);
    auto v_h = v.native().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 3.);
    }
}

TEST_F(ops_tpetra, column_expr_update2_a)
{
    auto v = pressio::column(*myMv_, 0);
    pressio::ops::fill(v, 1.);
    auto a = pressio::column(*myMv_, 1);
    pressio::ops::fill(a, 2.);
    auto b = pressio::column(*myMv_, 2);
    pressio::ops::fill(b, 3.);

    pressio::ops::update(v, 0., a, 1., b, 1.);
    auto v_h = v.native().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 5.);
    }
}

TEST_F(ops_tpetra, column_expr_update2_b)
{
    auto v = pressio::column(*myMv_, 0);
    pressio::ops::fill(v, 1.);
    auto a = pressio::column(*myMv_, 1);
    pressio::ops::fill(a, 2.);
    auto b = pressio::column(*myMv_, 2);
    pressio::ops::fill(b, 3.);

    pressio::ops::update(v, 1., a, 1., b, 1.);
    auto v_h = v.native().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 6.);
    }
}

TEST_F(ops_tpetra, column_expr_update3_a)
{
    auto v = pressio::column(*myMv_, 0);
    pressio::ops::fill(v, 1.);
    auto a = pressio::column(*myMv_, 1);
    pressio::ops::fill(a, 2.);
    auto b = pressio::column(*myMv_, 2);
    pressio::ops::fill(b, 3.);
    auto c = pressio::column(*myMv_, 3);
    pressio::ops::fill(c, 4.);

    pressio::ops::update(v, 0., a, 1., b, 1., c, 1.);
    auto v_h = v.native().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 9.);
    }
}

TEST_F(ops_tpetra, column_expr_update3_b)
{
    auto v = pressio::column(*myMv_, 0);
    pressio::ops::fill(v, 1.);
    auto a = pressio::column(*myMv_, 1);
    pressio::ops::fill(a, 2.);
    auto b = pressio::column(*myMv_, 2);
    pressio::ops::fill(b, 3.);
    auto c = pressio::column(*myMv_, 3);
    pressio::ops::fill(c, 4.);

    pressio::ops::update(v, 1., a, 1., b, 1., c, 1.);
    auto v_h = v.native().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 10.);
    }
}

TEST_F(ops_tpetra, column_expr_update_nan1)
{
    auto v = pressio::column(*myMv_, 0);
    auto v2_h = v.native().getLocalViewHost(Tpetra::Access::ReadOnly);
    auto v_h = Kokkos::subview(v2_h, Kokkos::ALL, 0);
    pressio::ops::fill(v, 1.);
    auto a = pressio::column(*myMv_, 1);
    pressio::ops::fill(a, 1.);
    auto vecOfNans = pressio::column(*myMv_, 2);
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
TEST_F(ops_tpetra, column_expr_update_nan2)
{
  const auto nan = std::nan("0");
  auto v = pressio::column(*myMv_, 0);
  auto v2_h = v.native().getLocalViewHost(Tpetra::Access::ReadOnly);
  auto v_h = Kokkos::subview(v2_h, Kokkos::ALL, 0);
  auto a = pressio::column(*myMv_, 1);
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

TEST_F(ops_tpetra, column_expr_elementwiseMultiply)
{
  auto shift = rank_*localSize_;

  // c = beta*c + alpha*x*b;

  // 1. create and fill with data
  mvec_t x(contigMap_, 3);
  auto xh = x.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  std::vector<ST> x0(numGlobalEntries_);
  for (int i=0; i<numGlobalEntries_; ++i) x0[i]= (double) i;
  for (int i=0; i<localSize_; ++i) {
    for (int j=0; j<3; ++j){
      xh(i,j) = x0[shift+i];
    }
  }

  // 2. create and fill with data
  mvec_t b(contigMap_, 3);
  auto bh = b.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  std::vector<ST> b0(numGlobalEntries_);
  for (int i=0; i<numGlobalEntries_; ++i) b0[i]= (double) i + 2.;
  for (int i=0; i<localSize_; ++i) {
    for (int j=0; j<3; ++j){
      bh(i,j) = b0[shift+i];
    }
  }

  // 3. compute ewise multiply
  const auto alpha1 = 1.;
  vec_t c(contigMap_);
  const double offset = -5.0;
  pressio::ops::fill(c, offset);
  auto x_expr = pressio::column(x, 0);
  auto b_expr = pressio::column(b, 0);
  pressio::ops::elementwise_multiply(alpha1, x_expr, b_expr, 1., c);

  auto ch = c.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  vec_t c0(contigMap_);
  pressio::ops::fill(c0, std::nan("0")); // test beta=0 with simulated NaN in uninitialized c
  pressio::ops::elementwise_multiply(alpha1, x_expr, b_expr, 0., c0);
  auto ch0 = c0.getLocalViewHost(Tpetra::Access::ReadWriteStruct());

  // 5. check correctness
  for (int i=0; i<localSize_; ++i){
    const auto gold = alpha1 * xh(i, 0) * bh(i, 0);
    EXPECT_DOUBLE_EQ(ch0(i,0), gold);
    EXPECT_DOUBLE_EQ(ch(i,0), gold + offset);
  }
}
