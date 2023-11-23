
#include "tpetra_block_only_fixtures.hpp"
#include "pressio/ops.hpp"

// convenient alias for nice test names
using ops_tpetra_block = tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture;

TEST_F(ops_tpetra_block, column_expr_extent)
{
  auto e = pressio::column(*myMv_, 0);
  ASSERT_TRUE(pressio::ops::extent(e,0) == numProc_ * 5.);
  ASSERT_TRUE(pressio::ops::extent(e,1) == 1); // check extent over the rank
}

TEST_F(ops_tpetra_block, column_expr_abs)
{
  const Tpetra::Access::ReadWriteStruct rw;

  myMv_->putScalar(-5);
  auto x_h = myMv_->getMultiVectorView().getLocalViewHost(rw);
  for (size_t k=0; k<x_h.extent(1); ++k){
    for (size_t i=0; i<x_h.extent(0); ++i){
      EXPECT_DOUBLE_EQ(x_h(i,k), -5.);
    }
  }

  for (int k=0; k<numVecs_; ++k){
    vec_t vtest(*contigMap_, blockSize_);
    vtest.putScalar(-11.);

    auto e = pressio::column(*myMv_, k);
    pressio::ops::abs(vtest, e);
    auto v_h = vtest.getVectorView().getLocalViewHost(rw);
    for (size_t i=0; i<v_h.extent(0); ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 5.);
    }
  }
}

TEST_F(ops_tpetra_block, column_expr_dot_with_vector)
{
  myMv_->putScalar(1.);
  auto a = pressio::column(*myMv_, 0);

  vec_t b(*contigMap_, blockSize_);
  b.putScalar(1.);

  {
    auto res = ::pressio::ops::dot(a, b);
    EXPECT_DOUBLE_EQ(res, numProc_ * 5. * blockSize_);
    res = 0.0;
    ::pressio::ops::dot(a, b, res);
    EXPECT_DOUBLE_EQ(res, numProc_ * 5. * blockSize_);
  }

  {
    auto res = ::pressio::ops::dot(b, a);
    EXPECT_DOUBLE_EQ(res, numProc_ * 5. * blockSize_);
    res = 0.0;
    ::pressio::ops::dot(a, b, res);
    EXPECT_DOUBLE_EQ(res, numProc_ * 5. * blockSize_);
  }
}

TEST_F(ops_tpetra_block, column_expr_dot_with_expr)
{
  myMv_->putScalar(1.);
  auto a = pressio::column(*myMv_, 0);
  auto b = pressio::column(*myMv_, 1);
  auto res = ::pressio::ops::dot(a, b);
  EXPECT_DOUBLE_EQ(res, numProc_ * 5. * blockSize_);
  res = 0.0;
  ::pressio::ops::dot(a, b, res);
  EXPECT_DOUBLE_EQ(res, numProc_ * 5. * blockSize_);
}

TEST_F(ops_tpetra_block, column_expr_setzero)
{
  const Tpetra::Access::ReadWriteStruct rw;

  myMv_->putScalar(-5);
  auto x_h = myMv_->getMultiVectorView().getLocalViewHost(rw);
  for (size_t k=0; k<x_h.extent(1); ++k){
    for (size_t i=0; i<x_h.extent(0); ++i){
      EXPECT_DOUBLE_EQ(x_h(i,k), -5.);
    }
  }

  for (int k=0; k<numVecs_; ++k){
    auto e = pressio::column(*myMv_, k);
    pressio::ops::set_zero(e);
    for (size_t i=0; i<x_h.extent(0); ++i){
      EXPECT_DOUBLE_EQ(x_h(i,0), 0.);
    }
  }
}

TEST_F(ops_tpetra_block, column_expr_fill)
{
  const Tpetra::Access::ReadWriteStruct rw;

  myMv_->putScalar(-5);
  auto x_h = myMv_->getMultiVectorView().getLocalViewHost(rw);
  for (size_t k=0; k<x_h.extent(1); ++k){
    for (size_t i=0; i<x_h.extent(0); ++i){
      EXPECT_DOUBLE_EQ(x_h(i,k), -5.);
    }
  }

  for (int k=0; k<numVecs_; ++k){
    auto e = pressio::column(*myMv_, k);
    pressio::ops::fill(e, 44.);
    for (size_t i=0; i<x_h.extent(0); ++i){
      EXPECT_DOUBLE_EQ(x_h(i,0), 44.);
    }
  }
}

TEST_F(ops_tpetra_block, column_expr_scale)
{
  const Tpetra::Access::ReadWriteStruct rw;

  myMv_->putScalar(-5);
  auto x_h = myMv_->getMultiVectorView().getLocalViewHost(rw);
  for (size_t k=0; k<x_h.extent(1); ++k){
    for (size_t i=0; i<x_h.extent(0); ++i){
      EXPECT_DOUBLE_EQ(x_h(i,k), -5.);
    }
  }

  for (int k=0; k<numVecs_; ++k){
    auto e = pressio::column(*myMv_, k);
    pressio::ops::scale(e, 2.);
    for (size_t i=0; i<x_h.extent(0); ++i){
      EXPECT_DOUBLE_EQ(x_h(i,0), -10.);
    }
  }
}

TEST_F(ops_tpetra_block, column_expr_minmax)
{
  for (int k=0; k<numVecs_; ++k){
    auto v_tp = myMv_->getMultiVectorView().getVectorNonConst(k);
    auto e = pressio::column(*myMv_, k);
    ASSERT_DOUBLE_EQ(pressio::ops::min(*v_tp), pressio::ops::min(e));
  }
}

TEST_F(ops_tpetra_block, column_expr_norm1)
{
  myMv_->putScalar(1);
  for (int k=0; k<numVecs_; ++k){
    auto e = pressio::column(*myMv_, k);
    ASSERT_DOUBLE_EQ(pressio::ops::norm1(e), localSize_*blockSize_*numProc_);
  }
}

TEST_F(ops_tpetra_block, column_expr_norm2)
{
  myMv_->putScalar(1);
  for (int k=0; k<numVecs_; ++k){
    auto e = pressio::column(*myMv_, k);
    ASSERT_DOUBLE_EQ(pressio::ops::norm2(e),
		     std::sqrt(localSize_*blockSize_*numProc_));
  }
}

TEST_F(ops_tpetra_block, column_expr_pow)
{
  const Tpetra::Access::ReadOnlyStruct ro;

  myMv_->putScalar(2.);
  for (int k=0; k<numVecs_; ++k){
    auto e = pressio::column(*myMv_, k);
    pressio::ops::pow(e, 2.);
  }

  auto x_h = myMv_->getMultiVectorView().getLocalViewHost(ro);
  for (size_t k=0; k<x_h.extent(1); ++k){
    for (size_t i=0; i<x_h.extent(0); ++i){
      EXPECT_DOUBLE_EQ(x_h(i,k), 4.);
    }
  }
}

TEST_F(ops_tpetra_block, column_expr_abs_pow_pos_expo)
{
  const Tpetra::Access::ReadOnlyStruct ro;

  myMv_->putScalar(-2.);
  auto a = pressio::column(*myMv_, 0);
  auto b = pressio::column(*myMv_, 1);
  pressio::ops::abs_pow(b, a, 3.);

  auto x_h = myMv_->getMultiVectorView().getLocalViewHost(ro);
  for (size_t k=0; k<x_h.extent(1); ++k){
    for (size_t i=0; i<x_h.extent(0); ++i){
      if (k==1){
	EXPECT_DOUBLE_EQ(x_h(i,k), 8.);
      }
      else{
	EXPECT_DOUBLE_EQ(x_h(i,k), -2.);
      }
    }
  }
}

TEST_F(ops_tpetra_block, column_expr_abs_pow_neg_expo)
{
  const Tpetra::Access::ReadOnlyStruct ro;

  myMv_->putScalar(-2.);
  auto a = pressio::column(*myMv_, 0);
  auto b = pressio::column(*myMv_, 1);
  pressio::ops::abs_pow(b, a, -3., 1e-6);

  auto x_h = myMv_->getMultiVectorView().getLocalViewHost(ro);
  for (size_t k=0; k<x_h.extent(1); ++k){
    for (size_t i=0; i<x_h.extent(0); ++i){
      if (k==1){
	EXPECT_DOUBLE_EQ(x_h(i,k), 1./8.);
      }
      else{
	EXPECT_DOUBLE_EQ(x_h(i,k), -2.);
      }
    }
  }
}

TEST_F(ops_tpetra_block, column_expr_update1_a)
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

TEST_F(ops_tpetra_block, column_expr_update1_b)
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

TEST_F(ops_tpetra_block, column_expr_update2_a)
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

TEST_F(ops_tpetra_block, column_expr_update2_b)
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

TEST_F(ops_tpetra_block, column_expr_update3_a)
{
    auto v = pressio::column(*myMv_, 0);
    pressio::ops::fill(v, 1.);
    auto a = pressio::column(*myMv_, 1);
    pressio::ops::fill(a, 2.);
    auto b = pressio::column(*myMv_, 2);
    pressio::ops::fill(b, 3.);
    vec_t c(*contigMap_, blockSize_);
    pressio::ops::fill(c, 4.);

    pressio::ops::update(v, 0., a, 1., b, 1., c, 1.);
    auto v_h = v.native().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 9.);
    }
}

TEST_F(ops_tpetra_block, column_expr_update3_b)
{
    auto v = pressio::column(*myMv_, 0);
    pressio::ops::fill(v, 1.);
    auto a = pressio::column(*myMv_, 1);
    pressio::ops::fill(a, 2.);
    auto b = pressio::column(*myMv_, 2);
    pressio::ops::fill(b, 3.);
    vec_t c(*contigMap_, blockSize_);
    pressio::ops::fill(c, 4.);

    pressio::ops::update(v, 1., a, 1., b, 1., c, 1.);
    auto v_h = v.native().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 10.);
    }
}

TEST_F(ops_tpetra_block, column_expr_elementwiseMultiply)
{
  // computing elementwise:  y = beta * y + alpha * x * z

  auto x = pressio::column(*myMv_, 0);
  pressio::ops::fill(x, 3.);

  auto z = pressio::column(*myMv_, 1);
  pressio::ops::fill(z, 2.);

  auto y = pressio::column(*myMv_, 2);
  pressio::ops::fill(y, 2.);

  pressio::ops::elementwise_multiply(2., x, z, 3., z);

  auto z_h = z.native().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    EXPECT_DOUBLE_EQ(z_h(i,0), 18.);
  }

  // test beta=0 with simulated NaN in uninitialized y
  pressio::ops::fill(x, 3.);
  pressio::ops::fill(z, 2.);
  pressio::ops::fill(y, std::nan("0"));
  pressio::ops::elementwise_multiply(2., x, z, 0., y);
  auto y_h = y.native().getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  for (int i=0; i<localSize_*blockSize_; ++i){
    EXPECT_DOUBLE_EQ(y_h(i,0), 12);
  }
}
