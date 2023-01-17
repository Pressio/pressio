
#include "epetra_only_fixtures.hpp"
#include "pressio/ops.hpp"

// convenient alias for nice test names
using ops_epetra = epetraVectorGlobSize15Fixture;

TEST_F(ops_epetra, vector_clone)
{
    auto a = pressio::ops::clone(*myVector_);
    ASSERT_TRUE(a.Values() != myVector_->Values());

    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(a[i], 0.0);
    }

    myVector_->PutScalar(23.);
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(a[i], 0.0);
    }
}

TEST_F(ops_epetra, vector_extent)
{
    ASSERT_TRUE(pressio::ops::extent(*myVector_,0) == numProc_ * 5.);
}

TEST_F(ops_epetra, vector_deep_copy)
{
    myVector_->PutScalar(-5.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::deep_copy(a, *myVector_);

    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(a[i], -5.);
    }
}

TEST_F(ops_epetra, vector_setzero)
{
    myVector_->PutScalar(23.);
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ((*myVector_)[i], 23.);
    }

    pressio::ops::set_zero(*myVector_);
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ((*myVector_)[i], 0.);
    }
}

TEST_F(ops_epetra, vector_fill)
{
    pressio::ops::fill(*myVector_, 55.);
    auto & x_h = *myVector_;
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(x_h[i], 55.);
    }
}

TEST_F(ops_epetra, vector_abs)
{
    pressio::ops::fill(*myVector_, -5.);
    auto & x_h = *myVector_;
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(x_h[i], -5.);
    }

    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::abs(a, *myVector_);
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(a[i], 5.);
        EXPECT_DOUBLE_EQ(x_h[i], -5.);
    }
}

TEST_F(ops_epetra, vector_dot)
{
  auto a = pressio::ops::clone(*myVector_);
  auto b = pressio::ops::clone(*myVector_);
  a.PutScalar(1.0);
  b.PutScalar(1.0);

  auto res = ::pressio::ops::dot(a, b);
  EXPECT_DOUBLE_EQ(res, numProc_ * 5.);

  res = 0.0;
  ::pressio::ops::dot(a, b, res);
  EXPECT_DOUBLE_EQ(res, numProc_ * 5.);
}

TEST_F(ops_epetra, vector_norm2)
{
  myVector_->PutScalar(1.0);
  auto mynorm = pressio::ops::norm2(*myVector_);
  EXPECT_NEAR(mynorm, std::sqrt(numProc_ * 5.0), 1e-15);
}

TEST_F(ops_epetra, vector_norm1)
{
  myVector_->PutScalar(1.0);
  auto mynorm = pressio::ops::norm1(*myVector_);
  EXPECT_DOUBLE_EQ(mynorm, numProc_ * 5.0);
}

TEST_F(ops_epetra, vector_pow)
{
  auto & x = *myVector_;
  x.PutScalar(2.);
  pressio::ops::pow(x, 2.);
  for (int i=0; i<localSize_; ++i){
    EXPECT_DOUBLE_EQ(x[i], 4.);
  }
}

TEST_F(ops_epetra, vector_absPowPos)
{
  auto & x = *myVector_;
  x.PutScalar(-2.);

  auto y = pressio::ops::clone(x);
  pressio::ops::abs_pow(y, x, 3.);
  for (int i=0; i<localSize_; ++i){
    EXPECT_DOUBLE_EQ(y[i], 8.);
  }
}
// TEST(ops_epetra, vector_absPowNeg)
// {
//   {
//     using tcomm = Teuchos::Comm<int>;
//     using map_t = Tpetra::Map<>;
//     using vec_t = Tpetra::Vector<>;
//     using ST = typename vec_t::scalar_type;

//     Teuchos::RCP<const tcomm> comm =
//       Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
//     auto rank    = comm->getRank();
//     auto numProc = comm->getSize();
//     EXPECT_EQ(numProc,3);

//     int numLocalEntries = 2;
//     auto numGlobalEntries = numProc * numLocalEntries;
//     Teuchos::RCP<const map_t> map =
//       Teuchos::rcp(new map_t(numGlobalEntries, 0, comm));

//     auto shift = rank*numLocalEntries;

//     // 1. create vector x and fill with data
//     vec_t y(map);
//     auto yh = y.getLocalViewHost();
//     vec_t x(map);
//     auto xh = x.getLocalViewHost();
//     std::vector<ST> x0(numGlobalEntries);
//     for (int i=0; i<6; ++i){
//       x0[i] = (double) i;
//       x0[i] = -x0[i];
//     }
//     for (int i=0; i<numLocalEntries; ++i) xh(i,0) = x0[shift+i];

//     // 2. compute
//     pressio::ops::abs_pow(y, x, -3., 0.001);

//     // 5. check correctness
//     Eigen::VectorXd g(numGlobalEntries);
//     g(0) = 1./0.001; //because we guard against div by 0
//     g(1) = 1.;
//     g(2) = 1./8.;
//     g(3) = 1./27.;
//     g(4) = 1./64.;
//     g(5) = 1./125.;
//     for (int i=0; i<numLocalEntries; ++i){
//       EXPECT_DOUBLE_EQ(yh(i,0), g(shift+i));
//     }
//   }
// }

TEST_F(ops_epetra, vector_update1_a)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    pressio::ops::update(v, 0., a, 1.);
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v[i], 2.);
    }
}

TEST_F(ops_epetra, vector_update1_b)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    pressio::ops::update(v, 1., a, 1.);
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v[i], 3.);
    }
}

TEST_F(ops_epetra, vector_update2_a)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    auto b = pressio::ops::clone(*myVector_);
    pressio::ops::fill(b, 3.);

    pressio::ops::update(v, 0., a, 1., b, 1.);
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v[i], 5.);
    }
}

TEST_F(ops_epetra, vector_update2_b)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    auto b = pressio::ops::clone(*myVector_);
    pressio::ops::fill(b, 3.);

    pressio::ops::update(v, 1., a, 1., b, 1.);
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v[i], 6.);
    }
}

TEST_F(ops_epetra, vector_update3_a)
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
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v[i], 9.);
    }
}

TEST_F(ops_epetra, vector_update3_b)
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
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v[i], 10.);
    }
}

TEST_F(ops_epetra, vector_update4_a)
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
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v[i], 14.);
    }
}

TEST_F(ops_epetra, vector_update4_b)
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
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v[i], 15.);
    }
}

TEST_F(ops_epetra, vector_update4_c)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 3.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, std::nan("0"));

    pressio::ops::update(v, 2., a, 0., a, 0., a, 0., a, 0.);
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v[i], 6.);
    }
}

TEST_F(ops_epetra, vector_update_nan)
{
    const auto nan = std::nan("0");
    auto v = pressio::ops::clone(*myVector_);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2);

    pressio::ops::fill(v, nan);
    pressio::ops::update(v, 0., a, 1.);
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v[i], 2.);
    }

    pressio::ops::fill(v, nan);
    pressio::ops::update(v, 0., a, 1., a, 2.);
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v[i], 6.);
    }

    pressio::ops::fill(v, nan);
    pressio::ops::update(v, 0., a, 1., a, 2., a, 3.);
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v[i], 12.);
    }

    pressio::ops::fill(v, nan);
    pressio::ops::update(v, 0., a, 1., a, 2., a, 3., a, 4.);
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v[i], 20.);
    }
}

TEST_F(ops_epetra, vector_elementwiseMultiply)
{
    auto y = pressio::ops::clone(*myVector_);
    pressio::ops::fill(y, 1.);
    auto x = pressio::ops::clone(*myVector_);
    pressio::ops::fill(x, 2.);
    auto z = pressio::ops::clone(*myVector_);
    pressio::ops::fill(z, 3.);

    // y = beta*y + alpha*x*z;
    pressio::ops::elementwise_multiply(1., x,z, 2., y);
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(y[i], 8.);
    }

    // test beta=0 with simulated NaN in uninitialized y
    y[0] = std::nan("0");
    pressio::ops::elementwise_multiply(1., x,z, 0., y);
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(y[i], 6.);
    }
}
