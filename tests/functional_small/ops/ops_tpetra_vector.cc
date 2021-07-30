
#include "tpetra_only_fixtures.hpp"
#include "pressio_ops.hpp"

TEST_F(tpetraVectorGlobSize15Fixture, vector_clone)
{
    auto a = pressio::ops::clone(*myVector_);
    ASSERT_TRUE(a.getLocalViewDevice().data() != myVector_->getLocalViewDevice().data());

    auto a_h = a.getLocalViewHost();
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(a_h(i,0), 0.0);
    }

    myVector_->putScalar(23.);
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(a_h(i,0), 0.0);
    }
}

TEST_F(tpetraVectorGlobSize15Fixture, vector_extent)
{
    ASSERT_TRUE(pressio::ops::extent(*myVector_,0) == 15);
}

TEST_F(tpetraVectorGlobSize15Fixture, vector_deep_copy)
{
    pressio::ops::fill(*myVector_, -5.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::deep_copy(a, *myVector_);

    auto a_h = a.getLocalViewHost();
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(a_h(i,0), -5.);
    }
}

TEST_F(tpetraVectorGlobSize15Fixture, vector_setzero)
{
    myVector_->putScalar(23.);
    auto x_h = myVector_->getLocalViewHost();
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(x_h(i,0), 23.);
    }

    pressio::ops::set_zero(*myVector_);
    auto x_h2 = myVector_->getLocalViewHost();
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(x_h2(i,0), 0.);
    }
}

TEST_F(tpetraVectorGlobSize15Fixture, vector_fill)
{
    pressio::ops::fill(*myVector_, 55.);
    auto x_h = myVector_->getLocalViewHost();
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(x_h(i,0), 55.);
    }
}

TEST_F(tpetraVectorGlobSize15Fixture, vector_abs)
{
    pressio::ops::fill(*myVector_, -5.);
    auto x_h = myVector_->getLocalViewHost();
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(x_h(i,0), -5.);
    }

    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::abs(a, *myVector_);
    auto a_h = a.getLocalViewHost();
    for (int i=0; i<localSize_; ++i){
        EXPECT_DOUBLE_EQ(a_h(i,0), 5.);
        EXPECT_DOUBLE_EQ(x_h(i,0), -5.);
    }
}

TEST_F(tpetraVectorGlobSize15Fixture, vector_dot)
{
  auto a = pressio::ops::clone(*myVector_);
  auto b = pressio::ops::clone(*myVector_);
  pressio::ops::fill(a, 1.0);
  pressio::ops::fill(b, 1.0);

  auto res = ::pressio::ops::dot(a, b);
  EXPECT_DOUBLE_EQ(res, 15.0);

  res = 0.0;
  ::pressio::ops::dot(a, b, res);
  EXPECT_DOUBLE_EQ(res, 15.0);
}

TEST_F(tpetraVectorGlobSize15Fixture, vector_norm2)
{
  pressio::ops::fill(*myVector_, 1.0);
  auto mynorm = pressio::ops::norm2(*myVector_);
  EXPECT_DOUBLE_EQ(mynorm, std::sqrt(15.0));
}

TEST_F(tpetraVectorGlobSize15Fixture, vector_norm1)
{
  pressio::ops::fill(*myVector_, 1.0);
  auto mynorm = pressio::ops::norm1(*myVector_);
  EXPECT_DOUBLE_EQ(mynorm, 15.0);
}

TEST(ops_tpetra, vector_pow)
{
  {
    using tcomm = Teuchos::Comm<int>;
    using map_t = Tpetra::Map<>;
    using vec_t = Tpetra::Vector<>;
    using ST = typename vec_t::scalar_type;

    Teuchos::RCP<const tcomm> comm =
      Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    auto rank    = comm->getRank();
    auto numProc = comm->getSize();
    EXPECT_EQ(numProc,3);

    int numLocalEntries = 2;
    auto numGlobalEntries = numProc * numLocalEntries;
    Teuchos::RCP<const map_t> map =
      Teuchos::rcp(new map_t(numGlobalEntries, 0, comm));

    auto shift = rank*numLocalEntries;

    // 1. create vector x and fill with data
    vec_t x(map);
    auto xh = x.getLocalViewHost();
    std::vector<ST> x0(numGlobalEntries);
    for (int i=0; i<6; ++i) x0[i]= (double) i;
    for (int i=0; i<numLocalEntries; ++i) xh(i,0) = x0[shift+i];

    pressio::ops::pow(x, 2.);

    // 5. check correctness
    Eigen::VectorXd g(numGlobalEntries);
    g(0) = 0.;
    g(1) = 1.;
    g(2) = 4.;
    g(3) = 9.;
    g(4) = 16.;
    g(5) = 25.;
    ASSERT_EQ( pressio::ops::extent(x,0), 6 );
    for (int i=0; i<numLocalEntries; ++i){
      EXPECT_DOUBLE_EQ(xh(i,0), g(shift+i));
    }
  }
}

TEST(ops_tpetra, vector_absPowPos)
{
  {
    using tcomm = Teuchos::Comm<int>;
    using map_t = Tpetra::Map<>;
    using vec_t = Tpetra::Vector<>;
    using ST = typename vec_t::scalar_type;

    Teuchos::RCP<const tcomm> comm =
      Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    auto rank    = comm->getRank();
    auto numProc = comm->getSize();
    EXPECT_EQ(numProc,3);

    int numLocalEntries = 2;
    auto numGlobalEntries = numProc * numLocalEntries;
    Teuchos::RCP<const map_t> map =
      Teuchos::rcp(new map_t(numGlobalEntries, 0, comm));

    auto shift = rank*numLocalEntries;

    // 1. create vector x and fill with data
    vec_t y(map);
    auto yh = y.getLocalViewHost();
    vec_t x(map);
    auto xh = x.getLocalViewHost();
    std::vector<ST> x0(numGlobalEntries);
    for (int i=0; i<6; ++i){
      x0[i] = (double) i;
      x0[i] = -x0[i];
    }
    for (int i=0; i<numLocalEntries; ++i) xh(i,0) = x0[shift+i];

    // 2. compute
    ::pressio::ops::abs_pow(y, x, 3.);

    // 5. check correctness
    Eigen::VectorXd g(numGlobalEntries);
    g(0) = 0.;
    g(1) = 1.;
    g(2) = 8.;
    g(3) = 27.;
    g(4) = 64.;
    g(5) = 125.;
    for (int i=0; i<numLocalEntries; ++i){
      EXPECT_DOUBLE_EQ(yh(i,0), g(shift+i));
    }
  }
}

TEST(ops_tpetra, vector_absPowNeg)
{
  {
    using tcomm = Teuchos::Comm<int>;
    using map_t = Tpetra::Map<>;
    using vec_t = Tpetra::Vector<>;
    using ST = typename vec_t::scalar_type;

    Teuchos::RCP<const tcomm> comm =
      Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    auto rank    = comm->getRank();
    auto numProc = comm->getSize();
    EXPECT_EQ(numProc,3);

    int numLocalEntries = 2;
    auto numGlobalEntries = numProc * numLocalEntries;
    Teuchos::RCP<const map_t> map =
      Teuchos::rcp(new map_t(numGlobalEntries, 0, comm));

    auto shift = rank*numLocalEntries;

    // 1. create vector x and fill with data
    vec_t y(map);
    auto yh = y.getLocalViewHost();
    vec_t x(map);
    auto xh = x.getLocalViewHost();
    std::vector<ST> x0(numGlobalEntries);
    for (int i=0; i<6; ++i){
      x0[i] = (double) i;
      x0[i] = -x0[i];
    }
    for (int i=0; i<numLocalEntries; ++i) xh(i,0) = x0[shift+i];

    // 2. compute
    pressio::ops::abs_pow(y, x, -3., 0.001);

    // 5. check correctness
    Eigen::VectorXd g(numGlobalEntries);
    g(0) = 1./0.001; //because we guard against div by 0
    g(1) = 1.;
    g(2) = 1./8.;
    g(3) = 1./27.;
    g(4) = 1./64.;
    g(5) = 1./125.;
    for (int i=0; i<numLocalEntries; ++i){
      EXPECT_DOUBLE_EQ(yh(i,0), g(shift+i));
    }
  }
}


TEST_F(tpetraVectorGlobSize15Fixture, vector_update1_a)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    pressio::ops::update(v, a, 1.);
    auto v_h = v.getLocalViewHost();
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 2.);
    }
}

TEST_F(tpetraVectorGlobSize15Fixture, vector_update1_b)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    pressio::ops::update(v, 1., a, 1.);
    auto v_h = v.getLocalViewHost();
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 3.);
    }
}

TEST_F(tpetraVectorGlobSize15Fixture, vector_update2_a)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    auto b = pressio::ops::clone(*myVector_);
    pressio::ops::fill(b, 3.);

    pressio::ops::update(v, a, 1., b, 1.);
    auto v_h = v.getLocalViewHost();
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 5.);
    }
}

TEST_F(tpetraVectorGlobSize15Fixture, vector_update2_b)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    auto b = pressio::ops::clone(*myVector_);
    pressio::ops::fill(b, 3.);

    pressio::ops::update(v, 1., a, 1., b, 1.);
    auto v_h = v.getLocalViewHost();
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 6.);
    }
}

TEST_F(tpetraVectorGlobSize15Fixture, vector_update3_a)
{
    auto v = pressio::ops::clone(*myVector_);
    pressio::ops::fill(v, 1.);
    auto a = pressio::ops::clone(*myVector_);
    pressio::ops::fill(a, 2.);
    auto b = pressio::ops::clone(*myVector_);
    pressio::ops::fill(b, 3.);
    auto c = pressio::ops::clone(*myVector_);
    pressio::ops::fill(c, 4.);

    pressio::ops::update(v, a, 1., b, 1., c, 1.);
    auto v_h = v.getLocalViewHost();
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 9.);
    }
}

TEST_F(tpetraVectorGlobSize15Fixture, vector_update3_b)
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
    auto v_h = v.getLocalViewHost();
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 10.);
    }
}

TEST_F(tpetraVectorGlobSize15Fixture, vector_update4_a)
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

    pressio::ops::update(v, a, 1., b, 1., c, 1., d, 1.);
    auto v_h = v.getLocalViewHost();
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 14.);
    }
}

TEST_F(tpetraVectorGlobSize15Fixture, vector_update4_b)
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
    auto v_h = v.getLocalViewHost();
    for (int i=0; i<localSize_; ++i){
      EXPECT_DOUBLE_EQ(v_h(i,0), 15.);
    }
}

TEST(ops_tpetra, vector_elementwiseMultiply)
{
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;
  using vec_t = Tpetra::Vector<>;
  using ST = typename vec_t::scalar_type;

  Teuchos::RCP<const tcomm> comm = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  auto rank    = comm->getRank();
  auto numProc = comm->getSize();
  EXPECT_EQ(numProc,3);

  int numLocalEntries = 2;
  auto numGlobalEntries = numProc * numLocalEntries;
  Teuchos::RCP<const map_t> map = Teuchos::rcp(new map_t(numGlobalEntries, 0, comm));

  auto shift = rank*numLocalEntries;

  // c = beta*c + alpha*x*b;

  // 1. create vector x and fill with data
  vec_t x(map);
  auto xh = x.getLocalViewHost();
  std::vector<ST> x0(numGlobalEntries);
  for (int i=0; i<6; ++i) x0[i]= (double) i;
  for (int i=0; i<numLocalEntries; ++i) xh(i,0) = x0[shift+i];

  // 2. create vector b and fill with data
  vec_t b(map);
  auto bh = b.getLocalViewHost();
  std::vector<ST> b0(numGlobalEntries);
  for (int i=0; i<6; ++i) b0[i]= (double) i + 2.;
  for (int i=0; i<numLocalEntries; ++i) bh(i,0) = b0[shift+i];

  // 3. compute ewise multiply
  vec_t c(map);
  pressio::ops::fill(c, 0.);
  pressio::ops::elementwise_multiply(1., x, b, 0., c);
  auto ch = c.getLocalViewHost();

  // 5. check correctness
  Eigen::VectorXd g(numGlobalEntries);
  g(0) = 0.;
  g(1) = 3.;
  g(2) = 8.;
  g(3) = 15.;
  g(4) = 24.;
  g(5) = 35.;

  for (int i=0; i<numLocalEntries; ++i){
    EXPECT_DOUBLE_EQ(ch(i,0), g(shift+i));
  }
}
