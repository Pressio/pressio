
#include "tpetra_block_only_fixtures.hpp"
#include "pressio/ops.hpp"
#include "ops_shared_level2.hpp"

//-------------------------------------------
// Test implementation and utilities
//-------------------------------------------

using vec_t = typename tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture::vec_t;
using mvec_t = typename tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture::mvec_t;

// convenient alias for nicer test names
using ops_tpetra_block = tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture;

template <typename MultiVectorType, typename CommType>
auto create_importer(const MultiVectorType &mv_src, const CommType &comm) {
  using map_t = typename MultiVectorType::map_type;
  using sc_t = typename MultiVectorType::scalar_type;
  using lo_t = typename MultiVectorType::local_ordinal_type;
  using go_t = typename MultiVectorType::global_ordinal_type;
  using node_t = typename MultiVectorType::node_type;
  using mvec_t = Tpetra::MultiVector<sc_t, lo_t, go_t, node_t>;
  using importer_t = Tpetra::Import<lo_t, go_t, node_t>;

  const auto numEntries = mv_src.getGlobalLength();
  const auto map0 = Teuchos::rcp(new map_t(numEntries, comm->getRank() == 0 ? numEntries : 0, 0, comm));
  return importer_t(mv_src.getMap(), map0);
}

// returns local copy containing all global data
template <typename T>
auto get_global_host_view(T &view, const tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture &test) {
  // need a deep copy here and pressio::ops::clone(view)
  // won't produce it for Eigen/Kokkos expressions
  auto size = ::pressio::ops::extent(view, 0);
  Eigen::VectorXd v_local(size);
  for (decltype(size) i = 0; i < size; ++i) {
    v_local(i) = view(i);
  }
  return v_local;
}

template <>
auto get_global_host_view<mvec_t>(mvec_t &mv, const tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture &test) {
  auto mv_src = mv.getMultiVectorView(); // cast to non-block view for importing
  auto importer = create_importer(mv_src, test.comm_);
  decltype(mv_src) mv_import(importer.getTargetMap(), mv.getNumVectors());
  mv_import.doImport(mv_src, importer, Tpetra::REPLACE);
  return mv_import.getLocalViewHost(Tpetra::Access::ReadWrite);
}

template <>
auto get_global_host_view<vec_t>(vec_t &v, const tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture &test) {
  auto v_src = v.getVectorView(); // cast to non-block view for importing
  auto importer = create_importer(v_src, test.comm_);
  decltype(v_src) v_import(importer.getTargetMap());
  v_import.doImport(v_src, importer, Tpetra::REPLACE);
  auto view2d = v_import.getLocalViewHost(Tpetra::Access::ReadWrite);
  return Kokkos::subview(view2d, Kokkos::ALL(), 0);
}

template <
    typename FixtureType,
    typename TransMode,
    typename AType,
    typename XType,
    typename YType,
    typename ScalarType
    >
void test_impl(FixtureType &test, TransMode trans, ScalarType alpha, AType A, XType x, ScalarType beta, YType &y) {
    Kokkos::Profiling::pushRegion("test_impl");
    // copy original values, fetch whole vector (from all ranks) to do the verification locally
    auto y_ref_h = get_global_host_view(y, test);

    // call tested routine on device
    ::pressio::ops::product(trans, alpha, A, x, beta, y);

    // Note: trick here is to fetch whole data (from all ranks)
    // and run the simplified reference routine locally
    auto A_h = get_global_host_view(A, test);
    auto x_h = get_global_host_view(x, test);
    auto y_h = get_global_host_view(y, test);
    if (test.rank_ == 0) {
        vanilla_gemv(trans, alpha, A_h, x_h, beta, y_ref_h);
        const size_t y_size = ::pressio::ops::extent(y_h, 0);
        for (size_t i = 0; i < y_size; ++i) {
          EXPECT_DOUBLE_EQ(y_h(i), y_ref_h(i));
        }
    }
    Kokkos::Profiling::popRegion();
}

template <
    typename FixtureType,
    typename TransMode,
    typename AType,
    typename XType,
    typename YType
    >
void test_impl(FixtureType &test, TransMode trans, AType A, XType x, YType y) {
    // non-zero alpha and beta
    test_impl(test, trans, 2., A, x, -1., y);
    // simulate beta=0 with uninitialized y containing NaN
    const auto nan = std::nan("0");
    ::pressio::ops::fill(y, nan);
    test_impl(test, trans, 2., A, x, 0., y);
    // simulate alpha=0 with NaN in input matrix
    auto A_h = A.getMultiVectorView().getLocalViewHost(Tpetra::Access::ReadWrite);
    auto a00 = A_h(0, 0);
    A_h(0, 0) = nan;
    test_impl(test, trans, 0., A, x, 0., y);
    // alpha=0, beta=1
    ::pressio::ops::fill(y, 2.);
    test_impl(test, trans, 0., A, x, 1., y);
    // restore original value
    A_h(0, 0) = a00;
}

//-------------------------------------------
// Test Kokkos x
//-------------------------------------------

TEST_F(ops_tpetra_block,
       mv_prod_kokkos_vector)
{
  Kokkos::View<double*> x_kokkos{"x", (size_t)numVecs_};
  auto x_h = Kokkos::create_mirror_view(Kokkos::HostSpace(), x_kokkos);
  for (int j = 0; j < numVecs_; ++j) {
    x_h(j) = (double)(j + 1.); // unique int values
  }
  Kokkos::deep_copy(x_kokkos, x_h);

  test_impl(*this, ::pressio::nontranspose{}, *myMv_, x_kokkos, *y_tpetra);
}

TEST_F(ops_tpetra_block,
       mv_prod_kokkos_span)
{
  Kokkos::View<double*> x0{ "x_span", numVecs_ + (size_t)2 };
  auto x_h = Kokkos::create_mirror_view(Kokkos::HostSpace(), x0);
  for (int i = 0; i < numVecs_ + 2; ++i) {
    x_h(i) = (double)(i + 1.); // unique int values
  }
  Kokkos::deep_copy(x0, x_h);
  auto x_kokkos_span = ::pressio::span(x0, 1, numVecs_);

  test_impl(*this, ::pressio::nontranspose{}, *myMv_, x_kokkos_span, *y_tpetra);
}

TEST_F(ops_tpetra_block,
       mv_prod_kokkos_diag)
{
  Kokkos::View<double**> x0{ "x_diag", (size_t)numVecs_, (size_t)numVecs_ };
  auto x_h = Kokkos::create_mirror_view(Kokkos::HostSpace(), x0);
  for (size_t i = 0; i < numVecs_; ++i) {
      for (size_t j = 0; j < numVecs_; ++j) {
          x_h(i, j) = (double)(i * numVecs_ + j + 1.0);
      }
  }
  Kokkos::deep_copy(x0, x_h);
  auto x_kokkos_diag = ::pressio::diag(x0);

  test_impl(*this, ::pressio::nontranspose{}, *myMv_, x_kokkos_diag, *y_tpetra);
}

//-------------------------------------------
// Test Kokkos y
//-------------------------------------------

TEST_F(ops_tpetra_block,
       mv_T_vector_storein_kokkos_vector)
{
  Kokkos::View<double*> y_kokkos("a", numVecs_);
  KokkosBlas::fill(y_kokkos, 1.);

  test_impl(*this, ::pressio::transpose{}, *myMv_, *x_tpetra, y_kokkos);
}

TEST_F(ops_tpetra_block,
       mv_T_vector_storein_kokkos_span)
{
    Kokkos::View<double*> y0{ "y_span", numVecs_ + 2 };
    auto y_h = Kokkos::create_mirror_view(Kokkos::HostSpace(), y0);
    for (int i = 0; i < numVecs_ + 2; ++i) {
      y_h(i) = (double)(i + 1.); // unique int values
    }
    Kokkos::deep_copy(y0, y_h);
    auto y_kokkos_span = ::pressio::span(y0, 1, numVecs_);

    test_impl(*this, ::pressio::transpose{}, *myMv_, *x_tpetra, y_kokkos_span);
}

TEST_F(ops_tpetra_block,
       mv_T_vector_storein_kokkos_diag)
{
    Kokkos::View<double**> y0{ "y_diag", numVecs_, numVecs_ };
    auto y_h = Kokkos::create_mirror_view(Kokkos::HostSpace(), y0);
    for (size_t i = 0; i < numVecs_; ++i) {
        for (size_t j = 0; j < numVecs_; ++j) {
            y_h(i, j) = (double)(i * numVecs_ + j + 1.0);
        }
    }
    Kokkos::deep_copy(y0, y_h);
    auto y_kokkos_diag = ::pressio::diag(y0);

    test_impl(*this, ::pressio::transpose{}, *myMv_, *x_tpetra, y_kokkos_diag);
}

#ifdef PRESSIO_ENABLE_TPL_EIGEN

//-------------------------------------------
// Test Eigen x
//-------------------------------------------

TEST_F(ops_tpetra_block,
       mv_prod_eigen_vector)
{
  Eigen::VectorXd x_eigen(numVecs_);
  x_eigen.setConstant(1.);

  test_impl(*this, ::pressio::nontranspose{}, *myMv_, x_eigen, *y_tpetra);
}

TEST_F(ops_tpetra_block,
       mv_prod_eigen_span)
{
  Eigen::VectorXd x0(numVecs_ + 3);
  x0.setConstant(1.);
  auto x_eigen_span = ::pressio::span(x0, 2, numVecs_);

  test_impl(*this, ::pressio::nontranspose{}, *myMv_, x_eigen_span, *y_tpetra);
}

TEST_F(ops_tpetra_block,
       mv_prod_eigen_diag)
{
  Eigen::MatrixXd M0(numVecs_, numVecs_);
  for (int i = 0; i < numVecs_; ++i) {
      M0(i, i) = 1.;
  }
  auto x_eigen_diag = ::pressio::diag(M0);

  test_impl(*this, ::pressio::nontranspose{}, *myMv_, x_eigen_diag, *y_tpetra);
}

//-------------------------------------------
// Test Eigen y
//-------------------------------------------

TEST_F(ops_tpetra_block,
       mv_T_vector_storein_eigen_vector)
{
  Eigen::VectorXd y_eigen(numVecs_);
  y_eigen.setConstant(1.);

  test_impl(*this, ::pressio::transpose{}, *myMv_, *x_tpetra, y_eigen);
}

TEST_F(ops_tpetra_block,
       mv_T_vector_storein_eigen_span)
{
  Eigen::VectorXd y0(numVecs_ + 3);
  y0.setConstant(1.);
  auto y_eigen_span = ::pressio::span(y0, 2, numVecs_);

  test_impl(*this, ::pressio::transpose{}, *myMv_, *x_tpetra, y_eigen_span);
}

TEST_F(ops_tpetra_block,
       mv_T_vector_storein_eigen_diag)
{
  Eigen::MatrixXd M0(numVecs_, numVecs_);
  for (int i = 0; i < numVecs_; ++i) {
      M0(i, i) = 1.;
  }
  auto y_eigen_diag = ::pressio::diag(M0);

  test_impl(*this, ::pressio::transpose{}, *myMv_, *x_tpetra, y_eigen_diag);
}
#endif
