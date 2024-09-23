
#include <gtest/gtest.h>
#include "pressio/ops.hpp"
#include "ops_shared_level2.hpp"

//-------------------------------------------
// Test implementation and utilities
//-------------------------------------------

// Note: Thanks to this wrapper we can create "double expressions"
// and pass them to testing routines same as Kokkos::DualView.
template <typename DualViewType, typename ExprGen>
class Expression2DualViewAdapter {
public:
  using dualview_type    = DualViewType;
  using device_view_type = typename dualview_type::t_dev;
  using host_view_type   = typename dualview_type::t_host;
  using device_expr_type = decltype((std::declval<ExprGen>())(std::declval<device_view_type&>()));
  using host_expr_type   = decltype((std::declval<ExprGen>())(std::declval<host_view_type&>()));

public:
  Expression2DualViewAdapter(DualViewType dual_view, ExprGen F)
    : base_view(dual_view),
      dev_view(base_view.view_device()),
      host_view(base_view.view_host()),
      // Important: expressions store references so we pass base views
      // that are stored in the adapter (as opposed to temp instances)
      expr_dev(F(dev_view)),
      expr_host(F(host_view))
    {}

  auto extent(size_t i) const { return expr_dev.extent(i); }

  auto& view_host()   { return expr_host; } // for modification on host
  auto& view_device() { return expr_dev; }  // for product() on device

  void sync_host()     { base_view.sync_host(); };
  void sync_device()   { base_view.sync_device(); };
  void modify_host()   { base_view.modify_host(); };
  void modify_device() { base_view.modify_device(); };

protected:
  dualview_type base_view;
  device_view_type dev_view;
  host_view_type host_view;
  device_expr_type expr_dev;
  host_expr_type expr_host;
};

template <typename DView, typename EGen>
static auto make_adapter(DView dual_view, EGen expr_gen) {
  return Expression2DualViewAdapter<DView, EGen>(dual_view, expr_gen);
}

template <typename DView>
auto span_adapter(DView dual_view, std::size_t index, std::size_t size) {
  return make_adapter(dual_view, [index, size](auto view) {
    return pressio::span(view, index, size);
  });
}

template <typename DView2D>
auto diag_adapter(DView2D matrix) {
  return make_adapter(matrix, [](auto mtx) {
    return pressio::diagonal(mtx);
  });
}

template <typename DView2D>
auto subspan_adapter(DView2D matrix, std::size_t i0, std::size_t ext0, std::size_t i1, std::size_t ext1) {
  using range_t = std::pair<size_t, size_t>;
  return make_adapter(matrix, [=](auto mtx) {
    return pressio::subspan(mtx, range_t{ i0, i0 + ext0 }, range_t{ i1, i1 + ext1 });
  });
}

struct kokkosFixture
  : public ::testing::Test {

  const double          NaN    = std::nan("0");
  static constexpr auto alpha0 = static_cast<double>(0);
  static constexpr auto alpha1 = static_cast<double>(1);
  static constexpr auto beta0  = alpha0;
  static constexpr auto beta1  = alpha1;

  const size_t x_size = 3;
  const size_t y_size = 4;
  // plain views
  Kokkos::DualView<double**> A{ "A", y_size, x_size };
  Kokkos::DualView<double*> x{ "x", x_size };
  Kokkos::DualView<double*> y{ "y", y_size };
  Kokkos::DualView<double*> xt{ "xt", y_size };
  Kokkos::DualView<double*> yt{ "yt", x_size };
  // expression base (data views)
  const size_t input_size_ext = (x_size > y_size ? x_size : y_size) + 2;
  Kokkos::DualView<double*> x_span_base{ "x_span",  input_size_ext };
  Kokkos::DualView<double*> y_span_base{ "y_span", input_size_ext };
  Kokkos::DualView<double**> x_diag_base{ "x_diag", x_size, x_size };
  Kokkos::DualView<double**> xt_diag_base{ "xt_diag", y_size, y_size };
  Kokkos::DualView<double**> y_diag_base{ "y_diag", y_size, y_size };
  Kokkos::DualView<double**> yt_diag_base{ "yt_diag", x_size, x_size };
  Kokkos::DualView<double**> A_subspan_base{ "A_subspan", y_size + 2, x_size + 2 };
  // expressions
  auto x_span()  { return span_adapter(x_span_base, 1, x_size); }
  auto xt_span() { return span_adapter(x_span_base, 1, y_size); }
  auto y_span()  { return span_adapter(y_span_base, 1, y_size); }
  auto yt_span() { return span_adapter(y_span_base, 1, x_size); }
  auto x_diagonal()  { return diag_adapter(x_diag_base); }
  auto xt_diagonal() { return diag_adapter(xt_diag_base); }
  auto y_diagonal()  { return diag_adapter(y_diag_base); }
  auto yt_diagonal() { return diag_adapter(yt_diag_base); }
  auto A_subspan() { return subspan_adapter(A_subspan_base, 1, y_size, 1, x_size); }

  virtual void SetUp(){
    auto A_h = A.view_host();
    A_h(0, 0) = 1.; A_h(0, 1) = 0.; A_h(0 ,2) = 2.;
    A_h(1, 0) = 2.; A_h(1, 1) = 1.; A_h(1, 2) = 3.;
    A_h(2, 0) = 0.; A_h(2, 1) = 0.; A_h(2, 2) = 1.;
    A_h(3, 0) = 2.; A_h(3, 1) = 3.; A_h(3, 2) = 4.;
    A.modify_host();
    set_input(x, { 2., 6., 4. });
    set_input(xt, { 4., 2., 6., 3. });
    // expressions
    set_input(x_span_base, { 1., 2., 3., 4., 5., 6. });
    set_matrix(x_diag_base);
    set_matrix(xt_diag_base);
    set_matrix(y_diag_base);
    set_matrix(yt_diag_base);
    set_matrix(A_subspan_base);
  }

  virtual void TearDown(){}

private:

  template <typename ...ViewProps>
  static void set_input(Kokkos::View<ViewProps...> x, const std::vector<double> &values) {
    assert(x.extent(0) == values.size());
    auto x_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
    for (size_t i = 0; i < x.extent(0); ++i) {
      x_h(i) = values[i];
    }
    Kokkos::deep_copy(x, x_h);
  }

  template <typename ...ViewProps>
  static void set_input(Kokkos::DualView<ViewProps...> x, const std::vector<double> &values) {
    set_input(x.view_device(), values);
    x.modify_device();
    x.sync_host();
  }

  // populates input matrix with unique integer values
  template <typename ...ViewProps>
  static void set_matrix(Kokkos::DualView<ViewProps...> mtx) {
    auto mtx_h = mtx.view_host();
    size_t ex0 = mtx.extent(0), ex1 = mtx.extent(1);
    for (size_t i = 0; i < ex0; ++i) {
      for (size_t j = 0; j < ex1; ++j) {
        mtx_h(i, j) = (double)(i * ex1 + j + 1.0);
      }
    }
    mtx.modify_host();
  }

};

using ops_kokkos = kokkosFixture; // alias for nicer test naming

template <typename TransMode, typename AType, typename XType, typename YType, typename ScalarType>
void test_impl(TransMode trans, ScalarType alpha, AType A, XType x, ScalarType beta, YType y) {
  // copy original values
  Kokkos::View<double*, Kokkos::HostSpace> y_ref("y_ref", y.extent(0));
  y.sync_host();
  auto y_h = y.view_host(); // can't use deep_copy() because y can be [wrapped] Pressio expression
  for (size_t i = 0; i < y.extent(0); ++i) {
    y_ref(i) = y_h(i);
  }

  // call tested routine on device
  A.sync_device();
  x.sync_device();
  y.sync_device();
  // note: explicit instance needed here because we take ref in ::pressio::ops::product()
  auto y_d = y.view_device();
  pressio::ops::product(trans, alpha, A.view_device(), x.view_device(), beta, y_d);
  y.modify_device();

  // call reference gemv() on host
  vanilla_gemv(trans, alpha, A.view_host(), x.view_host(), beta, y_ref);
  y.sync_host();
  for (size_t i = 0; i < y_h.extent(0); ++i) {
    EXPECT_DOUBLE_EQ(y_h(i), y_ref(i));
  }
}

// Important (Thread Safety): modifications on A and y inputs for NaN injection simulation
// breaks thread-safety and assumes sequential execution. It's fine until we introduce threads.
template <typename FixtureType, typename TransMode, typename AType, typename XType, typename YType>
void test_impl(const FixtureType &test, TransMode trans, AType A, XType x, YType y) {
  // alpha = 1, beta = 0, simulate NaN injection in uninitialized y
  y.sync_device();
  ::pressio::ops::fill(y.view_device(), test.NaN);
  y.modify_device();
  test_impl(trans, test.alpha1, test.A, x, test.beta0, y);

  // alpha = 1, beta = 1, reuse values in y
  test_impl(trans, test.alpha1, test.A, x, test.beta1, y);

  // simulate NaN in input
  A.sync_host();
  auto A_h = A.view_host();
  const auto original = A_h(0, 0);
  A_h(0, 0) = test.NaN;
  A.modify_host();

  // alpha = 0, beta = 1, simulate NaN in input
  test_impl(trans, test.alpha0, test.A, x, test.beta1, y);

  // alpha = 0, beta = 0, NaN in both input and result
  ::pressio::ops::fill(y.view_device(), test.NaN);
  y.modify_device();
  test_impl(trans, test.alpha0, test.A, x, test.beta0, y);

  // restore original A
  A_h(0, 0) = original;
  A.modify_host();
}

//-------------------------------------------
// Test plain Kokkos views
//-------------------------------------------

TEST_F(ops_kokkos, dense_mat_vec_vec_NT)
{
  test_impl(*this, pressio::nontranspose(), A, x, y);
}

TEST_F(ops_kokkos, dense_mat_vec_vec_T)
{
  test_impl(*this, pressio::transpose(), A, xt, yt);
}

//-------------------------------------------
// Test plain matrix, x as expression and plain y vector
//-------------------------------------------

TEST_F(ops_kokkos, dense_mat_span_vec_NT)
{
  test_impl(*this, pressio::nontranspose(), A, x_span(), y);
}

TEST_F(ops_kokkos, dense_mat_span_vec_T)
{
  test_impl(*this, pressio::transpose(), A, xt_span(), yt);
}

TEST_F(ops_kokkos, dense_mat_diag_vec_NT)
{
  test_impl(*this, pressio::nontranspose(), A, x_diagonal(), y);
}

TEST_F(ops_kokkos, dense_mat_diag_vec_T)
{
  test_impl(*this, pressio::transpose(), A, xt_diagonal(), yt);
}

//-------------------------------------------
// Test plain matrix, plain x vector and y as expression
//-------------------------------------------

TEST_F(ops_kokkos, dense_mat_vec_span_NT)
{
  test_impl(*this, pressio::nontranspose(), A, x, y_span());
}

TEST_F(ops_kokkos, dense_mat_vec_span_T)
{
  test_impl(*this, pressio::transpose(), A, xt, yt_span());
}

TEST_F(ops_kokkos, dense_mat_vec_diag_NT)
{
  test_impl(*this, pressio::nontranspose(), A, x, y_diagonal());
}

TEST_F(ops_kokkos, dense_mat_vec_diag_T)
{
  test_impl(*this, pressio::transpose(), A, xt, yt_diagonal());
}

//-------------------------------------------
// Test plain matrix with x and y as expressions
//-------------------------------------------

TEST_F(ops_kokkos, dense_mat_span_span_NT)
{
  test_impl(*this, pressio::nontranspose(), A, x_span(), y_span());
}

TEST_F(ops_kokkos, dense_mat_span_span_T)
{
  test_impl(*this, pressio::transpose(), A, xt_span(), yt_span());
}

TEST_F(ops_kokkos, dense_mat_span_diag_NT)
{
  test_impl(*this, pressio::nontranspose(), A, x_span(), y_diagonal());
}

TEST_F(ops_kokkos, dense_mat_span_diag_T)
{
  test_impl(*this, pressio::transpose(), A, xt_span(), yt_diagonal());
}

TEST_F(ops_kokkos, dense_mat_diag_span_NT)
{
  test_impl(*this, pressio::nontranspose(), A, x_diagonal(), y_span());
}

TEST_F(ops_kokkos, dense_mat_diag_span_T)
{
  test_impl(*this, pressio::transpose(), A, xt_diagonal(), yt_span());
}

TEST_F(ops_kokkos, dense_mat_diag_diag_NT)
{
  test_impl(*this, pressio::nontranspose(), A, x_diagonal(), y_diagonal());
}

TEST_F(ops_kokkos, dense_mat_diag_diag_T)
{
  test_impl(*this, pressio::transpose(), A, xt_diagonal(), yt_diagonal());
}

//-------------------------------------------
// Test A matrix as expression with x, y as plain vectors
//-------------------------------------------

TEST_F(ops_kokkos, dense_submat_vec_vec_NT)
{
  test_impl(*this, pressio::nontranspose(), A_subspan(), x, y);
}

TEST_F(ops_kokkos, dense_submat_vec_vec_T)
{
  test_impl(*this, pressio::transpose(), A_subspan(), xt, yt);
}

//-------------------------------------------
// Test A matrix as expression with x as expression and plain y vector
//-------------------------------------------

TEST_F(ops_kokkos, dense_submat_span_vec_NT)
{
  test_impl(*this, pressio::nontranspose(), A_subspan(), x_span(), y);
}

TEST_F(ops_kokkos, dense_submat_span_vec_T)
{
  test_impl(*this, pressio::transpose(), A_subspan(), xt_span(), yt);
}

TEST_F(ops_kokkos, dense_submat_diag_vec_NT)
{
  test_impl(*this, pressio::nontranspose(), A_subspan(), x_diagonal(), y);
}

TEST_F(ops_kokkos, dense_submat_diag_vec_T)
{
  test_impl(*this, pressio::transpose(), A_subspan(), xt_diagonal(), yt);
}

//-------------------------------------------
// Test A matrix as expression with plain x vector and y as expression
//-------------------------------------------

TEST_F(ops_kokkos, dense_submat_vec_span_NT)
{
  test_impl(*this, pressio::nontranspose(), A_subspan(), x, y_span());
}

TEST_F(ops_kokkos, dense_submat_vec_span_T)
{
  test_impl(*this, pressio::transpose(), A_subspan(), xt, yt_span());
}

TEST_F(ops_kokkos, dense_submat_vec_diag_NT)
{
  test_impl(*this, pressio::nontranspose(), A_subspan(), x, y_diagonal());
}

TEST_F(ops_kokkos, dense_submat_vec_diag_T)
{
  test_impl(*this, pressio::transpose(), A_subspan(), xt, yt_diagonal());
}

//-------------------------------------------
// Test A matrix, x and y as expressions
//-------------------------------------------

TEST_F(ops_kokkos, dense_submat_span_span_NT)
{
  test_impl(*this, pressio::nontranspose(), A_subspan(), x_span(), y_span());
}

TEST_F(ops_kokkos, dense_submat_span_span_T)
{
  test_impl(*this, pressio::transpose(), A_subspan(), xt_span(), yt_span());
}

TEST_F(ops_kokkos, dense_submat_span_diag_NT)
{
  test_impl(*this, pressio::nontranspose(), A_subspan(), x_span(), y_diagonal());
}

TEST_F(ops_kokkos, dense_submat_span_diag_T)
{
  test_impl(*this, pressio::transpose(), A_subspan(), xt_span(), yt_diagonal());
}

TEST_F(ops_kokkos, dense_submat_diag_span_NT)
{
  test_impl(*this, pressio::nontranspose(), A_subspan(), x_diagonal(), y_span());
}

TEST_F(ops_kokkos, dense_submat_diag_span_T)
{
  test_impl(*this, pressio::transpose(), A_subspan(), xt_diagonal(), yt_span());
}

TEST_F(ops_kokkos, dense_submat_diag_diag_NT)
{
  test_impl(*this, pressio::nontranspose(), A_subspan(), x_diagonal(), y_diagonal());
}

TEST_F(ops_kokkos, dense_submat_diag_diag_T)
{
  test_impl(*this, pressio::transpose(), A_subspan(), xt_diagonal(), yt_diagonal());
}
