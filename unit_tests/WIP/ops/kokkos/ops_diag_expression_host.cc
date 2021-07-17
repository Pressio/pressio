
#include <gtest/gtest.h>
#include "pressio_ops.hpp"

template<typename T>
void fillMatrix(T & A)
{
  A(0,0) = 1.; A(0,1) = 2.; A(0,2) = 3.; A(0,3) = 4.;
  A(1,0) = 1.; A(1,1) = 2.; A(1,2) = 3.; A(1,3) = 4.;
  A(2,0) = 1.; A(2,1) = 2.; A(2,2) = 3.; A(2,3) = 4.;
  A(3,0) = 1.; A(3,1) = 2.; A(3,2) = 3.; A(3,3) = 4.;
}

TEST(ops_kokkos, diag_scale)
{
  using layout = Kokkos::LayoutLeft;
  using kv_t   = Kokkos::View<double**, layout, Kokkos::HostSpace>;
  using pc_t   = pressio::containers::DenseMatrix<kv_t>;

  pc_t A(4,4); fillMatrix(A);
  auto diag = pressio::containers::diag(A);
  pressio::ops::scale(diag, 3.);
  EXPECT_DOUBLE_EQ(diag(0), 3.);
  EXPECT_DOUBLE_EQ(diag(1), 6.);
  EXPECT_DOUBLE_EQ(diag(2), 9.);
  EXPECT_DOUBLE_EQ(diag(3), 12.);

  // diag modifies the matrix itself in this case
  EXPECT_DOUBLE_EQ(A(0,0), 3.);
  EXPECT_DOUBLE_EQ(A(1,1), 6.);
  EXPECT_DOUBLE_EQ(A(2,2), 9.);
  EXPECT_DOUBLE_EQ(A(3,3), 12.);
}

TEST(ops_kokkos, diag_dot_regular_wrapper)
{
  using layout = Kokkos::LayoutLeft;
  using kv_t   = Kokkos::View<double**, layout, Kokkos::HostSpace>;
  using pc_t   = pressio::containers::DenseMatrix<kv_t>;
  pc_t A(4,4); fillMatrix(A);
  auto d1 = pressio::containers::diag(A);

  using vec_t =pressio::containers::Vector<Kokkos::View<double*, Kokkos::HostSpace>>;
  vec_t v(4);
  pressio::ops::fill(v, 1.);

  const auto val = pressio::ops::dot(d1,v);
  EXPECT_DOUBLE_EQ(val, 10.);
}

TEST(ops_kokkos, diag_dot_span)
{
  using layout = Kokkos::LayoutLeft;
  using kv_t   = Kokkos::View<double**, layout, Kokkos::HostSpace>;
  using pc_t   = pressio::containers::DenseMatrix<kv_t>;
  pc_t A(4,4); fillMatrix(A);
  auto d1 = pressio::containers::diag(A);

  using vec_t =pressio::containers::Vector<Kokkos::View<double*, Kokkos::HostSpace>>;
  vec_t v(8);
  pressio::ops::fill(v, 1.);
  // start at 3 and span 4 elements
  auto sp = pressio::containers::span(v, 3, 4);
  pressio::ops::fill(sp, 2.);

  const auto val = pressio::ops::dot(d1, sp);
  EXPECT_DOUBLE_EQ(val, 20.);
}

TEST(ops_kokkos, mat_prod_diag)
{
  using layout = Kokkos::LayoutLeft;
  using kv_t   = Kokkos::View<double**, layout, Kokkos::HostSpace>;
  using pc_t   = pressio::containers::DenseMatrix<kv_t>;
  pc_t A(4,4); fillMatrix(A);
  auto d1 = pressio::containers::diag(A);

  pc_t M(4,4);
  pressio::ops::fill(M, 2.);

  using vec_t =pressio::containers::Vector<Kokkos::View<double*, Kokkos::HostSpace>>;
  vec_t y(4);
  pressio::ops::fill(y, 0.);

  //y = 0*y + 1*M*d1
  pressio::ops::product(pressio::nontranspose(), 1., M, d1, 0., y);
  EXPECT_DOUBLE_EQ(y(0), 20.);
  EXPECT_DOUBLE_EQ(y(1), 20.);
  EXPECT_DOUBLE_EQ(y(2), 20.);
  EXPECT_DOUBLE_EQ(y(3), 20.);
}
