
#include <gtest/gtest.h>
#include "pressio/ops.hpp"

using vec_t = Kokkos::View<double*>;
using mat_t = Kokkos::View<double**>;

TEST(ops_kokkos, subspan_clone)
{
  const int m = 6, n = 8;
  Kokkos::View<double**> A("matrix", m + 2, n + 2);
  auto A_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  const auto ex0 = ::pressio::ops::extent(A_h, 0);
  const auto ex1 = ::pressio::ops::extent(A_h, 1);
  for (int i = 0; i < ex0; ++i) {
    for (int j = 0; j < ex1; ++j) {
      A_h(i, j) = (double)(i * ex1 + j + 1); // unique int values from 1
    }
  }
  Kokkos::deep_copy(A, A_h);
  auto ex = pressio::subspan(A, std::make_pair(1, 1 + m), std::make_pair(1, 1 + n));
  auto B = pressio::ops::clone(ex);
  ASSERT_EQ(B.extent(0), m);
  ASSERT_EQ(B.extent(1), n);
  const auto ex_view = ::pressio::ops::impl::get_native(ex);
  const auto B_view = ::pressio::ops::impl::get_native(B);
  ASSERT_FALSE(ex_view.data() == B_view.data());

  const auto ex_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ex_view);
  const auto B_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), B_view);
  for (int i = 0; i < m; ++i){
    for (int j = 0; j < n; ++j){
      ASSERT_DOUBLE_EQ(B_h(i, j), ex_h(i, j));
    }
  }
}

TEST(ops_kokkos, subspan_extent)
{
  mat_t A("A",5,5);
  std::pair<int,int> r(1,3);
  std::pair<int,int> c(2,5);

  auto ex = pressio::subspan(A,r,c);
  ASSERT_TRUE(pressio::ops::extent(ex,0)==2);
  ASSERT_TRUE(pressio::ops::extent(ex,1)==3);
}

TEST(ops_kokkos, subspan_scale)
{
  mat_t A("A",4,5);
  KokkosBlas::fill(A,1.);

  std::pair<int,int> r(1,3);
  std::pair<int,int> c(2,4);
  auto sp = pressio::subspan(A,r,c);
  pressio::ops::scale(sp, 3.);

  auto A_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  ASSERT_DOUBLE_EQ(A_h(0,0),1.);
  ASSERT_DOUBLE_EQ(A_h(0,1),1.);
  ASSERT_DOUBLE_EQ(A_h(0,2),1.);
  ASSERT_DOUBLE_EQ(A_h(0,3),1.);
  ASSERT_DOUBLE_EQ(A_h(0,4),1.);

  ASSERT_DOUBLE_EQ(A_h(1,0),1.);
  ASSERT_DOUBLE_EQ(A_h(1,1),1.);
  ASSERT_DOUBLE_EQ(A_h(1,2),3.);
  ASSERT_DOUBLE_EQ(A_h(1,3),3.);
  ASSERT_DOUBLE_EQ(A_h(1,4),1.);

  ASSERT_DOUBLE_EQ(A_h(2,0),1.);
  ASSERT_DOUBLE_EQ(A_h(2,1),1.);
  ASSERT_DOUBLE_EQ(A_h(2,2),3.);
  ASSERT_DOUBLE_EQ(A_h(2,3),3.);
  ASSERT_DOUBLE_EQ(A_h(2,4),1.);

  ASSERT_DOUBLE_EQ(A_h(3,0),1.);
  ASSERT_DOUBLE_EQ(A_h(3,1),1.);
  ASSERT_DOUBLE_EQ(A_h(3,2),1.);
  ASSERT_DOUBLE_EQ(A_h(3,3),1.);
  ASSERT_DOUBLE_EQ(A_h(3,4),1.);
}

TEST(ops_kokkos, subspan_set_zero)
{
  mat_t A("A",4,5);
  KokkosBlas::fill(A,1.);

  std::pair<int,int> r(1,3);
  std::pair<int,int> c(2,4);
  auto sp = pressio::subspan(A,r,c);
  pressio::ops::set_zero(sp);

  auto A_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  ASSERT_DOUBLE_EQ(A_h(0,0),1.);
  ASSERT_DOUBLE_EQ(A_h(0,1),1.);
  ASSERT_DOUBLE_EQ(A_h(0,2),1.);
  ASSERT_DOUBLE_EQ(A_h(0,3),1.);
  ASSERT_DOUBLE_EQ(A_h(0,4),1.);

  ASSERT_DOUBLE_EQ(A_h(1,0),1.);
  ASSERT_DOUBLE_EQ(A_h(1,1),1.);
  ASSERT_DOUBLE_EQ(A_h(1,2),0.);
  ASSERT_DOUBLE_EQ(A_h(1,3),0.);
  ASSERT_DOUBLE_EQ(A_h(1,4),1.);

  ASSERT_DOUBLE_EQ(A_h(2,0),1.);
  ASSERT_DOUBLE_EQ(A_h(2,1),1.);
  ASSERT_DOUBLE_EQ(A_h(2,2),0.);
  ASSERT_DOUBLE_EQ(A_h(2,3),0.);
  ASSERT_DOUBLE_EQ(A_h(2,4),1.);

  ASSERT_DOUBLE_EQ(A_h(3,0),1.);
  ASSERT_DOUBLE_EQ(A_h(3,1),1.);
  ASSERT_DOUBLE_EQ(A_h(3,2),1.);
  ASSERT_DOUBLE_EQ(A_h(3,3),1.);
  ASSERT_DOUBLE_EQ(A_h(3,4),1.);
}

TEST(ops_kokkos, subspan_fill)
{
  mat_t A("A",4,5);
  KokkosBlas::fill(A,1.);

  std::pair<int,int> r(1,3);
  std::pair<int,int> c(2,4);
  auto sp = pressio::subspan(A,r,c);
  pressio::ops::fill(sp, 44.);

  auto A_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  ASSERT_DOUBLE_EQ(A_h(0,0),1.);
  ASSERT_DOUBLE_EQ(A_h(0,1),1.);
  ASSERT_DOUBLE_EQ(A_h(0,2),1.);
  ASSERT_DOUBLE_EQ(A_h(0,3),1.);
  ASSERT_DOUBLE_EQ(A_h(0,4),1.);

  ASSERT_DOUBLE_EQ(A_h(1,0),1.);
  ASSERT_DOUBLE_EQ(A_h(1,1),1.);
  ASSERT_DOUBLE_EQ(A_h(1,2),44.);
  ASSERT_DOUBLE_EQ(A_h(1,3),44.);
  ASSERT_DOUBLE_EQ(A_h(1,4),1.);

  ASSERT_DOUBLE_EQ(A_h(2,0),1.);
  ASSERT_DOUBLE_EQ(A_h(2,1),1.);
  ASSERT_DOUBLE_EQ(A_h(2,2),44.);
  ASSERT_DOUBLE_EQ(A_h(2,3),44.);
  ASSERT_DOUBLE_EQ(A_h(2,4),1.);

  ASSERT_DOUBLE_EQ(A_h(3,0),1.);
  ASSERT_DOUBLE_EQ(A_h(3,1),1.);
  ASSERT_DOUBLE_EQ(A_h(3,2),1.);
  ASSERT_DOUBLE_EQ(A_h(3,3),1.);
  ASSERT_DOUBLE_EQ(A_h(3,4),1.);
}

TEST(ops_kokkos, subspan_min_max)
{
  mat_t A("A", 6, 6);
  auto A_h = Kokkos::create_mirror_view(Kokkos::HostSpace(), A);
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      A_h(i, j) = 100 - i * 6 - j;
    }
  }
  Kokkos::deep_copy(A, A_h);
  auto A1 = pressio::subspan(A, {1, 3}, {2, 4});

  ASSERT_DOUBLE_EQ(pressio::ops::min(A1), 85.);
  ASSERT_DOUBLE_EQ(pressio::ops::max(A1), 92.);
}
