
#include <gtest/gtest.h>
#include "pressio/ops.hpp"

TEST(ops_kokkos, dense_matrix_clone)
{
  Kokkos::View<double**> A("A", 6,8);
  auto A_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  int c=0;
  for (int i=0; i<6; ++i){
    for (int j=0; j<8; ++j){
     A_h(i,j)= (double) ++c;
   }
  }

  auto B = pressio::ops::clone(A);
  ASSERT_EQ(B.extent(0), 6);
  ASSERT_EQ(B.extent(1), 8);
  ASSERT_FALSE(A.data()==B.data());

  auto B_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), B);
  for (int i=0; i<6; ++i){
    for (int j=0; j<8; ++j){
      ASSERT_DOUBLE_EQ(B_h(i,j), A_h(i,j));
   }
  }
}

TEST(ops_kokkos, dense_matrix_extent)
{
  Kokkos::View<double**> A("A", 6,8);
  ASSERT_TRUE(pressio::ops::extent(A,0) == 6);
  ASSERT_TRUE(pressio::ops::extent(A,1) == 8);
  ASSERT_TRUE(pressio::ops::extent(A,2) == 1); // check extent over the rank
}

TEST(ops_kokkos, dense_matrix_scale)
{
  Kokkos::View<double**> A("A", 6,8);
  KokkosBlas::fill(A, 1.);

  pressio::ops::scale(A, 3.);
  auto A_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  for (int i=0; i<6; ++i){
    for (int j=0; j<8; ++j){
      ASSERT_DOUBLE_EQ(A_h(i,j), 3.);
   }
  }
}


TEST(ops_kokkos, dense_matrix_setzero)
{
  Kokkos::View<double**> A("A", 6,8);
  KokkosBlas::fill(A, 1.);

  pressio::ops::set_zero(A);
  auto A_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  for (int i=0; i<6; ++i){
    for (int j=0; j<8; ++j){
      ASSERT_DOUBLE_EQ(A_h(i,j), 0.);
   }
  }
}

TEST(ops_kokkos, dense_matrix_fill)
{
  Kokkos::View<double**> A("A", 6,8);
  pressio::ops::fill(A, 44.);
  auto A_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  for (int i=0; i<6; ++i){
    for (int j=0; j<8; ++j){
      ASSERT_DOUBLE_EQ(A_h(i,j), 44.);
   }
  }
}

TEST(ops_kokkos, dense_matrix_resize)
{
  Kokkos::View<double**> A("A", 6,8);
  ASSERT_EQ(A.extent(0), 6);
  ASSERT_EQ(A.extent(1), 8);

  pressio::ops::resize(A,3,4);
  ASSERT_EQ(A.extent(0), 3);
  ASSERT_EQ(A.extent(1), 4);
}

TEST(ops_kokkos, dense_matrix_deep_copy)
{
  Kokkos::View<double**> A("A", 6,8);
  pressio::ops::fill(A,44.);

  Kokkos::View<double**> B("B", 6,8);
  pressio::ops::deep_copy(B,A);
  auto B_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), B);
  for (int i=0; i<6; ++i){
   for (int j=0; j<5; ++j){
    ASSERT_DOUBLE_EQ(B_h(i,j),44.);
   }
  }
}

TEST(ops_kokkos, dense_matrix_min_max)
{
  Kokkos::View<double**> A("A", 6, 3);
  auto A_h = Kokkos::create_mirror_view(Kokkos::HostSpace(), A);
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 3; ++j) {
      A_h(i, j) = 100 - i * 6 - j;
    }
  }
  Kokkos::deep_copy(A, A_h);
  ASSERT_DOUBLE_EQ(pressio::ops::min(A), 68.);
  ASSERT_DOUBLE_EQ(pressio::ops::max(A), 100.);
}

TEST(ops_kokkos, dense_matrix_update)
{
  Kokkos::View<double**> M("A", 2, 2);
  Kokkos::View<double**> A("A", 2, 2);
  pressio::ops::fill(M, 1.);
  pressio::ops::fill(A, 2.);

  pressio::ops::update(M, 2., A, 3.);
  auto M_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), M);
  EXPECT_DOUBLE_EQ(M_h(0, 0), 8.);
  EXPECT_DOUBLE_EQ(M_h(0, 1), 8.);
  EXPECT_DOUBLE_EQ(M_h(1, 0), 8.);
  EXPECT_DOUBLE_EQ(M_h(1, 1), 8.);

  // NaN injection through alpha=0
  const auto nan = std::nan("0");
  pressio::ops::fill(M, nan);
  pressio::ops::update(M, 0., A, 2.);
  M_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), M);
  EXPECT_DOUBLE_EQ(M_h(0, 0), 4.);
  EXPECT_DOUBLE_EQ(M_h(0, 1), 4.);
  EXPECT_DOUBLE_EQ(M_h(1, 0), 4.);
  EXPECT_DOUBLE_EQ(M_h(1, 1), 4.);

  // NaN injection through beta=0
  pressio::ops::fill(A, nan);
  pressio::ops::update(M, -1., A, 0.);
  M_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), M);
  EXPECT_DOUBLE_EQ(M_h(0, 0), -4.);
  EXPECT_DOUBLE_EQ(M_h(0, 1), -4.);
  EXPECT_DOUBLE_EQ(M_h(1, 0), -4.);
  EXPECT_DOUBLE_EQ(M_h(1, 1), -4.);

  // alpha=beta=0 corner case
  pressio::ops::fill(M, nan);
  pressio::ops::fill(A, nan);
  pressio::ops::update(M, 0., A, 0.);
  M_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), M);
  EXPECT_DOUBLE_EQ(M_h(0, 0), 0.);
  EXPECT_DOUBLE_EQ(M_h(0, 1), 0.);
  EXPECT_DOUBLE_EQ(M_h(1, 0), 0.);
  EXPECT_DOUBLE_EQ(M_h(1, 1), 0.);
}

TEST(ops_kokkos, dense_matrix_update_epxr)
{
  Kokkos::View<double**> M0("M", 4, 4);
  Kokkos::View<double**> A0("A", 4, 4);
  pressio::ops::fill(M0, 1);
  pressio::ops::fill(A0, 2);
  auto M = pressio::subspan(M0, {1, 3}, {1, 3});
  auto A = pressio::subspan(A0, {1, 3}, {1, 3});

  pressio::ops::update(M, 2., A, 3.);
  auto M0_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), M0);
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      const bool sub = i > 0 && i < 3 && j > 0 && j < 3;
      EXPECT_DOUBLE_EQ(M0_h(i, j), sub ? 8.   // updated M
                                       : 1.); // unmodified part of M0
    }
  }
}