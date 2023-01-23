
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
}