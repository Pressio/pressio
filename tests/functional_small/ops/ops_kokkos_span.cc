
#include <gtest/gtest.h>
#include "pressio/ops.hpp"

using vec_t = Kokkos::View<double*>;
using mat_t = Kokkos::View<double**>;

TEST(ops_kokkos, span_extent)
{
  vec_t a("a", 8);
  auto ex = pressio::span(a,5,2);
  ASSERT_TRUE(pressio::ops::extent(ex,0)==2);
}

TEST(ops_kokkos, span_abs)
{
  vec_t a("a", 8);
  KokkosBlas::fill(a, -1.);
  auto ex = pressio::span(a,5,2);

  vec_t y("y", 2);
  pressio::ops::abs(y,ex);
  auto y_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y);
  ASSERT_DOUBLE_EQ(y_h(0),1.);
  ASSERT_DOUBLE_EQ(y_h(1),1.);
}

TEST(ops_kokkos, span_scale)
{
  vec_t a("a",6);
  KokkosBlas::fill(a, 1.);

  auto sp = pressio::span(a, 1,3);
  pressio::ops::scale(sp, 3.);
  auto a_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
  ASSERT_DOUBLE_EQ(a_h(0),1.);
  ASSERT_DOUBLE_EQ(a_h(1),3.);
  ASSERT_DOUBLE_EQ(a_h(2),3.);
  ASSERT_DOUBLE_EQ(a_h(3),3.);
  ASSERT_DOUBLE_EQ(a_h(4),1.);
  ASSERT_DOUBLE_EQ(a_h(5),1.);
}

TEST(ops_kokkos, span_set_zero)
{
  vec_t a("a",6);
  KokkosBlas::fill(a, 1.2);

  auto sp = pressio::span(a, 1,3);
  pressio::ops::set_zero(sp);
  auto a_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
  ASSERT_DOUBLE_EQ(a_h(0),1.2);
  ASSERT_DOUBLE_EQ(a_h(1),0.);
  ASSERT_DOUBLE_EQ(a_h(2),0.);
  ASSERT_DOUBLE_EQ(a_h(3),0.);
  ASSERT_DOUBLE_EQ(a_h(4),1.2);
  ASSERT_DOUBLE_EQ(a_h(5),1.2);
}

TEST(ops_kokkos, span_fill)
{
  vec_t a("a",6);
  KokkosBlas::fill(a, 1.2);

  auto sp = pressio::span(a, 1,3);
  pressio::ops::fill(sp, 44.);
  auto a_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
  ASSERT_DOUBLE_EQ(a_h(0),1.2);
  ASSERT_DOUBLE_EQ(a_h(1),44.);
  ASSERT_DOUBLE_EQ(a_h(2),44.);
  ASSERT_DOUBLE_EQ(a_h(3),44.);
  ASSERT_DOUBLE_EQ(a_h(4),1.2);
  ASSERT_DOUBLE_EQ(a_h(5),1.2);
}

TEST(ops_kokkos, span_norms)
{
  vec_t a("a",6);
  auto a_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
  for (int i=0; i<6; ++i){
    a_h(i) = (double) i;
  }
  Kokkos::deep_copy(a, a_h);

  auto sp = pressio::span(a, 1,3);
  ASSERT_DOUBLE_EQ(pressio::ops::norm1(sp), 6.);
  ASSERT_DOUBLE_EQ(pressio::ops::norm2(sp), std::sqrt(14.));
}

TEST(ops_kokkos, span_dot)
{
  vec_t a("a",6);
  auto a_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
  for (int i=0; i<6; ++i){
    a_h(i) = (double) i;
  }
  Kokkos::deep_copy(a, a_h);
  auto sp = pressio::span(a, 1,3);

  vec_t b("b",3);
  KokkosBlas::fill(b,2.);
  ASSERT_DOUBLE_EQ(pressio::ops::dot(sp, b), 12.);
}

TEST(ops_kokkos, span_pow)
{
  vec_t a("a",8);
  auto a_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
  for (int i=2; i<8; ++i){
    a_h(i) = (double) (i-2);
  }
  Kokkos::deep_copy(a, a_h);

  auto sp = pressio::span(a, 2,6);
  pressio::ops::pow(sp, 2.);
  Kokkos::deep_copy(a, a_h);

  EXPECT_DOUBLE_EQ(a_h(0), 0);
  EXPECT_DOUBLE_EQ(a_h(1), 0);
  EXPECT_DOUBLE_EQ(a_h(2), 0);
  EXPECT_DOUBLE_EQ(a_h(3), 1);
  EXPECT_DOUBLE_EQ(a_h(4), 4);
  EXPECT_DOUBLE_EQ(a_h(5), 9);
  EXPECT_DOUBLE_EQ(a_h(6), 16);
  EXPECT_DOUBLE_EQ(a_h(7), 25);
}

namespace {
vec_t createVectorForUpdate()
{
  vec_t a("a", 6);
  auto a_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
  a_h(2) = 1.;
  a_h(3) = 2.;
  a_h(4) = 3.;
  Kokkos::deep_copy(a, a_h);
  return a;
}
}

TEST(ops_kokkos, span_update1)
{
  auto M1 = createVectorForUpdate();
  auto d1 = pressio::span(M1,2,3);
  auto M2 = createVectorForUpdate();
  auto d2 = pressio::span(M2,2,3);

  pressio::ops::update(d1, 1., d2, 1.);
  auto M1_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), M1);
  EXPECT_DOUBLE_EQ( M1_h(0), 0.0);
  EXPECT_DOUBLE_EQ( M1_h(1), 0.0);
  EXPECT_DOUBLE_EQ( M1_h(2), 2.0);
  EXPECT_DOUBLE_EQ( M1_h(3), 4.0);
  EXPECT_DOUBLE_EQ( M1_h(4), 6.0);
  EXPECT_DOUBLE_EQ( M1_h(5), 0.0);

  pressio::ops::update(d1, 0., d2, 1.);
  Kokkos::deep_copy(M1_h, M1);
  EXPECT_DOUBLE_EQ( M1_h(0), 0.0);
  EXPECT_DOUBLE_EQ( M1_h(1), 0.0);
  EXPECT_DOUBLE_EQ( M1_h(2), 1.0);
  EXPECT_DOUBLE_EQ( M1_h(3), 2.0);
  EXPECT_DOUBLE_EQ( M1_h(4), 3.0);
  EXPECT_DOUBLE_EQ( M1_h(5), 0.0);
}

TEST(ops_kokkos, span_update2)
{
  auto M1 = createVectorForUpdate();
  auto v = pressio::span(M1,2,3);
  auto M2 = createVectorForUpdate();
  auto a = pressio::span(M2,2,3);
  auto M3 = createVectorForUpdate();
  auto b = pressio::span(M3,2,3);

  pressio::ops::update(v, 1., a, 1., b, 1.);
  auto M1_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), M1);
  EXPECT_DOUBLE_EQ( M1_h(2), 3.0);
  EXPECT_DOUBLE_EQ( M1_h(3), 6.0);
  EXPECT_DOUBLE_EQ( M1_h(4), 9.0);

  pressio::ops::update(v, 0., a, 1., b, 1.);
  Kokkos::deep_copy(M1_h, M1);
  EXPECT_DOUBLE_EQ( M1_h(2), 2.0);
  EXPECT_DOUBLE_EQ( M1_h(3), 4.0);
  EXPECT_DOUBLE_EQ( M1_h(4), 6.0);
}

TEST(ops_kokkos, span_update3)
{
  auto M1 = createVectorForUpdate();
  auto v = pressio::span(M1,2,3);
  auto M2 = createVectorForUpdate();
  auto a = pressio::span(M2,2,3);
  auto M3 = createVectorForUpdate();
  auto b = pressio::span(M3,2,3);
  auto M4 = createVectorForUpdate();
  auto c = pressio::span(M4,2,3);

  pressio::ops::update(v, 1., a, 1., b, 1., c, 1.);
  auto M1_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), M1);
  EXPECT_DOUBLE_EQ( M1_h(2), 4.0);
  EXPECT_DOUBLE_EQ( M1_h(3), 8.0);
  EXPECT_DOUBLE_EQ( M1_h(4), 12.0);

  pressio::ops::update(v, 0., a, 1., b, 1., c, 1.);
  Kokkos::deep_copy(M1_h, M1);
  EXPECT_DOUBLE_EQ( M1_h(2), 3.0);
  EXPECT_DOUBLE_EQ( M1_h(3), 6.0);
  EXPECT_DOUBLE_EQ( M1_h(4), 9.0);
}

TEST(ops_kokkos, span_update4)
{
  auto M1 = createVectorForUpdate();
  auto v = pressio::span(M1,2,3);
  auto M2 = createVectorForUpdate();
  auto a = pressio::span(M2,2,3);
  auto M3 = createVectorForUpdate();
  auto b = pressio::span(M3,2,3);
  auto M4 = createVectorForUpdate();
  auto c = pressio::span(M4,2,3);
  auto M5 = createVectorForUpdate();
  auto d = pressio::span(M5,2,3);

  pressio::ops::update(v, 1., a, 1., b, 1., c, 1., d, 1.);
  auto M1_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), M1);
  EXPECT_DOUBLE_EQ( M1_h(2), 5.0);
  EXPECT_DOUBLE_EQ( M1_h(3), 10.0);
  EXPECT_DOUBLE_EQ( M1_h(4), 15.0);

  pressio::ops::update(v, 0., a, 1., b, 1., c, 1., d, 1.);
  Kokkos::deep_copy(M1_h, M1);
  EXPECT_DOUBLE_EQ( M1_h(2), 4.0);
  EXPECT_DOUBLE_EQ( M1_h(3), 8.0);
  EXPECT_DOUBLE_EQ( M1_h(4), 12.0);
}

TEST(ops_kokkos, span_elementwiseMultiply)
{
  vec_t M1("M1",6);
  vec_t M2("M2",7);
  vec_t M3("M3",8);
  auto M1_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), M1);
  auto M2_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), M2);
  auto M3_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), M3);

  M1_h(3)=1.; M1_h(4)=2.; M1_h(5)=3.;
  M2_h(2)=2.; M2_h(3)=3.; M2_h(4)=4.;
  M3_h(4)=3.; M3_h(5)=4.; M3_h(6)=5.;
  Kokkos::deep_copy(M1, M1_h);
  Kokkos::deep_copy(M2, M2_h);
  Kokkos::deep_copy(M3, M3_h);

  auto y = pressio::span(M1,3,3);
  const auto x = pressio::span(M2,2,3);
  const auto z = pressio::span(M3,4,3);

  pressio::ops::elementwise_multiply(1., x, z, 1., y);
  Kokkos::deep_copy(M1_h, M1);
  EXPECT_DOUBLE_EQ( M1_h(3), 7.0);
  EXPECT_DOUBLE_EQ( M1_h(4), 14.0);
  EXPECT_DOUBLE_EQ( M1_h(5), 23.0);
}
