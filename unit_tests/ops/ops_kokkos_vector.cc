
#include <gtest/gtest.h>
#include "pressio_ops.hpp"

namespace {

template<class T>
void fill1(T v)
{
  using sc_t = typename T::value_type;
  auto v_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
  for (int i=0; i<6; ++i){
   v_h(i) = (sc_t) i;
  }
}

template<class T>
void fill2(T v)
{
  using sc_t = typename T::value_type;
  auto v_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
  for (int i=0; i<6; ++i){
   v_h(i) = -(sc_t) i;
  }
}

template<class T>
void fill3(T v)
{
  auto v_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
  v_h(0) = 1.;
  v_h(1) = 2.;
  v_h(2) = 3.;
}
}//end namespace

TEST(ops_kokkos, vector_clone)
{
  Kokkos::View<double*> a("a", 6);
  fill1(a);
  auto b = pressio::ops::clone(a);
  ASSERT_FALSE( b.data()==a.data());
  ASSERT_EQ(b.extent(0), 6);
  auto a_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
  auto b_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), b);
  for (int i=0; i<6; ++i){
    ASSERT_DOUBLE_EQ(b(i),a_h(i));
  }
}

TEST(ops_kokkos, vector_extent)
{
  Kokkos::View<double*> a("a", 6);
  ASSERT_TRUE(pressio::ops::extent(a,0)== 6);
}

TEST(ops_kokkos, vector_abs)
{
  Kokkos::View<double*> x("a", 6);
  fill2(x);

  Kokkos::View<double*> y("y",6);
  pressio::ops::abs(y,x);

  auto y_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y);
  EXPECT_DOUBLE_EQ(y_h(0), 0.);
  EXPECT_DOUBLE_EQ(y_h(1), 1.);
  EXPECT_DOUBLE_EQ(y_h(2), 2.);
  EXPECT_DOUBLE_EQ(y_h(3), 3.);
  EXPECT_DOUBLE_EQ(y_h(4), 4.);
  EXPECT_DOUBLE_EQ(y_h(5), 5.);
}

TEST(ops_kokkos, vector_setzero)
{
  Kokkos::View<double*> x("a", 6);
  fill1(x);

  pressio::ops::set_zero(x);
  auto x_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
  for (int i=0; i<6; ++i){
    ASSERT_DOUBLE_EQ(x_h(i),0.);
  }
}

TEST(ops_kokkos, vector_scale)
{
  Kokkos::View<double*> x("a", 6);
  fill1(x);

  pressio::ops::scale(x, 2.);
  auto x_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
  EXPECT_DOUBLE_EQ(x_h(0), 0.);
  EXPECT_DOUBLE_EQ(x_h(1), 2.);
  EXPECT_DOUBLE_EQ(x_h(2), 4.);
  EXPECT_DOUBLE_EQ(x_h(3), 6.);
  EXPECT_DOUBLE_EQ(x_h(4), 8.);
  EXPECT_DOUBLE_EQ(x_h(5), 10.);
}

TEST(ops_kokkos, vector_fill)
{
  Kokkos::View<double*> x("a", 6);
  pressio::ops::fill(x, 44.44);
  auto x_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
  EXPECT_DOUBLE_EQ(x_h(0), 44.44);
  EXPECT_DOUBLE_EQ(x_h(1), 44.44);
  EXPECT_DOUBLE_EQ(x_h(2), 44.44);
  EXPECT_DOUBLE_EQ(x_h(3), 44.44);
  EXPECT_DOUBLE_EQ(x_h(4), 44.44);
  EXPECT_DOUBLE_EQ(x_h(5), 44.44);
}

TEST(ops_kokkos, vector_resize)
{
  Kokkos::View<double*> x("a", 6);
  pressio::ops::resize(x,3);
  ASSERT_EQ(x.extent(0), 3);
}

TEST(ops_kokkos, vector_deep_copy)
{
  Kokkos::View<double*> x("a", 6);
  pressio::ops::fill(x, 44.);

  Kokkos::View<double*> b("b", 6);
  pressio::ops::deep_copy(b,x);
  auto b_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), b);
  for (int i=0; i<6; ++i){
    ASSERT_DOUBLE_EQ(b(i),44.);
  }
}

// TEST(ops_kokkos, vector_min_max)
// {
//   using T = Eigen::VectorXd;
//   T a(5);
//   for (int i=0; i<5; ++i){
//    a(i)= (double) i;
//   }

//   ASSERT_DOUBLE_EQ(pressio::ops::min(a), 0.);
//   ASSERT_DOUBLE_EQ(pressio::ops::max(a), 4.);
// }

TEST(ops_kokkos, vector_norm1)
{
  Kokkos::View<double*> x("a", 6);
  fill1(x);
  ASSERT_DOUBLE_EQ(pressio::ops::norm1(x), 15);
}

TEST(ops_kokkos, vector_norm2)
{
  Kokkos::View<double*> x("a", 6);
  fill1(x);
  ASSERT_DOUBLE_EQ(pressio::ops::norm2(x), 7.416198487095663);
}

TEST(ops_kokkos, vector_dot)
{
  Kokkos::View<double*> a("a", 6);
  pressio::ops::fill(a, 1.);
  Kokkos::View<double*> b("b", 6);
  pressio::ops::fill(b, 2.);
  ASSERT_DOUBLE_EQ(pressio::ops::dot(a,b), 12.);
}

TEST(ops_kokkos, vector_pow)
{
  Kokkos::View<double*> a("a", 6);
  fill1(a);

  ::pressio::ops::pow(a, 2.);
  auto a_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
  EXPECT_DOUBLE_EQ(a_h(0), 0.);
  EXPECT_DOUBLE_EQ(a_h(1), 1.);
  EXPECT_DOUBLE_EQ(a_h(2), 4.);
  EXPECT_DOUBLE_EQ(a_h(3), 9.);
  EXPECT_DOUBLE_EQ(a_h(4), 16.);
  EXPECT_DOUBLE_EQ(a_h(5), 25.);
}

TEST(ops_kokkos, vector_absPowPos)
{
  Kokkos::View<double*> x("x", 6);
  fill2(x);

  Kokkos::View<double*> y("y", 6);
  ::pressio::ops::abs_pow(y, x, 3.);

  auto y_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y);
  EXPECT_DOUBLE_EQ(y_h(0), 0.);
  EXPECT_DOUBLE_EQ(y_h(1), 1.);
  EXPECT_DOUBLE_EQ(y_h(2), 8.);
  EXPECT_DOUBLE_EQ(y_h(3), 27.);
  EXPECT_DOUBLE_EQ(y_h(4), 64.);
  EXPECT_DOUBLE_EQ(y_h(5), 125.);
}

TEST(ops_kokkos, vector_absPowNeg)
{
  Kokkos::View<double*> x("x", 6);
  fill2(x);
  Kokkos::View<double*> y("y", 6);
  ::pressio::ops::abs_pow(y, x, -3., 0.00001);

  auto y_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y);
  EXPECT_DOUBLE_EQ(y_h(0), 1./0.00001); // because we guard against diving by zero0.);
  EXPECT_DOUBLE_EQ(y_h(1), 1.);
  EXPECT_DOUBLE_EQ(y_h(2), 1./8.);
  EXPECT_DOUBLE_EQ(y_h(3), 1./27.);
  EXPECT_DOUBLE_EQ(y_h(4), 1./64.);
  EXPECT_DOUBLE_EQ(y_h(5), 1./125.);
}

TEST(ops_kokkos, vector_update1)
{
  Kokkos::View<double*> v("v", 3);
  Kokkos::View<double*> a("a", 3);
  fill3(v);
  fill3(a);

  pressio::ops::update(v, 1., a, 1.);
  auto v_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
  EXPECT_DOUBLE_EQ( v_h(0), 2.0);
  EXPECT_DOUBLE_EQ( v_h(1), 4.0);
  EXPECT_DOUBLE_EQ( v_h(2), 6.0);

  pressio::ops::update(v, a, 1.);
  auto v_h2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
  EXPECT_DOUBLE_EQ( v_h2(0), 1.0);
  EXPECT_DOUBLE_EQ( v_h2(1), 2.0);
  EXPECT_DOUBLE_EQ( v_h2(2), 3.0);
}

TEST(ops_kokkos, vector_update2)
{
  Kokkos::View<double*> v("v", 3);
  Kokkos::View<double*> a("a", 3);
  Kokkos::View<double*> b("b", 3);
  fill3(v);
  fill3(a);
  fill3(b);

  pressio::ops::update(v, 1., a, 1., b, 1.);
  auto v_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
  EXPECT_DOUBLE_EQ( v_h(0), 3.0);
  EXPECT_DOUBLE_EQ( v_h(1), 6.0);
  EXPECT_DOUBLE_EQ( v_h(2), 9.0);

  pressio::ops::update(v, a, 1., b, 1.);
  auto v_h2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
  EXPECT_DOUBLE_EQ( v_h2(0), 2.0);
  EXPECT_DOUBLE_EQ( v_h2(1), 4.0);
  EXPECT_DOUBLE_EQ( v_h2(2), 6.0);
}

TEST(ops_kokkos, vector_update3)
{
  Kokkos::View<double*> v("v", 3);
  Kokkos::View<double*> a("a", 3);
  Kokkos::View<double*> b("b", 3);
  Kokkos::View<double*> c("b", 3);
  fill3(v);
  fill3(a);
  fill3(b);
  fill3(c);

  pressio::ops::update(v, 1., a, 1., b, 1., c, 1.);
  auto v_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
  EXPECT_DOUBLE_EQ( v_h(0), 4.0);
  EXPECT_DOUBLE_EQ( v_h(1), 8.0);
  EXPECT_DOUBLE_EQ( v_h(2), 12.0);

  pressio::ops::update(v, a, 1., b, 1., c, 1.);
  auto v_h2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
  EXPECT_DOUBLE_EQ( v_h2(0), 3.0);
  EXPECT_DOUBLE_EQ( v_h2(1), 6.0);
  EXPECT_DOUBLE_EQ( v_h2(2), 9.0);
}

TEST(ops_kokkos, vector_update4)
{
  Kokkos::View<double*> v("v", 3);
  Kokkos::View<double*> a("a", 3);
  Kokkos::View<double*> b("b", 3);
  Kokkos::View<double*> c("b", 3);
  Kokkos::View<double*> d("b", 3);
  fill3(v);
  fill3(a);
  fill3(b);
  fill3(c);
  fill3(d);

  pressio::ops::update(v, 1., a, 1., b, 1., c, 1., d, 1.);
  auto v_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
  EXPECT_DOUBLE_EQ( v_h(0), 5.0);
  EXPECT_DOUBLE_EQ( v_h(1), 10.0);
  EXPECT_DOUBLE_EQ( v_h(2), 15.0);

  pressio::ops::update(v, a, 1., b, 1., c, 1., d, 1.);
  auto v_h2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
  EXPECT_DOUBLE_EQ( v_h2(0), 4.0);
  EXPECT_DOUBLE_EQ( v_h2(1), 8.0);
  EXPECT_DOUBLE_EQ( v_h2(2), 12.0);
}

TEST(ops_kokkos, vector_elementwiseMultiply)
{
  Kokkos::View<double*> y("y", 3);
  Kokkos::View<double*> x("x", 3);
  Kokkos::View<double*> z("z", 3);
  auto y_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y);
  auto x_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
  auto z_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), z);
  y_h(0) = 1.; y_h(1) = 2.; y_h(2) = 3.;
  x_h(0) = 2.; x_h(1) = 3.; x_h(2) = 4.;
  z_h(0) = 3.; z_h(1) = 4.; z_h(2) = 5.;

  Kokkos::deep_copy(y, y_h);
  Kokkos::deep_copy(x, x_h);
  Kokkos::deep_copy(z, z_h);

  pressio::ops::elementwise_multiply(1., x, z, 1., y);
  Kokkos::deep_copy(y_h, y);
  EXPECT_DOUBLE_EQ( y_h(0), 7.0);
  EXPECT_DOUBLE_EQ( y_h(1), 14.0);
  EXPECT_DOUBLE_EQ( y_h(2), 23.0);
}
