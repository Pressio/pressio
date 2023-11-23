#include <gtest/gtest.h>
#include "pressio/expressions.hpp"

namespace
{
template <typename T>
void fillVector(T & a)
{
  a(0) = 1.;
  a(1) = 5.;
  a(2) = 9.;
  a(3) = 13.;
  a(4) = 17.;
  a(5) = 21.;
}
}

TEST(expressions_kokkos, span0)
{
  using kv_t = Kokkos::View<double*, Kokkos::HostSpace>;
  kv_t a("a",5);
  a(0)=1.1;
  a(1)=2.1;
  a(2)=3.1;
  a(3)=4.1;
  a(4)=5.1;

  auto s = pressio::span(a, 2, 2);
  s(0)=10.;
  s(1)=11.;
  EXPECT_EQ(s.extent(0), 2);
  EXPECT_DOUBLE_EQ( a(0), 1.1 );
  EXPECT_DOUBLE_EQ( a(1), 2.1 );
  EXPECT_DOUBLE_EQ( a(2), 10. );
  EXPECT_DOUBLE_EQ( a(3), 11. );
  EXPECT_DOUBLE_EQ( a(4), 5.1 );

  const auto s2 = pressio::span(a, 3, 2);
  s2(0)=20.;
  s2(1)=21.;
  EXPECT_EQ(s2.extent(0), 2);
  EXPECT_DOUBLE_EQ( a(0), 1.1 );
  EXPECT_DOUBLE_EQ( a(1), 2.1 );
  EXPECT_DOUBLE_EQ( a(2), 10. );
  EXPECT_DOUBLE_EQ( a(3), 20. );
  EXPECT_DOUBLE_EQ( a(4), 21. );
}

TEST(expressions_kokkos, span1)
{
  using kv_t = Kokkos::View<double*, Kokkos::HostSpace>;
  kv_t a("a",5);
  a(0)=1.1;
  a(1)=2.1;
  a(2)=3.1;
  a(3)=4.1;
  a(4)=5.1;

  using kv2_t = Kokkos::View<const double*, Kokkos::HostSpace>;
  kv2_t a2 = a;
  auto s = pressio::span(a2, 2, 2);
  static_assert
    (std::is_const<
     typename std::remove_reference<decltype(s(0))>::type
     >::value, "");

  const auto s2 = pressio::span(a2, 2, 2);
  static_assert
    (std::is_const<
     typename std::remove_reference<decltype(s2(0))>::type
     >::value, "");
}

TEST(expressions_kokkos, span2)
{
  using kv_t = Kokkos::View<double*, Kokkos::HostSpace>;
  kv_t a("a",5);
  a(0)=1.1;
  a(1)=2.1;
  a(2)=3.1;
  a(3)=4.1;
  a(4)=5.1;

  const kv_t b = a;
  auto s = pressio::span(b, 2, 2);
  s(0)=10.;
  s(1)=11.;
  EXPECT_EQ(s.extent(0), 2);
  EXPECT_DOUBLE_EQ( a(0), 1.1 );
  EXPECT_DOUBLE_EQ( a(1), 2.1 );
  EXPECT_DOUBLE_EQ( a(2), 10. );
  EXPECT_DOUBLE_EQ( a(3), 11. );
  EXPECT_DOUBLE_EQ( a(4), 5.1 );
  EXPECT_DOUBLE_EQ( b(0), 1.1 );
  EXPECT_DOUBLE_EQ( b(1), 2.1 );
  EXPECT_DOUBLE_EQ( b(2), 10. );
  EXPECT_DOUBLE_EQ( b(3), 11. );
  EXPECT_DOUBLE_EQ( b(4), 5.1 );
}

TEST(expressions_kokkos, span_traits)
{
  {
    using T = Kokkos::View<double*, Kokkos::HostSpace>;
    T o("o", 10);
    using expr_t = decltype(pressio::span(o, 0, 1));
    static_assert(pressio::Traits<expr_t>::rank == 1);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::scalar_type, double>);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::reference_type, double &>);
  }

  {
    using T = Kokkos::View<double[4], Kokkos::HostSpace>;
    T o("o");
    using expr_t = decltype(pressio::span(o, 0, 1));
    static_assert(pressio::Traits<expr_t>::rank == 1);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::scalar_type, double>);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::reference_type, double &>);
  }

  {
    using T = Kokkos::View<double*, Kokkos::HostSpace>;
    T o("o", 10);
    typename T::const_type o2 = o;
    using expr_t = decltype(pressio::span(o2, 0, 1));
    static_assert(pressio::Traits<expr_t>::rank == 1);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::scalar_type, const double>);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::reference_type, const double &>);
  }
}
