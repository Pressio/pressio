#include <gtest/gtest.h>
#include "test_helpers.hpp"

#include "pressio/expressions.hpp"

namespace
{

using kv_t = Kokkos::View<double*, Kokkos::HostSpace>;

kv_t create_view(){
  kv_t a("a",5);
  a(0)=1.1;
  a(1)=2.1;
  a(2)=3.1;
  a(3)=4.1;
  a(4)=5.1;
  return a;
}

TEST(expressions_kokkos, span0)
{
  auto a = create_view();

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
  auto a = create_view();

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
  auto a = create_view();
  const auto b = a;
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
    kv_t o("o", 10);
    check_span_traits<double>(o);
  }

  {
    Kokkos::View<double[4], Kokkos::HostSpace> o("o");
    check_span_traits<double>(o);
  }

  {
    kv_t o("o", 10);
    typename kv_t::const_type o2 = o;
    check_span_traits<const double>(o2);
  }
}

} // namespace
