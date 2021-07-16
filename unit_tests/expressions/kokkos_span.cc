#include <gtest/gtest.h>
#include "pressio_expressions.hpp"

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
};

TEST(containers_expressions_kokkos, span0)
{
  using kv_t        = Kokkos::View<double *, Kokkos::HostSpace>;
  using pressio_v_t = pressio::containers::Vector<kv_t>;
  kv_t a("a",5);
  const pressio_v_t aw(a);
  auto s = pressio::containers::span(aw, 3, 2);
  static_assert
    (std::is_const<
     typename std::remove_reference<decltype(s(0))>::type
     >::value, "");
}

TEST(containers_expressions_kokkos, span1)
{
  using n_t = Kokkos::View<double *, Kokkos::HostSpace>;
  using my_t = pressio::containers::Vector<n_t>;

  my_t a(6);
  fillVector(a);

  const auto sp = pressio::containers::span(a, 3, 2);
  EXPECT_EQ( sp.extent(0), 2 );
  EXPECT_DOUBLE_EQ( sp(0), 13. );
  EXPECT_DOUBLE_EQ( sp(1), 17. );
  EXPECT_DOUBLE_EQ( a(3), 13. );
  EXPECT_DOUBLE_EQ( a(4), 17. );

  //change a, should change span too
  a(3)=43.;
  EXPECT_DOUBLE_EQ( sp(0), 43. );
}

TEST(containers_expressions_kokkos, span2)
{
  using n_t = Kokkos::View<double *, Kokkos::HostSpace>;
  using my_t = pressio::containers::Vector<n_t>;
  my_t a(6);
  fillVector(a);
  EXPECT_DOUBLE_EQ( a(3), 13. );
  EXPECT_DOUBLE_EQ( a(4), 17. );

  auto sp = pressio::containers::span(a, 3, 2);
  // sp can be assigned to
  static_assert
    (!std::is_const<
     typename std::remove_reference<decltype(sp(0))>::type
     >::value, "");

  // change values
  sp(0) = 44.;
  sp(1) = 55.;

  // both sp and a should be changed
  EXPECT_DOUBLE_EQ( sp(0), 44. );
  EXPECT_DOUBLE_EQ( sp(1), 55. );
  EXPECT_DOUBLE_EQ( a(3), 44. );
  EXPECT_DOUBLE_EQ( a(4), 55. );
}

TEST(containers_expressions_kokkos, span3)
{
  using n_t = Kokkos::View<double *, Kokkos::HostSpace>;
  using my_t = pressio::containers::Vector<n_t>;
  n_t a0("a0", 6);
  fillVector(a0);
  const my_t a(a0);
  static_assert
    (std::is_const<
     typename std::remove_reference<decltype(a(0))>::type
     >::value, "");

  EXPECT_DOUBLE_EQ( a(3), 13. );
  EXPECT_DOUBLE_EQ( a(4), 17. );

  auto sp = pressio::containers::span(a, 3, 2);
  // since a is const, sp should be read-only
  static_assert
    (std::is_const<
     typename std::remove_reference<decltype(sp(0))>::type
     >::value, "");

  EXPECT_DOUBLE_EQ( sp(0), 13. );
  EXPECT_DOUBLE_EQ( sp(1), 17. );
  EXPECT_DOUBLE_EQ( a(3),  13. );
  EXPECT_DOUBLE_EQ( a(4),  17. );
}
