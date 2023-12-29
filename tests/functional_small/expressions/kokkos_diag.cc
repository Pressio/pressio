
#include <gtest/gtest.h>
#include "test_helpers.hpp"

#include "pressio/expressions.hpp"

namespace
{

TEST(expressions_kokkos, diag1)
{
  using n_t = Kokkos::View<double **, Kokkos::HostSpace>;
  n_t A("A", 4,4);
  helpers::fill_matrix(A);

  auto d = pressio::diagonal(A);
  d(0) = 5.;
  EXPECT_EQ( d.extent(0), 4 );
  EXPECT_DOUBLE_EQ( d(0), 5. );
  EXPECT_DOUBLE_EQ( d(1), 6.2 );
  EXPECT_DOUBLE_EQ( d(2), 11.2 );
  EXPECT_DOUBLE_EQ( d(3), 16. );
  EXPECT_DOUBLE_EQ( A(0,0), 5. );
  EXPECT_DOUBLE_EQ( A(1,1), 6.2 );
  EXPECT_DOUBLE_EQ( A(2,2), 11.2 );
  EXPECT_DOUBLE_EQ( A(3,3), 16. );

  const auto d2 = pressio::diagonal(A);
  d2(2) = 55.;
  EXPECT_DOUBLE_EQ( d2(0), 5. );
  EXPECT_DOUBLE_EQ( d2(1), 6.2 );
  EXPECT_DOUBLE_EQ( d2(2), 55. );
  EXPECT_DOUBLE_EQ( d2(3), 16. );
  EXPECT_DOUBLE_EQ( A(0,0), 5. );
  EXPECT_DOUBLE_EQ( A(1,1), 6.2 );
  EXPECT_DOUBLE_EQ( A(2,2), 55. );
  EXPECT_DOUBLE_EQ( A(3,3), 16. );
}

TEST(expressions_kokkos, diag2)
{
  using n_t = Kokkos::View<double**, Kokkos::HostSpace>;
  n_t A("A", 4,4);
  helpers::fill_matrix(A);

  using T = Kokkos::View<const double**, Kokkos::HostSpace>;
  T B = A;

  auto d = pressio::diagonal(B);
  static_assert
    (std::is_const<
     typename std::remove_reference<decltype(d(0))>::type
     >::value, "");
}

TEST(expressions_kokkos, diag_traits)
{
  {
    using T = Kokkos::View<double**, Kokkos::HostSpace>;
    T o("o", 10, 10);
    using expr_t = decltype(pressio::diagonal(o));
    static_assert(pressio::Traits<expr_t>::rank == 1);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::scalar_type, double>);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::reference_type, double &>);
  }

  {
    using T = Kokkos::View<double[4][4], Kokkos::HostSpace>;
    T o("o");
    using expr_t = decltype(pressio::diagonal(o));
    static_assert(pressio::Traits<expr_t>::rank == 1);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::scalar_type, double>);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::reference_type, double &>);
  }

  {
    using T = Kokkos::View<double**, Kokkos::HostSpace>;
    T o("o", 10, 10);
    typename T::const_type o2 = o;
    using expr_t = decltype(pressio::diagonal(o2));
    static_assert(pressio::Traits<expr_t>::rank == 1);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::scalar_type, const double>);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::reference_type, const double &>);
  }
}

} // namespace
