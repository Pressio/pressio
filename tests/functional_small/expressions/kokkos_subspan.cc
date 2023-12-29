#include <gtest/gtest.h>
#include "test_helpers.hpp"

#include "pressio/expressions.hpp"

namespace
{

template <typename T>
void fillMatrix(T & A)
{
  A(0,0) = 1.2; A(0,1) = 2.;  A(0,2) = 3.;   A(0,3) = 4.;
  A(1,0) = 5.;  A(1,1) = 6.2; A(1,2) = 7.;   A(1,3) = 8.;
  A(2,0) = 9.;  A(2,1) = 10.; A(2,2) = 11.2; A(2,3) = 12.;
  A(3,0) = 13.; A(3,1) = 14.; A(3,2) = 15.;  A(3,3) = 16.;
  A(4,0) = 12.; A(4,1) = 14.; A(4,2) = 5.;   A(4,3) = 6.;
}

using n_t  = Kokkos::View<double**, Kokkos::HostSpace>;

TEST(expressions_kokkos, subspan0)
{
  const n_t A("A",5,4);
  auto ss = pressio::subspan(A,
			     std::make_pair(0,2),
			     std::make_pair(1,2) );
  static_assert
    (!std::is_const<
     typename std::remove_reference<decltype(ss(0,0))>::type
     >::value, "");
}

TEST(expressions_kokkos, subspan1)
{
  using T  = Kokkos::View<double **, Kokkos::HostSpace>;
  T A("A",5,4);
  fillMatrix(A);

  const auto ss = pressio::subspan(A,
				   std::make_pair(0,2),
				   std::make_pair(2,4) );
  EXPECT_EQ( ss.extent(0), 2 );
  EXPECT_EQ( ss.extent(1), 2 );

  EXPECT_DOUBLE_EQ( ss(0,0), 3.); EXPECT_DOUBLE_EQ( ss(0,1), 4.);
  EXPECT_DOUBLE_EQ( ss(1,0), 7);  EXPECT_DOUBLE_EQ( ss(1,1), 8.);

  //change A, should change expr too
  A(0,3)=43.;
  EXPECT_DOUBLE_EQ( ss(0,1), 43. );
}


TEST(expressions_kokkos, subspan2)
{
  n_t A("A",5,4);
  fillMatrix(A);

  auto ss = pressio::subspan(A,
			     std::make_pair(2,4),
			     std::make_pair(2,4) );
  EXPECT_EQ( ss.extent(0), 2 );
  EXPECT_EQ( ss.extent(1), 2 );
  ss(0,0) = -1.;
  ss(1,1) = -2.;
  EXPECT_DOUBLE_EQ( A(2,2), -1.);
  EXPECT_DOUBLE_EQ( A(3,3), -2.);

  const auto ss2 = pressio::subspan(A,
				    std::make_pair(2,3),
				    std::make_pair(2,4) );
  EXPECT_EQ( ss2.extent(0), 1);
  EXPECT_EQ( ss2.extent(1), 2);
  ss2(0,0) = -3.;
  ss2(0,1) = -4.;
  EXPECT_DOUBLE_EQ( A(2,2), -3.);
  EXPECT_DOUBLE_EQ( A(2,3), -4.);
  EXPECT_DOUBLE_EQ( A(3,3), -2.);
}


TEST(expressions_kokkos, span_traits)
{
  {
    n_t o("o", 10, 10);
    helpers::check_subspan_traits<double>(o);
  }

  {
    Kokkos::View<double[4][4], Kokkos::HostSpace> o("o");
    helpers::check_subspan_traits<double>(o);
  }

  {
    using pair_t = std::pair<std::size_t, std::size_t>;
    n_t o("o",10, 10);
    typename n_t::const_type o2 = o;
    using expr_t = decltype(pressio::subspan(o2, pair_t{0, 1}, pair_t{0,1}));

    static_assert(pressio::Traits<expr_t>::rank == 2);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::scalar_type, const double>);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::reference_type, const double &>);
  }
}

} // namespace
