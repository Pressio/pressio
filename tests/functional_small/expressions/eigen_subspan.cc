
#include <gtest/gtest.h>
#include "test_helpers.hpp"

#include "pressio/expressions.hpp"

namespace
{

template <typename T>
void test1(T & A)
{
  {
    const auto sspan = pressio::subspan(A,
					std::make_pair(0,2),
					std::make_pair(1,2) );
    EXPECT_EQ( sspan.extent(0), 2 ); EXPECT_EQ( sspan.extent(1), 1 );
    EXPECT_DOUBLE_EQ( sspan(0,0), 2. );
    EXPECT_DOUBLE_EQ( sspan(1,0), 6. );
  }
  {
    const auto sspan = pressio::subspan(A,
					std::make_pair(0,3),
					std::make_pair(1,3) );
    EXPECT_EQ( sspan.extent(0), 3 ); EXPECT_EQ( sspan.extent(1), 2 );
    EXPECT_DOUBLE_EQ( sspan(0,0), 2. );  EXPECT_DOUBLE_EQ( sspan(0,1), 3. );
    EXPECT_DOUBLE_EQ( sspan(1,0), 6. );  EXPECT_DOUBLE_EQ( sspan(1,1), 7. );
    EXPECT_DOUBLE_EQ( sspan(2,0), 10. ); EXPECT_DOUBLE_EQ( sspan(2,1), 11. );
  }

  {
    const auto sspan = pressio::subspan(A,
					std::make_pair(2,4),
					std::make_pair(1,3) );
    EXPECT_EQ( sspan.extent(0), 2 ); EXPECT_EQ( sspan.extent(1), 2 );
    EXPECT_DOUBLE_EQ( sspan(0,0), 10. );  EXPECT_DOUBLE_EQ( sspan(0,1), 11. );
    EXPECT_DOUBLE_EQ( sspan(1,0), 14. );  EXPECT_DOUBLE_EQ( sspan(1,1), 15. );
  }
  {
    const auto sspan = pressio::subspan(A,
					std::make_pair(2,3),
					std::make_pair(0,1) );
    EXPECT_EQ( sspan.extent(0), 1 ); EXPECT_EQ( sspan.extent(1), 1 );
    EXPECT_DOUBLE_EQ( sspan(0,0), 9. );
  }
}

template <typename T>
void test2(T & A)
{
  {
    // change some entries
    auto sspan = pressio::subspan(A,
				  std::make_pair(2,4),
				  std::make_pair(1,3) );
    EXPECT_EQ( sspan.extent(0), 2 ); EXPECT_EQ( sspan.extent(1), 2 );

    // before changing it
    EXPECT_DOUBLE_EQ( sspan(0,0), 10. );  EXPECT_DOUBLE_EQ( sspan(0,1), 11. );
    EXPECT_DOUBLE_EQ( sspan(1,0), 14. );  EXPECT_DOUBLE_EQ( sspan(1,1), 15. );
    // modify
    sspan(0,0) = 44.; sspan(0,1) = 33.;
    // after
    EXPECT_DOUBLE_EQ( sspan(0,0), 44. );  EXPECT_DOUBLE_EQ( sspan(0,1), 33. );
    EXPECT_DOUBLE_EQ( sspan(1,0), 14. );  EXPECT_DOUBLE_EQ( sspan(1,1), 15. );
  }

  {
    // get the native expression
    const auto sspan = pressio::subspan(A,
					std::make_pair(2,4),
					std::make_pair(1,3) );
    auto &natEx = sspan.native();
    EXPECT_EQ( natEx.rows(), 2 ); EXPECT_EQ( natEx.cols(), 2 );
    EXPECT_DOUBLE_EQ( natEx(0,0), 44. );  EXPECT_DOUBLE_EQ( natEx(0,1), 33. );
    EXPECT_DOUBLE_EQ( natEx(1,0), 14. );  EXPECT_DOUBLE_EQ( natEx(1,1), 15. );
  }
}

template <typename T>
void testConst(const T & A)
{
  const auto sspan = pressio::subspan(A,
				      std::make_pair(2,4),
				      std::make_pair(1,3) );
  EXPECT_EQ( sspan.extent(0), 2 ); EXPECT_EQ( sspan.extent(1), 2 );
  EXPECT_DOUBLE_EQ( sspan(0,0), 44. );  EXPECT_DOUBLE_EQ( sspan(0,1), 33. );
  EXPECT_DOUBLE_EQ( sspan(1,0), 14. );  EXPECT_DOUBLE_EQ( sspan(1,1), 15. );

  auto & natEx = sspan.native();
  EXPECT_DOUBLE_EQ( natEx(0,0), 44. );  EXPECT_DOUBLE_EQ( natEx(0,1), 33. );
  EXPECT_DOUBLE_EQ( natEx(1,0), 14. );  EXPECT_DOUBLE_EQ( natEx(1,1), 15. );
}

template <typename T>
class EigenSubspanTest : public testing::Test {};

using TestingTypes = ::testing::Types<
    Eigen::MatrixXd,
    Eigen::Matrix<double, -1, -1, Eigen::RowMajor>
>;
TYPED_TEST_SUITE(EigenSubspanTest, TestingTypes);

TYPED_TEST(EigenSubspanTest, baseline)
{
  TypeParam A(6,4);
  helpers::fill_matrix_seq(A);

  test1(A);
  test2(A);
  testConst(A);
}

TEST(expressions_eigen, subspan_traits)
{
  {
    Eigen::MatrixXd o(10,10);
    helpers::check_subspan_traits<double>(o);
  }

  {
    Eigen::Matrix<int,-1,-1> o(10,10);
    helpers::check_subspan_traits<int>(o);
  }

  {
    using pair_t = std::pair<std::size_t, std::size_t>;
    const Eigen::Matrix<int,-1,-1> o(10,10);
    using expr_t = decltype(pressio::subspan(o, pair_t{0, 1}, pair_t{0,1}));

    static_assert(pressio::Traits<expr_t>::rank == 2);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::scalar_type, int>);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::reference_type, int const &>);
  }
}

} // namespace
