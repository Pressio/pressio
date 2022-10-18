
#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"
#include "pressio/rom_subspaces.hpp"

constexpr int m = 10;
constexpr int n = 3;

TEST(rom, linear_col_subspace_constructor_1)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  using space_t = LinearSubspace<basis_t>;
  basis_t A = basis_t::Random(m, n);
  space_t space(A, space_t::SpanningSet::Columns);

  const auto & B = space.basis();
  EXPECT_TRUE(B.data() != A.data());
  EXPECT_TRUE( B.isApprox(A) );
  EXPECT_TRUE( space.dimension() ==n );
  EXPECT_TRUE(  space.isColumnSpace() );
  EXPECT_FALSE( space.isRowSpace() );
}

TEST(rom, linear_col_subspace_constructor_2)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  using space_t = LinearSubspace<basis_t>;
  basis_t A = basis_t::Random(m, n);
  auto Adata = A.data();
  space_t space(std::move(A), space_t::SpanningSet::Columns);

  const auto & B = space.basis();
  EXPECT_TRUE(B.data() == Adata);
  EXPECT_TRUE( space.dimension() ==n );
  EXPECT_TRUE(  space.isColumnSpace() );
  EXPECT_FALSE( space.isRowSpace() );
}

TEST(rom, linear_row_subspace_constructor_1)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  using space_t = LinearSubspace<basis_t>;
  basis_t A = basis_t::Random(m, n);
  space_t space(A, space_t::SpanningSet::Rows);

  const auto & B = space.basis();
  EXPECT_TRUE(B.data() != A.data());
  EXPECT_TRUE( B.isApprox(A) );
  EXPECT_TRUE( space.dimension() ==m );
  EXPECT_FALSE(  space.isColumnSpace() );
  EXPECT_TRUE( space.isRowSpace() );
}

TEST(rom, linear_row_subspace_constructor_2)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  using space_t = LinearSubspace<basis_t>;
  basis_t A = basis_t::Random(m, n);
  auto Adata = A.data();
  space_t space(std::move(A), space_t::SpanningSet::Rows);

  const auto & B = space.basis();
  EXPECT_TRUE(B.data() == Adata);
  EXPECT_TRUE( space.dimension() ==m );
  EXPECT_FALSE(  space.isColumnSpace() );
  EXPECT_TRUE( space.isRowSpace() );
}
