
#include <gtest/gtest.h>
#include "pressio_expressions.hpp"

TEST(expressions_eigen, span_traits)
{
  using T = Eigen::VectorXd;
  T o(10);
  using expr_t = decltype(pressio::expressions::span(o, 0, 1));
  static_assert(pressio::traits<expr_t>::is_static, "");
}

TEST(expressions_eigen, subspan_traits)
{
  using T = Eigen::MatrixXd;
  T o(10,10);
  using pair_t = std::pair<int,int>;
  using expr_t = decltype(pressio::expressions::subspan(o, pair_t{0, 1}, pair_t{0,1}));
  static_assert(pressio::traits<expr_t>::is_static, "");
}

TEST(expressions_eigen, diag_traits)
{
  using T = Eigen::MatrixXd;
  T o(10,10);
  using expr_t = decltype(pressio::expressions::diag(o));
  static_assert(pressio::traits<expr_t>::is_static, "");
}

TEST(expressions_eigen, asDiagMatrix_traits)
{
  using T = Eigen::VectorXd;
  T o(10);
  using expr_t = decltype(pressio::expressions::asDiagonalMatrix(o));
  static_assert(pressio::traits<expr_t>::is_static, "");
}
