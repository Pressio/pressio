
#include <gtest/gtest.h>
#include "pressio_expressions.hpp"

TEST(expressions_eigen, span_traits)
{
  using T = Eigen::VectorXd;
  T o(10);
  using expr_t = decltype(pressio::span(o, 0, 1));
  static_assert(pressio::Traits<expr_t>::is_static, "");

  const T o1(10);
  using expr1_t = decltype(pressio::span(o1, 0, 1));
  static_assert(pressio::Traits<expr1_t>::is_static, "");
}

TEST(expressions_eigen, subspan_traits)
{
  using T = Eigen::MatrixXd;
  T o(10,10);
  using pair_t = std::pair<int,int>;
  using expr_t = decltype(pressio::subspan(o, pair_t{0, 1}, pair_t{0,1}));
  static_assert(pressio::Traits<expr_t>::is_static, "");

  const T o1(10,10);
  using expr1_t = decltype(pressio::subspan(o1, pair_t{0, 1}, pair_t{0,1}));
  static_assert(pressio::Traits<expr1_t>::is_static, "");
}

TEST(expressions_eigen, diag_traits)
{
  using T = Eigen::MatrixXd;
  T o(10,10);
  using expr_t = decltype(pressio::diag(o));
  static_assert(pressio::Traits<expr_t>::is_static, "");

  const T o1(10,10);
  using expr1_t = decltype(pressio::diag(o1));
  static_assert(pressio::Traits<expr1_t>::is_static, "");
}

TEST(expressions_eigen, asDiagMatrix_traits)
{
  using T = Eigen::VectorXd;
  T o(10);
  using expr_t = decltype(pressio::asDiagonalMatrix(o));
  static_assert(pressio::Traits<expr_t>::is_static, "");

  const T o1(10);
  using expr1_t = decltype(pressio::asDiagonalMatrix(o1));
  static_assert(pressio::Traits<expr1_t>::is_static, "");
}
