
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

TEST(containers_expressions, spanIsStatic)
{
  using eig_t  = Eigen::VectorXd;
  using my_t   = pressio::containers::Vector<eig_t>;
  using expr_t = pressio::containers::expressions::SpanExpr<my_t>;
  static_assert(expr_t::traits::is_static, "");
}

TEST(containers_expressions, subspanIsStatic)
{
  using eig_t   = Eigen::MatrixXd;
  using my_t    = pressio::containers::DenseMatrix<eig_t>;
  using expr_t= pressio::containers::expressions::SubspanExpr<my_t>;
  static_assert(expr_t::traits::is_static, "");
}

TEST(containers_expressions, diagIsStatic)
{
  using eig_t = Eigen::MatrixXd;
  using my_t = pressio::containers::DenseMatrix<eig_t>;
  using expr_t= pressio::containers::expressions::DiagExpr<my_t>;
  static_assert(expr_t::traits::is_static, "");
}

TEST(containers_expressions, asDiagMatrixIsStatic)
{
  using eig_t = Eigen::VectorXd;
  using my_t  = pressio::containers::Vector<eig_t>;
  using expr_t= pressio::containers::expressions::AsDiagonalMatrixExpr<my_t>;
  static_assert(expr_t::traits::is_static, "");
}
