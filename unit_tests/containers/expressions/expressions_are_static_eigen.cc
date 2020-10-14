
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

TEST(containers_expressions, spanIsStatic)
{
  using eig_t = Eigen::VectorXd;
  using my_t = pressio::containers::Vector<eig_t>;
  using span_t= pressio::containers::expressions::SpanExpr<my_t>;
  static_assert
    (pressio::containers::predicates::is_static_vector_wrapper_eigen<span_t>::value, "");
  static_assert
    (!pressio::containers::predicates::is_dynamic_vector_wrapper_eigen<span_t>::value, "");
}

TEST(containers_expressions, subspanIsStatic)
{
  using eig_t = Eigen::MatrixXd;
  using my_t = pressio::containers::DenseMatrix<eig_t>;
  using subspan_t= pressio::containers::expressions::SubspanExpr<my_t>;
  static_assert
    (pressio::containers::predicates::is_static_dense_matrix_wrapper_eigen<subspan_t>::value, "");
  static_assert
    (!pressio::containers::predicates::is_dynamic_dense_matrix_wrapper_eigen<subspan_t>::value, "");
}

TEST(containers_expressions, diagIsStatic)
{
  using eig_t = Eigen::MatrixXd;
  using my_t = pressio::containers::DenseMatrix<eig_t>;
  using diag_t= pressio::containers::expressions::DiagExpr<my_t>;
  static_assert
    (pressio::containers::predicates::is_static_vector_wrapper_eigen<diag_t>::value, "");
  static_assert
    (!pressio::containers::predicates::is_dynamic_vector_wrapper_eigen<diag_t>::value, "");
}
