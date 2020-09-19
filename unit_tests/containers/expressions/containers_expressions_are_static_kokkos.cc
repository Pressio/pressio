
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

TEST(containers_expressions, spanIsStaticKokkos)
{
  using namespace pressio::containers;
  using v_t = Kokkos::View<double*>;
  using my_t = Vector<v_t>;
  using span_t= expressions::SpanExpr<my_t>;
  static_assert
    (predicates::is_static_vector_wrapper_kokkos<span_t>::value, "");
  static_assert
    (!predicates::is_dynamic_vector_wrapper_kokkos<span_t>::value, "");
}

TEST(containers_expressions, spanIsStaticKokkos2)
{
  using namespace pressio::containers;
  using v_t = Kokkos::View<double[4]>;
  using my_t = Vector<v_t>;
  using span_t= expressions::SpanExpr<my_t>;
  static_assert
    (predicates::is_static_vector_wrapper_kokkos<span_t>::value, "");
  static_assert
    (!predicates::is_dynamic_vector_wrapper_kokkos<span_t>::value, "");
}

TEST(containers_expressions, subspanIsStaticKokkos)
{
  using namespace pressio::containers;
  using m_t = Kokkos::View<double**>;
  using my_t = DenseMatrix<m_t>;
  using subspan_t= expressions::SubspanExpr<my_t>;
  static_assert
    (predicates::is_static_dense_matrix_wrapper_kokkos<subspan_t>::value, "");
  static_assert
    (!predicates::is_dynamic_dense_matrix_wrapper_kokkos<subspan_t>::value, "");
}

TEST(containers_expressions, subspanIsStaticKokkos2)
{
  using namespace pressio::containers;
  using m_t = Kokkos::View<double*[4]>;
  using my_t = DenseMatrix<m_t>;
  using subspan_t= expressions::SubspanExpr<my_t>;
  static_assert
    (predicates::is_static_dense_matrix_wrapper_kokkos<subspan_t>::value, "");
  static_assert
    (!predicates::is_dynamic_dense_matrix_wrapper_kokkos<subspan_t>::value, "");
}
