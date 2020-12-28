
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

TEST(containers_expressions, spanIsStaticKokkos)
{
  using namespace pressio::containers;
  using v_t = Kokkos::View<double*>;
  using my_t = Vector<v_t>;
  using span_t= expressions::SpanExpr<my_t>;
  static_assert(span_t::traits::is_static,"");
}

TEST(containers_expressions, spanIsStaticKokkos2)
{
  using namespace pressio::containers;
  using v_t = Kokkos::View<double[4]>;
  using my_t = Vector<v_t>;
  using span_t= expressions::SpanExpr<my_t>;
  static_assert(span_t::traits::is_static,"");
}

TEST(containers_expressions, subspanIsStaticKokkos)
{
  using namespace pressio::containers;
  using m_t = Kokkos::View<double**>;
  using my_t = DenseMatrix<m_t>;
  using subspan_t= expressions::SubspanExpr<my_t>;
  static_assert(subspan_t::traits::is_static,"");
}

TEST(containers_expressions, subspanIsStaticKokkos2)
{
  using namespace pressio::containers;
  using m_t = Kokkos::View<double*[4]>;
  using my_t = DenseMatrix<m_t>;
  using subspan_t= expressions::SubspanExpr<my_t>;
  static_assert(subspan_t::traits::is_static,"");
}
