
#include <gtest/gtest.h>
#include "pressio_expressions.hpp"

TEST(expressions, spanIsStaticKokkos)
{
  using v_t = Kokkos::View<double*>;
  using span_t= pressio::expressions::SpanExpr<v_t>;
  static_assert(span_t::traits::is_static,"");
}

TEST(expressions, spanIsStaticKokkos2)
{
  using v_t = Kokkos::View<double[4]>;
  using span_t= pressio::expressions::SpanExpr<v_t>;
  static_assert(span_t::traits::is_static,"");
}

TEST(expressions, subspanIsStaticKokkos)
{
  using m_t = Kokkos::View<double**>;
  using subspan_t= pressio::expressions::SubspanExpr<m_t>;
  static_assert(subspan_t::traits::is_static,"");
}

TEST(expressions, subspanIsStaticKokkos2)
{
  using m_t = Kokkos::View<double*[4]>;
  using subspan_t= pressio::expressions::SubspanExpr<m_t>;
  static_assert(subspan_t::traits::is_static,"");
}
