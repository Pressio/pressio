
#include <gtest/gtest.h>
#include "pressio_ops.hpp"

TEST(ops_vector_kokkos, elementWiseMultiply)
{
  using kv = Kokkos::View<double*, Kokkos::HostSpace>;
  using wt = pressio::containers::Vector<kv>;

  wt a("a",5); pressio::ops::fill(a, 1.);
  wt b("b",5); pressio::ops::fill(b, 2.);
  wt c("c",5); pressio::ops::fill(c, 3.);

  pressio::ops::elementwise_multiply(1.0, a,b, 2., c);

  for (std::size_t i=0; i<5; ++i)
    EXPECT_DOUBLE_EQ(c(i), 8.);
}
