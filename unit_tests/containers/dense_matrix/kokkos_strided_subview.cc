
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

TEST(containers_sharedmem_kokkos, stridedSubviewLayoutLeft)
{
  using m_t = Kokkos::View<double**, Kokkos::LayoutLeft,  Kokkos::HostSpace>;

  m_t A("A",5,5);
  double val = 0.0;
  for (auto j=0; j<A.extent(1); ++j)
    for (auto i=0; i<A.extent(0); ++i)
      A(i,j) = (val+=1.);

  /*
    1 6  11  16 21
    2 7  12  17 22
    3 8  13  18 23
    4 9  14  19 24
    5 10 15  20 25
   */

  Kokkos::LayoutStride lo(5, A.stride(0)+A.stride(1));
  Kokkos::View<double*,Kokkos::LayoutStride> diag(A.data(), lo);

  for (auto i=0; i<diag.extent(0); ++i)
    std::cout << diag(i) << std::endl;
  ASSERT_EQ(diag(0), 1.);
  ASSERT_EQ(diag(1), 7.);
  ASSERT_EQ(diag(2), 13.);
  ASSERT_EQ(diag(3), 19.);
  ASSERT_EQ(diag(4), 25.);
}


TEST(containers_sharedmem_kokkos, stridedSubviewLayoutRight)
{
  using m_t = Kokkos::View<double**, Kokkos::LayoutRight,  Kokkos::HostSpace>;
  m_t A("A",5,5);
  double val = 0.0;
  for (auto j=0; j<A.extent(1); ++j)
    for (auto i=0; i<A.extent(0); ++i)
      A(i,j) = (val+=1.);

  /*
    1 6  11  16 21
    2 7  12  17 22
    3 8  13  18 23
    4 9  14  19 24
    5 10 15  20 25
   */

  Kokkos::LayoutStride lo(5, A.stride(0)+A.stride(1));
  Kokkos::View<double*,Kokkos::LayoutStride> diag(A.data(), lo);

  for (auto i=0; i<diag.extent(0); ++i)
    std::cout << diag(i) << std::endl;
  ASSERT_EQ(diag(0), 1.);
  ASSERT_EQ(diag(1), 7.);
  ASSERT_EQ(diag(2), 13.);
  ASSERT_EQ(diag(3), 19.);
  ASSERT_EQ(diag(4), 25.);
}
