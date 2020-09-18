
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

namespace{

  template <typename T>
  void test1(T & a)
  {
    {
      const auto sp = ::pressio::containers::span(a, 3, 2);
      EXPECT_EQ( sp.extent(0), 2 );
      EXPECT_DOUBLE_EQ( sp[0], 13. );
      EXPECT_DOUBLE_EQ( sp[1], 17. );
    }
    {
      const auto sp = ::pressio::containers::span(a, 2, 1);
      EXPECT_EQ( sp.extent(0), 1 );
      EXPECT_DOUBLE_EQ( sp[0], 9. );
    }
  }

  template <typename T>
  void test2(T & a)
  {
    {
      // change some entries
      auto sp = pressio::containers::span(a, 2, 3);
      EXPECT_EQ( sp.extent(0), 3 );

      // before changing it
      EXPECT_DOUBLE_EQ( sp[0], 9. );
      EXPECT_DOUBLE_EQ( sp[1], 13. );
      EXPECT_DOUBLE_EQ( sp[2], 17. );
      // modify
      sp[0] = 44.;
      // after
      EXPECT_DOUBLE_EQ( sp[0], 44. );
      EXPECT_DOUBLE_EQ( sp[1], 13. );
      EXPECT_DOUBLE_EQ( sp[2], 17. );
    }

    {
      // get the native expression
      auto sp = pressio::containers::span(a, 2, 3);
      auto & natEx = *sp.data();
      EXPECT_EQ( natEx.size(), 3 );
      EXPECT_DOUBLE_EQ( natEx(0), 44. );
      EXPECT_DOUBLE_EQ( natEx(1), 13. );
      EXPECT_DOUBLE_EQ( natEx(2), 17. );
    }
  }

  template <typename T>
  void testConst(const T & a){
    auto sp = pressio::containers::span(a, 2, 3);
    EXPECT_EQ( sp.extent(0), 3 );
    EXPECT_DOUBLE_EQ( sp[0], 44. );
    EXPECT_DOUBLE_EQ( sp[1], 13. );
    EXPECT_DOUBLE_EQ( sp[2], 17. );
    const auto &natEx = *sp.data();
    EXPECT_DOUBLE_EQ( natEx(0), 44. );
    EXPECT_DOUBLE_EQ( natEx(1), 13. );
    EXPECT_DOUBLE_EQ( natEx(2), 17. );
  }
};

TEST(containers_vector_eigen, span)
{
  using eigv_t = Eigen::VectorXd;
  using my_t = pressio::containers::Vector<eigv_t>;

  my_t a(6);
  a(0) = 1.;
  a(1) = 5.;
  a(2) = 9.;
  a(3) = 13.;
  a(4) = 17.;
  a(5) = 21.;

  test1(a);
  test2(a);
  testConst(a);
}

TEST(containers_vector_eigen, spanConstructor){
  using eigv_t = Eigen::VectorXd;
  using my_t = pressio::containers::Vector<eigv_t>;

  my_t a(6);
  a(0) = 1.;
  a(1) = 5.;
  a(2) = 9.;
  a(3) = 13.;
  a(4) = 17.;
  a(5) = 21.;

  // span elements a[2],a[3],a[4]
  const auto sp = pressio::containers::span(a, 2, 3);
  EXPECT_EQ( sp.extent(0), 3 );

  const auto sp2 = pressio::containers::span(a, std::make_pair(2,5));
  EXPECT_EQ( sp2.extent(0), 3 );

  EXPECT_DOUBLE_EQ(sp[0], sp2[0]);
  EXPECT_DOUBLE_EQ(sp[1], sp2[1]);
  EXPECT_DOUBLE_EQ(sp[2], sp2[2]);
}
