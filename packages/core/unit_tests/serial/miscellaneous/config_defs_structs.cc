#include <gtest/gtest.h>
#include "CORE_BASIC"

TEST(core_basic, plusOpStructForExpTemplates){
  using namespace rompp::core;

  exprtemplates::plus_ OP;
  double a = 1.0;
  double b = 1.0;
  auto c = OP(a,b);

  using c_t = decltype(c);
  ::testing::StaticAssertTypeEq<c_t, double>();
  EXPECT_NEAR(c,2.0,1e-12);
}


TEST(core_basic, subtractOpStructForExpTemplates){
  using namespace rompp::core;

  exprtemplates::subtract_ OP;
  double a = 1.1;
  double b = 1.0;
  auto c = OP(a,b);

  using c_t = decltype(c);
  ::testing::StaticAssertTypeEq<c_t, double>();
  EXPECT_NEAR(c,0.1,1e-12);
}


TEST(core_basic, prodOpStructForExpTemplates){
  using namespace rompp::core;

  exprtemplates::times_ OP;
  double a = 1.1;
  double b = 2.0;
  auto c = OP(a,b);

  using c_t = decltype(c);
  ::testing::StaticAssertTypeEq<c_t, double>();
  EXPECT_NEAR(c,2.2,1e-12);
}


