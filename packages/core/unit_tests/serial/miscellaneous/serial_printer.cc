#include <gtest/gtest.h>
#include "CORE_BASIC"

TEST(core_basic, serialPrint){
  auto a1 = "fr";
  int b1 = 2;
  double b2 = 44.4;
  // tellp output position indicator of the current associated streambuf object.

  int n0 = std::cout.tellp();
  // prints with spaces between args
  rompp::core::debug::print(a1,b1,b2);
  int n1 = std::cout.tellp();
  std::cout << n0 << " " << n1 << '\n';
  EXPECT_EQ(n1, n0+9);

}
