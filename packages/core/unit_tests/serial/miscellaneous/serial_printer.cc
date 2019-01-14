#include <gtest/gtest.h>
#include "CORE_BASIC"

TEST(core_basic, serialPrint){
  auto a1 = "fr";
  int b1 = 2;
  double b2 = 44.4;
  // tellp output position indicator
  //of the current associated streambuf object.

  int n0 = std::cout.tellp();
  // prints with spaces between args
  rompp::core::io::print_stdout(a1,b1,b2);
  int n1 = std::cout.tellp();
  std::cout << n0 << " " << n1 << '\n';
  EXPECT_EQ(n1, n0+9);
}


TEST(core_basic, serialPrintColor){
  using namespace rompp;

  auto a1 = "fr";
  int b1 = 2;
  double b2 = 44.4;

  auto bg = core::io::bg_grey();
  auto col1 = core::io::green();
  auto reset = core::io::reset();
  core::io::print_stdout(bg, col1, a1, b1, b2,
			 reset, "\n");
}
