#include <gtest/gtest.h>
#include "CORE_BASIC"

TEST(core_basic, serialPrint){
  auto a1 = "fr";
  int b1 = 2;
  double b2 = 44.4;
  rompp::core::io::print_stdout(a1,b1,b2);
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
