#include <gtest/gtest.h>
#include "ALGEBRA_BASIC"

TEST(algebra_basic, serialPrint){
  auto a1 = "fr";
  int b1 = 2;
  double b2 = 44.4;
  rompp::utils::io::print_stdout(a1,b1,b2);
}

TEST(algebra_basic, serialPrintColor){
  using namespace rompp;

  auto a1 = "fr";
  int b1 = 2;
  double b2 = 44.4;

  auto bg = utils::io::bg_grey();
  auto col1 = utils::io::green();
  auto reset = utils::io::reset();
  utils::io::print_stdout(bg, col1, a1, b1, b2,
			 reset, "\n");
}
