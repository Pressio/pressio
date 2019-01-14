#include <gtest/gtest.h>
#include "CORE_BASIC"

TEST(core_basic, mpiPrint){
  auto a1 = "fr";
  int b1 = 2;
  double b2 = 44.4;
  rompp::core::io::print_stdout(a1,b1,b2);
}
