
#include <gtest/gtest.h>

int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc,argv);
  int err = 0;
  {
  	err = RUN_ALL_TESTS();
  }
  return err;
}
