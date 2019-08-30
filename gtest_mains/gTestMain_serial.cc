
#include <gtest/gtest.h>
//#include <gmock/gmock.h>

int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc,argv);
  //::testing::InitGoogleMock(&argc,argv);
  int err = 0;
  {
  	err = RUN_ALL_TESTS();
  }
  return err;
}
