
#include <gtest/gtest.h>
#include <gmock/gmock.h>

int main(int argc, char *argv[])
{
  ::testing::InitGoogleMock(&argc,argv);
  ::testing::InitGoogleTest(&argc,argv);
  int err = 0;
  {
  	err = RUN_ALL_TESTS();
  }
  return err;
}
