

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc,argv);
  ::testing::InitGoogleMock(&argc,argv);

  MPI_Init(&argc, &argv);  
  RUN_ALL_TESTS();

  MPI_Finalize();
  return 0;
}
