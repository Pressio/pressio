

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <mpi.h>

class MPIEnvironment : public ::testing::Environment
{
public:
  int rank_;
  
  virtual void SetUp() {
    char** argv;
    int argc = 0;
    int mpiError = MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    ASSERT_FALSE(mpiError);
  }
  virtual void TearDown() {
    int mpiError = MPI_Finalize();
    ASSERT_FALSE(mpiError);
  }
  virtual ~MPIEnvironment() {}
};


int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc,argv);
  ::testing::InitGoogleMock(&argc,argv);
  ::testing::AddGlobalTestEnvironment(new MPIEnvironment);
  
  return RUN_ALL_TESTS();
}
