
#ifdef HAVE_GTEST
#include <gtest/gtest.h>
//#include <gmock/gmock.h>
#include <mpi.h>

struct MPIEnv : public ::testing::Environment{

  int rank_;
  int argc_;
  char ** argv_;

  MPIEnv(int argc, char **argv)
    : argc_(argc), argv_(argv){}

  virtual void SetUp() {
    int mpiError = MPI_Init(&argc_, &argv_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    ASSERT_FALSE(mpiError);
  }

  virtual void TearDown() {
    int mpiError = MPI_Finalize();
    ASSERT_FALSE(mpiError);
  }

  virtual ~MPIEnv() {}
};


int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc,argv);
  //::testing::InitGoogleMock(&argc,argv);
  ::testing::AddGlobalTestEnvironment(new MPIEnv(argc, argv));

  return RUN_ALL_TESTS();
}
#endif
