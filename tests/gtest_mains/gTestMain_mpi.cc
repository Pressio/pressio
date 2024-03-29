
#include <gtest/gtest.h>
#include <mpi.h>

// struct MPIEnv : public ::testing::Environment{

//   int rank_;
//   int argc_;
//   char ** argv_;

//   MPIEnv(int argc, char **argv)
//     : argc_(argc), argv_(argv){}

//   virtual void SetUp() {
//     int mpiError = MPI_Init(&argc_, &argv_);
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
//     ASSERT_FALSE(mpiError);
//   }

//   virtual void TearDown() {
//     int mpiError = MPI_Finalize();
//     ASSERT_FALSE(mpiError);
//   }

//   virtual ~MPIEnv() {}
// };


int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc,argv);
  int err = 0;
  {
    // std::unique_ptr<MPIEnv> envPtr(new MPIEnv(argc, argv));
    // ::testing::AddGlobalTestEnvironment(envPtr.get());

    int rank = {};
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // ASSERT_FALSE(mpiError);

    err = RUN_ALL_TESTS();

    MPI_Finalize();
    //ASSERT_FALSE(mpiError);
  }
  return err;
}
