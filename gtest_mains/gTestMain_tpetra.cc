
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <mpi.h>
#include <Tpetra_Core.hpp>

class MPIEnvironment : public ::testing::Environment
{
public:
  
  virtual void SetUp() {
    //    char** argv;
    // int argc = 0;
    // int mpiError = MPI_Init(&argc, &argv);
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    // ASSERT_FALSE(mpiError);

  }
  virtual void TearDown() {
    // int mpiError = MPI_Finalize();
    // ASSERT_FALSE(mpiError);
  }
  virtual ~MPIEnvironment() {}
};


int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc,argv);
  ::testing::InitGoogleMock(&argc,argv);
  //  ::testing::AddGlobalTestEnvironment(new MPIEnvironment);
  int ws_;
  int err;
  
/* from web:
 * https://trilinos.org/docs/dev/packages/tpetra/doc/html/Tpetra_Lesson01.htmlneeded
 * looks like the scope is needed only if one uses mpi via trilinos only
 * if one uses also native mpi, we dont need it
 * the scope guard initializes both MPI and kokkos
 */
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    MPI_Comm_size(MPI_COMM_WORLD, &ws_);
    assert(ws_ > 1);

    err = RUN_ALL_TESTS();
  }
  return err;
}
