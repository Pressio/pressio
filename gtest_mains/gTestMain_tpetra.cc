
#ifdef HAVE_GTEST
#include <gtest/gtest.h>
//#include <gmock/gmock.h>
#include <mpi.h>
#include <Tpetra_Core.hpp>

int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc,argv);
  //  ::testing::InitGoogleMock(&argc,argv);
  int ws_;
  int err;

/* from web:
 * https://trilinos.org/docs/dev/packages/tpetra/doc/html/Tpetra_Lesson01.htmlneeded
 * the tpetra scope is needed if MPI is initialized within trilinos
 * the scope guard initializes both MPI and kokkos
 * if one explicity calls MPI_Init we dont need tpetraScope
 */
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    MPI_Comm_size(MPI_COMM_WORLD, &ws_);
    assert(ws_ > 1);
    err = RUN_ALL_TESTS();
  }
  return err;
}
#endif
