
#include <gtest/gtest.h>
//#include <gmock/gmock.h>
#include <Kokkos_Core.hpp>
#include <memory>

struct KokkosEnv : public ::testing::Environment{
  int argc_;
  char ** argv_;

  KokkosEnv(int argc, char **argv)
    : argc_(argc), argv_(argv){}

  virtual void SetUp() {
    Kokkos::initialize (argc_, argv_);
  }

  virtual void TearDown() {
    Kokkos::finalize();
  }

  virtual ~KokkosEnv() {}
};


int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc,argv);
  //  ::testing::InitGoogleMock(&argc,argv);
  int err = 0;
  {
    std::unique_ptr<KokkosEnv> envPtr(new KokkosEnv(argc, argv));
    ::testing::AddGlobalTestEnvironment(envPtr.get());
    err = RUN_ALL_TESTS();
  }
  return err;
}
