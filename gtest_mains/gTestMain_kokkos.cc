
#include <gtest/gtest.h>
//#include <gmock/gmock.h>
#include <Kokkos_Core.hpp>

class KokkosEnvironment : public ::testing::Environment
{
public:
  int rank_;
  
  virtual void SetUp() {
    char**argv;
    int argc = 0;
    Kokkos::initialize (argc, argv);
  }
  
  virtual void TearDown() {
    Kokkos::finalize();
  }

  virtual ~KokkosEnvironment() {}
};


int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc,argv);
  //  ::testing::InitGoogleMock(&argc,argv);
  ::testing::AddGlobalTestEnvironment(new KokkosEnvironment);
  
  return RUN_ALL_TESTS();
}
