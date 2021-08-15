
#include "pressio_apps.hpp"
#include "utils_tpetra.hpp"
#include "../../helpers/wls_burgers_driver_mpi.hpp"

int main(int argc, char *argv[])
{
  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;
  using fom_t		= pressio::apps::Burgers1dTpetraBlock;
  using scalar_t        = typename fom_t::scalar_type;
  using rom_data_t	= romDataTypeKokkos<scalar_t>;
  using ode_tag_euler   = ::pressio::ode::BDF1;
  using ode_tag_BDF2    = ::pressio::ode::ode::BDF2;
  using lowTri		= pressio::matrixLowerTriangular;
  using phi_native_t = Tpetra::BlockMultiVector<>;

  // scope guard needed for tpetra
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rcpcomm_t Comm = Teuchos::rcp (new tcomm_t(MPI_COMM_WORLD));
    std::string checkStr = "PASSED";

    const auto checkStr1 = pressio::testing::wls::doRun
      <fom_t, rom_data_t, tcomm_t, lowTri, phi_native_t, ode_tag_euler>(Comm, rank);
    const auto checkStr2 = pressio::testing::wls::doRun
      <fom_t, rom_data_t, tcomm_t, lowTri, phi_native_t, ode_tag_BDF2>(Comm, rank);

   if (checkStr1 == "FAILED"){
      std::cout << "WLS failed on implicit Euler" << std::endl;
      checkStr = "FAILED";
    }
    if (checkStr2 == "FAILED"){
      std::cout << "WLS failed on implicit BDF2" << std::endl;
      checkStr = "FAILED";
    }
    std::cout << checkStr << std::endl;
  }
  return 0;
}
