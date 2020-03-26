
#include "pressio_rom.hpp"
#include "pressio_apps.hpp"
#include "utils_tpetra.hpp"
#include "../../wls/wls_driver.hpp"
//namespace {

int main(int argc, char *argv[])
{
  using tcomm_t		       = Teuchos::MpiComm<int>;
  using rcpcomm_t	       = Teuchos::RCP<const tcomm_t>;
  using fom_t     = pressio::apps::Burgers1dTpetraBlock;
  using scalar_t         = typename fom_t::scalar_type;
  using rom_data_t_eigen = romDataTypeEigen<scalar_t>;
  using ode_tag_euler      = typename ::pressio::ode::implicitmethods::Euler;
  using ode_tag_BDF2       = ::pressio::ode::implicitmethods::BDF2;

  // scope guard needed for tpetra
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rcpcomm_t Comm = Teuchos::rcp (new tcomm_t(MPI_COMM_WORLD));
    std::string checkStr = "PASSED";

    std::string checkStr1 = doRun< fom_t, rom_data_t_eigen, tcomm_t,  pressio::matrixLowerTriangular ,ode_tag_euler>(Comm, rank);
    std::string checkStr2 = doRun< fom_t, rom_data_t_eigen, tcomm_t,  pressio::matrixLowerTriangular ,ode_tag_BDF2>(Comm, rank);

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

//}
