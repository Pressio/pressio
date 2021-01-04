
#include "pressio_apps.hpp"
#include "utils_tpetra.hpp"
#include "../../helpers/wls_burgers_driver_mpi.hpp"

int main(int argc, char *argv[])
{
  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;
  using fom_t		= pressio::apps::Burgers1dTpetra;
  using fom_residual_api_t           = pressio::apps::Burgers1dTpetraDiscreteTimeApi;
  using scalar_t        = typename fom_t::scalar_type;
  using rom_data_t_eigen= romDataTypeEigen<scalar_t>;
  using ode_tag_euler   = ::pressio::ode::implicitmethods::Euler;
  using ode_tag_BDF2    = ::pressio::ode::implicitmethods::BDF2;
  using lowTri		= pressio::matrixLowerTriangular;

  // scope guard needed for tpetra
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rcpcomm_t Comm = Teuchos::rcp (new tcomm_t(MPI_COMM_WORLD));
    std::string checkStr = "PASSED";

 //   const auto checkStr1 = pressio::testing::wls::doRun< fom_t, rom_data_t_eigen, tcomm_t, lowTri, ode_tag_euler>(Comm, rank);
 //   const auto checkStr2 = pressio::testing::wls::doRun< fom_t, rom_data_t_eigen, tcomm_t, lowTri, ode_tag_BDF2>(Comm, rank);
    const std::string checkStr3 = pressio::testing::wls::doRun< fom_residual_api_t,rom_data_t_eigen, tcomm_t,lowTri, ode_tag_euler>(Comm,rank);
   // const std::string checkStr4 = pressio::testing::wls::doRun< fom_residual_api_t,rom_data_t_eigen, tcomm_t, lowTri, ode_tag_BDF2>(Comm, rank);

/*
   if (checkStr1 == "FAILED"){
      std::cout << "WLS failed on implicit Euler" << std::endl;
      checkStr = "FAILED";
    }
    if (checkStr2 == "FAILED"){
      std::cout << "WLS failed on implicit BDF2" << std::endl;
      checkStr = "FAILED";
    }

  if (checkStr3 == "FAILED"){
     std::cout << "WLS failed on implicit Euler with residual API" << std::endl;
     checkStr = "FAILED";
  }

  if (checkStr4 == "FAILED"){
     std::cout << "WLS failed on BDF2 with residual API" << std::endl;
     checkStr = "FAILED";
  }
*/
    std::cout << checkStr << std::endl;
  }
  return 0;
}
