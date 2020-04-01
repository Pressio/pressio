
#include "pressio_rom.hpp"
#include "pressio_apps.hpp"
#include "utils_tpetra.hpp"
#include "../../wls/wls_burgers_driver_serial.hpp"
//namespace {

int main(int argc, char *argv[])
{
  using fom_t              = pressio::apps::Burgers1dKokkos;
  using scalar_t           = typename fom_t::scalar_type;
  using rom_data_t         = romDataTypeKokkos<scalar_t>;
  using ode_tag_euler      = ::pressio::ode::implicitmethods::Euler;
  using ode_tag_BDF2       = ::pressio::ode::implicitmethods::BDF2;

  // scope guard needed for tpetra
  Kokkos::initialize (argc, argv);
  {
    std::string checkStr = "PASSED";

    std::string checkStr1 = doRun< fom_t, rom_data_t, pressio::matrixLowerTriangular, ode_tag_euler>();
    std::string checkStr2 = doRun< fom_t, rom_data_t, pressio::matrixLowerTriangular, ode_tag_BDF2>();

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
