
#include "pressio_rom.hpp"
#include "pressio_apps.hpp"
#include "../../wls/wls_driver.hpp"


//namespace {

int main(int argc, char *argv[])
{
  std::string checkStr = "PASSED";
  using fom_t            = pressio::apps::Burgers1dEigen;
  using scalar_t         = typename fom_t::scalar_type;
  using rom_data_t       = romDataTypeEigen<scalar_t>;
  using ode_tag_euler    = typename ::pressio::ode::implicitmethods::Euler;
  using ode_tag_BDF2     = ::pressio::ode::implicitmethods::BDF2;
  std::string checkStr1 = doRun< fom_t,rom_data_t, pressio::matrixLowerTriangular, ode_tag_euler >();
  std::string checkStr2 = doRun< fom_t,rom_data_t, pressio::matrixLowerTriangular, ode_tag_BDF2 >();

  if (checkStr1 == "FAILED"){
     std::cout << "WLS failed on implicit Euler" << std::endl;
     checkStr = "FAILED";
  }
  if (checkStr2 == "FAILED"){
     std::cout << "WLS failed on implicit BDF2" << std::endl;
     checkStr = "FAILED";
  }

  std::cout << checkStr << std::endl;
  return 0;
}

//}
