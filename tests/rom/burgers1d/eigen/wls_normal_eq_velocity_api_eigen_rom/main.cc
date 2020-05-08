
#include "../../wls/wls_burgers_driver_serial.hpp"

int main(int argc, char *argv[])
{
  std::string checkStr = "PASSED";
  using fom_t           = pressio::apps::Burgers1dEigen;
  using fom_residual_api_t           = pressio::apps::Burgers1dEigenResidualApi;
  using scalar_t        = typename fom_t::scalar_type;
  using rom_data_t      = romDataTypeEigen<scalar_t>;
  using ode_tag_euler   = ::pressio::ode::implicitmethods::Euler;
  using ode_tag_BDF2    = ::pressio::ode::implicitmethods::BDF2;
  using lowTri		= pressio::matrixLowerTriangular;

  const std::string checkStr1 = pressio::testing::wls::doRun< fom_t,rom_data_t, lowTri, ode_tag_euler>();
  const std::string checkStr2 = pressio::testing::wls::doRun< fom_t,rom_data_t, lowTri, ode_tag_BDF2>();
  const std::string checkStr3 = pressio::testing::wls::doRun< fom_residual_api_t,rom_data_t, lowTri, ode_tag_euler>();
  const std::string checkStr4 = pressio::testing::wls::doRun< fom_residual_api_t,rom_data_t, lowTri, ode_tag_BDF2>();

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



  std::cout << checkStr << std::endl;
  return 0;
}
