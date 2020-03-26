
#include "pressio_rom.hpp"
#include "pressio_apps.hpp"
#include "../../wls/wls_driver.hpp"
//namespace {

int main(int argc, char *argv[])
{
  using tcomm_t		       = Teuchos::MpiComm<int>;
  using rcpcomm_t	       = Teuchos::RCP<const tcomm_t>;
  using fom_t            = pressio::apps::Burgers1dEigen;
  using scalar_t         = typename fom_t::scalar_type;
  using rom_data_t       = romDataTypeEigen<scalar_t>;
  using ode_tag_euler    = typename ::pressio::ode::implicitmethods::Euler;
  using ode_tag_BDF2     = ::pressio::ode::implicitmethods::BDF2;

  std::string checkStr = doRun< fom_t, rom_data_t, tcomm_t,  pressio::matrixLowerTriangular, ode_tag_euler >();
  std::string checkStr = doRun< fom_t, rom_data_t, tcomm_t,  pressio::matrixLowerTriangular, ode_tag_BDF2 >();

  return 0;
}

//}
