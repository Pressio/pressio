
#include "../../helpers/wls_burgers_driver_serial.hpp"

int main(int argc, char *argv[])
{
  using fom_t              = pressio::apps::Burgers1dKokkos;
  using scalar_t           = typename fom_t::scalar_type;
  using rom_data_t         = romDataTypeKokkos<scalar_t>;
  using ode_tag_euler      = ::pressio::ode::BDF1;
  using ode_tag_BDF2       = ::pressio::ode::ode::BDF2;
  using lowTri		= pressio::matrixLowerTriangular;
  using phi_native_mat_t = typename fom_t::mv_d;

  Kokkos::initialize (argc, argv);
  {
    std::string checkStr = "PASSED";
    const std::string checkStr1 = pressio::testing::wls::doRun
      < fom_t, rom_data_t, lowTri, phi_native_mat_t, ode_tag_euler>();
    const std::string checkStr2 = pressio::testing::wls::doRun
      < fom_t, rom_data_t, lowTri, phi_native_mat_t, ode_tag_BDF2>();

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
  Kokkos::finalize();
  return 0;
}
