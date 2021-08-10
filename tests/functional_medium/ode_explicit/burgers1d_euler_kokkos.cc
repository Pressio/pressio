
#include "pressio_ode_explicit.hpp"
#include "../testing_apps/apps.hpp"
#include <array>

constexpr double eps = 1e-12;
std::string checkStr {"PASSED"};

template <typename T>
void checkSol(const T & y, const std::vector<double> & trueS){
  for (size_t i=0; i< trueS.size(); i++){
    std::cout << std::setprecision(14) << y(i) << " " << trueS[i] << std::endl;
    if (std::abs(y(i) - trueS[i]) > eps) checkStr = "FAILED";
  }
}

int main(int argc, char *argv[]){
  Kokkos::initialize (argc, argv);
  {
    using app_t			= pressio::apps::Burgers1dKokkos;
    using scalar_t		= typename app_t::scalar_type;
    using state_t		= typename app_t::state_type;

    std::array< scalar_t, 3> mu({5.0, 0.02, 0.02});
    const int Ncell = 20;
    app_t appObj(mu, Ncell);
    auto y0n = appObj.getInitialState();

    state_t y(y0n);
    auto stepperObj = pressio::ode::create_forward_euler_stepper(appObj, y);

    scalar_t fint = 35;
    scalar_t dt = 0.01;
    auto Nsteps = static_cast<::pressio::ode::step_count_type>(fint/dt);
    pressio::ode::advance_n_steps(stepperObj, y, 0.0, dt, Nsteps);

    using state_type_h = typename state_t::HostMirror;
    state_type_h yH("yH", Ncell);
    Kokkos::deep_copy(yH, y);

    {
      using namespace pressio::apps::test;
      checkSol(yH, Burgers1dExpGoldStatesEuler::get(Ncell, dt, fint));
    }
    std::cout << checkStr << std::endl;
  }
  Kokkos::finalize();
  return 0;
}
