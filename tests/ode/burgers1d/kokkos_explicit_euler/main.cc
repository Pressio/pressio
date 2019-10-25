
#include "UTILS_ALL"
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "APPS_UNSTEADYBURGERS1D"
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
    using app_state_t		= typename app_t::state_type;
    using app_rhs_t		= typename app_t::velocity_type;

    //-------------------------------
    // create app object
    std::array< scalar_t, 3> mu({5.0, 0.02, 0.02});
    const int Ncell = 20;
    app_t appObj(mu, Ncell);
    appObj.setup();
    auto y0n = appObj.getInitialState();

    // types for ode
    using ode_state_t = pressio::containers::Vector<app_state_t>;
    using ode_res_t   = pressio::containers::Vector<app_rhs_t>;

    ode_state_t y(y0n);
    constexpr auto ode_case = pressio::ode::ExplicitEnum::Euler;
    using stepper_t = pressio::ode::ExplicitStepper
      <ode_case, ode_state_t, app_t, ode_res_t, scalar_t>;
    stepper_t stepperObj(y, appObj);

    // integrate in time
    scalar_t fint = 35;
    scalar_t dt = 0.01;
    auto Nsteps = static_cast<::pressio::ode::types::step_t>(fint/dt);
    pressio::ode::integrateNSteps(stepperObj, y, 0.0, dt, Nsteps);

    using state_type_h = typename app_state_t::HostMirror;
    state_type_h yH("yH", Ncell);
    Kokkos::deep_copy(yH, *y.data());

    {
      using namespace pressio::apps::test;
      checkSol(yH, Burgers1dExpGoldStatesEuler::get(Ncell, dt, fint));
    }
    std::cout << checkStr << std::endl;
  }
  Kokkos::finalize();
  return 0;
}
