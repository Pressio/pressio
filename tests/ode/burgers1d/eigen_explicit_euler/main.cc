
#include "pressio_ode.hpp"
#include "pressio_apps.hpp"
#include <array>

constexpr double eps = 1e-12;
std::string checkStr {"PASSED"};

template <typename T>
void checkSol(const T & y, const std::vector<double> & trueS){
  for (size_t i=0; i< trueS.size(); i++){
    if (std::abs(y[i] - trueS[i]) > eps) checkStr = "FAILED";
  }
}

int main(int argc, char *argv[]){
  using app_t		= pressio::apps::Burgers1dEigen;
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;

  //-------------------------------
  // create app object
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  const int Ncell = 20;
  app_t appObj(mu, Ncell);
  auto & y0n = appObj.getInitialState();

  auto j0n = appObj.jacobian(y0n, static_cast<scalar_t>(0));
  std::cout << std::fixed << std::setprecision(15) << y0n << std::endl;
  std::cout << " --- " << std::endl;
  std::cout << std::fixed << std::setprecision(15) << j0n << std::endl;
  std::cout << " --- " << std::endl;

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
  {
    using namespace pressio::apps::test;
    checkSol(y, Burgers1dExpGoldStatesEuler::get(Ncell, dt, fint));
  }
  std::cout << std::fixed << std::setprecision(15) << *y.data() << std::endl;
  std::cout << checkStr << std::endl;
  return 0;
}
