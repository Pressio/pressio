
#include "CORE_ALL"
#include "ODE_ALL"
#include "APPS_BURGERS1D"
#include "../burgers1d_gold_states.hpp"

constexpr double eps = 1e-12;

template <typename T>
void checkSol(const T & y, const std::vector<double> & trueS){
  for (size_t i=0; i< trueS.size(); i++)
    assert(std::abs(y[i] - trueS[i]) < eps);
}

int main(int argc, char *argv[]){
  using app_t		= rompp::apps::Burgers1dEigen;
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::residual_type;

  //-------------------------------
  // create app object
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  app_t appObj(mu, 20);
  appObj.setup();
  auto & y0n = appObj.getInitialState();
  auto r0n = appObj.residual(y0n, static_cast<scalar_t>(0));

  // types for ode
  using ode_state_t = rompp::core::Vector<app_state_t>;
  using ode_res_t   = rompp::core::Vector<app_rhs_t>;

  ode_state_t y(y0n);
  ode_res_t r(r0n);
  using stepper_t = rompp::ode::ExplicitStepper
    <rompp::ode::ExplicitEnum::Euler, ode_state_t, app_t, ode_res_t>;
  stepper_t stepperObj(appObj, y, r);

  // integrate in time
  scalar_t fint = 35;
  scalar_t dt = 0.01;
  auto Nsteps = static_cast<unsigned int>(fint/dt);
  rompp::ode::integrateNSteps(stepperObj, y, 0.0, dt, Nsteps);
  checkSol(y, rompp::apps::test::Burg1DtrueExpEulerN20t35);
  std::cout << std::setprecision(13) << *y.data();

  return 0;
}
