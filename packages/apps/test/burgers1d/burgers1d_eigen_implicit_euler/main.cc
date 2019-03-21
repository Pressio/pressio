
#include "CORE_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
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
  using app_jacob_t	= typename app_t::jacobian_type;

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
  using ode_jac_t   = rompp::core::Matrix<app_jacob_t>;

  ode_state_t y(y0n);
  using stepper_t = rompp::ode::ImplicitStepper<
    rompp::ode::ImplicitEnum::Euler,
    ode_state_t, ode_res_t, ode_jac_t, app_t>;
  stepper_t stepperObj(appObj, y);

  // define solver
  using lin_solver_t = rompp::solvers::EigenIterative<
    rompp::solvers::linear::Bicgstab, ode_jac_t>;
  rompp::solvers::NewtonRaphson<scalar_t, lin_solver_t> solverO;

  // integrate in time
  scalar_t fint = 35;
  scalar_t dt = 0.01;
  auto Nsteps = static_cast<unsigned int>(fint/dt);
  rompp::ode::integrateNSteps(stepperObj, y, 0.0, dt, Nsteps, solverO);

  checkSol(y, rompp::apps::test::Burg1DtrueImpEulerN20t35);
  std::cout << std::setprecision(14) << *y.data();

  return 0;
}
