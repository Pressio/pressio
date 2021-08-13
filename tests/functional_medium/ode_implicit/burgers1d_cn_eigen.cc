
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "../testing_apps/apps.hpp"

constexpr double eps = 1e-12;
std::string checkStr {"PASSED"};

std::vector<double> gold=
  {
    5.0209814000128,
    5.044067908724,
    5.0694601439537,
    5.0973757621557,
    5.1280505160954,
    5.1617393080994,
    5.1987172232059,
    5.2392805183926,
    5.2837475208045,
    5.332459322153,
    5.3857799831395,
    5.4440955404205,
    5.5078101629081,
    5.5773358408742,
    5.6530682313685,
    5.7353347085903,
    5.8242904267329,
    5.9197248504917,
    6.0207295624717,
    6.12518281164653};

template <typename T>
void checkSol(const T & y)
{
  for (size_t i=0; i< gold.size(); i++)
    if (std::abs(y(i) - gold[i]) > eps) {
      std::cout << i << std::endl;
      checkStr = "FAILED";
    }
}

int main(int argc, char *argv[]){
  using app_t		= pressio::apps::Burgers1dEigen;
  using scalar_t	= typename app_t::scalar_type;
  using ode_state_t	= typename app_t::state_type;

  //-------------------------------
  // create app object
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  constexpr int Ncell = 20;
  app_t appObj(mu, Ncell);
  auto & y0n = appObj.getInitialState();

  ode_state_t y(y0n);
  auto stepperObj = pressio::ode::create_cranknicolson_stepper(y,appObj);

  using ode_jac_t = typename app_t::jacobian_type;
  using lin_solver_t = pressio::linearsolvers::Solver<
    pressio::linearsolvers::iterative::Bicgstab, ode_jac_t>;
  lin_solver_t linSolverObj;

  auto NonLinSolver=
    pressio::nonlinearsolvers::create_newton_raphson(stepperObj, y, linSolverObj);
  NonLinSolver.setTolerance(1e-11);
  // integrate in time
  scalar_t fint = 35;
  scalar_t dt = 0.01;
  auto Nsteps = static_cast<::pressio::ode::step_count_type>(fint/dt);
  pressio::ode::advance_n_steps(stepperObj, y, 0.0, dt, Nsteps, NonLinSolver);
  std::cout << std::setprecision(14) << y;
  checkSol(y);

  std::cout << checkStr << std::endl;
  return 0;
}
