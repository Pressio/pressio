
#include "pressio_ode.hpp"
#include "pressio_apps.hpp"

constexpr double eps = 1e-12;
std::string checkStr {"PASSED"};

template <typename T>
void checkSol(const T & y, const std::vector<double> & trueS){
  for (size_t i=0; i< trueS.size(); i++)
    if (std::abs(y(i) - trueS[i]) > eps) checkStr = "FAILED";
}

int main(int argc, char *argv[]){
  using app_t		= pressio::apps::Burgers1dEigen;
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;
  using app_jacob_t	= typename app_t::jacobian_type;

  //-------------------------------
  // create app object
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  const int Ncell = 20;
  app_t appObj(mu, Ncell);
  auto & y0n = appObj.getInitialState();

  // types for ode
  using ode_state_t = pressio::containers::Vector<app_state_t>;
  using ode_res_t   = pressio::containers::Vector<app_rhs_t>;
  using ode_jac_t   = pressio::containers::SparseMatrix<app_jacob_t>;

  ode_state_t y(y0n);

  // define auxiliary stepper
  using aux_stepper_t = pressio::ode::ImplicitStepper<
    pressio::ode::implicitmethods::Euler,
    ode_state_t, ode_res_t, ode_jac_t, app_t>;
  aux_stepper_t stepperAux(y, appObj);

  // nonimal stepper
  using ode_tag = pressio::ode::implicitmethods::BDF2;
  using stepper_t = pressio::ode::ImplicitStepper<
    ode_tag, ode_state_t, ode_res_t, ode_jac_t, app_t, aux_stepper_t>;
  stepper_t stepperObj(y, appObj, stepperAux);

  // define solver
  using lin_solver_t = pressio::solvers::linear::Solver<
    pressio::solvers::linear::iterative::Bicgstab, ode_jac_t>;
  lin_solver_t linSolverObj;

  // using nl_solver_t = pressio::solvers::nonlinear::composeNewtonRaphson_t<
  //   stepper_t, lin_solver_t>;
  // nl_solver_t NonLinSolver(stepperObj, y, linSolverObj);  
  auto NonLinSolver= 
    pressio::solvers::nonlinear::createNewtonRaphson(stepperObj, y, linSolverObj);

  // integrate in time
  scalar_t fint = 10;
  scalar_t dt = 0.01;
  auto Nsteps = static_cast<::pressio::ode::types::step_t>(fint/dt);
  pressio::ode::advanceNSteps(stepperObj, y, 0.0, dt, Nsteps, NonLinSolver);
  {
    using namespace pressio::apps::test;
    checkSol(y, Burgers1dImpGoldStatesBDF2::get(Ncell, dt, fint));
  }

  std::cout << checkStr << std::endl;
  return 0;
}
