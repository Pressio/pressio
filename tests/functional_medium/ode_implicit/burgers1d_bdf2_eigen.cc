
#include "pressio/ode_implicit.hpp"
#include "../testing_apps/apps.hpp"

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
  using ode_state_t = typename app_t::state_type;

  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  const int Ncell = 20;
  app_t appObj(mu, Ncell);
  auto & y0n = appObj.getInitialState();
  ode_state_t y(y0n);

  auto stepperObj = pressio::ode::create_bdf2_stepper(appObj, y);

  using ode_jac_t = typename app_t::jacobian_type;
  using lin_solver_t = pressio::linearsolvers::Solver<
    pressio::linearsolvers::iterative::Bicgstab, ode_jac_t>;
  lin_solver_t linSolverObj;

  auto NonLinSolver= pressio::nonlinearsolvers::create_newton_raphson(
      stepperObj, y, linSolverObj);

  // integrate in time
  scalar_t fint = 10;
  scalar_t dt = 0.01;
  auto Nsteps = static_cast<::pressio::ode::step_count_type>(fint/dt);
  pressio::ode::advance_n_steps(stepperObj, y, 0.0, dt, Nsteps, NonLinSolver);
  {
    using namespace pressio::apps::test;
    checkSol(y, Burgers1dImpGoldStatesBDF2::get(Ncell, dt, fint));
  }

  std::cout << checkStr << std::endl;
  return 0;
}
