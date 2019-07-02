
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "APPS_UNSTEADYNONLINADVDIFFREACTION2D"
#include "../gold_states_implicit.hpp"

constexpr double eps = 1e-12;
std::string checkStr {"PASSED"};

template <typename T>
void checkSol(const T & y,
	      const std::vector<double> & trueS){
  if (trueS.empty()) {
    std::cout << " true solution not found, empty " << std::endl;
    checkStr = "FAILED";
  }
  for (size_t i=0; i<trueS.size(); i++){
    if (std::abs(y[i] - trueS[i]) > eps or
	std::isnan(y[i])){
      checkStr = "FAILED";
      break;
    }
  }
}

int main(int argc, char *argv[]){
  using app_t		= pressio::apps::UnsteadyNonLinAdvDiffReac2dEigen;
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;
  using app_jacob_t	= typename app_t::jacobian_type;

  using ode_state_t = pressio::containers::Vector<app_state_t>;
  using ode_res_t   = pressio::containers::Vector<app_rhs_t>;
  using ode_jac_t   = pressio::containers::Matrix<app_jacob_t>;

  constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();

  constexpr int Nx = 11, Ny = Nx*2-1;
  app_t appobj(Nx, Ny);
  appobj.setup();
  const auto y0n = appobj.getInitialState();

  ode_state_t y(y0n);

  // define auxiliary stepper
  using aux_stepper_t = pressio::ode::ImplicitStepper<
    pressio::ode::ImplicitEnum::Euler,
    ode_state_t, ode_res_t, ode_jac_t, app_t>;
  aux_stepper_t stepperAux(y, appobj);

  // the target BDF2 stepper
  constexpr auto ode_case = pressio::ode::ImplicitEnum::BDF2;
  using stepper_t = pressio::ode::ImplicitStepper<
    ode_case, ode_state_t, ode_res_t, ode_jac_t, app_t, aux_stepper_t>;
  stepper_t stepperObj(y, appobj, stepperAux);

  // define solver
  using lin_solver_t = pressio::solvers::iterative::EigenIterative<
    pressio::solvers::linear::iterative::Bicgstab, ode_jac_t>;
  pressio::solvers::NewtonRaphson<scalar_t, lin_solver_t> solverO;
  solverO.setTolerance(1e-14);
  solverO.setMaxIterations(200);

  // integrate in time
  constexpr scalar_t dt = 0.1;
  constexpr auto Nsteps = static_cast<unsigned int>(10);
  constexpr scalar_t fint = Nsteps*dt;
  pressio::ode::integrateNSteps(stepperObj, y, 0.0, dt, Nsteps, solverO);
  std::cout << std::fixed << std::setprecision(14) << *y.data() << std::endl;
  {
    using namespace pressio::apps::test;
    checkSol(y,
  	     NonLinAdvDiffReac2dImpGoldStates<ode_case>::get(Nx, Ny, dt, fint));
  }

  std::cout << checkStr << std::endl;
  return 0;
}
