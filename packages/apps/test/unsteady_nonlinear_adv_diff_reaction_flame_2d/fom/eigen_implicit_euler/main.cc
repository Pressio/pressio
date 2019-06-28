
#include "ALGEBRA_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "APPS_UNSTEADYNONLINADVDIFFREACTIONFLAME2D"
#include "../gold_states_implicit.hpp"

constexpr double eps = 1e-10;
std::string checkStr {"PASSED"};

template <typename T>
void checkSol(const T & y,
	      const std::vector<double> & trueS){
  if (trueS.empty()) {
    std::cout << " true solution not found, empty " << std::endl;
    checkStr = "FAILED";
  }
  for (size_t i=0; i<trueS.size(); i++){
    const auto err = std::abs(y[i] - trueS[i]);
    std::cout << std::fixed << std::setprecision(15)
	      << " true = " << trueS[i]
	      << " y = " << y[i]
	      << " err = " << err
	      << std::endl;
    if ( err > eps or std::isnan(y[i])) checkStr = "FAILED";
  }
}

constexpr bool do_print = false;

struct Observer{
  Observer() = default;

  template <typename T>
  void operator()(size_t step, double t, const T & y)
  {
    if (do_print){
      if (step % 5 == 0){
    	std::ofstream file;
    	file.open( "sol_" + std::to_string(step) + ".txt" );
    	for(auto i=0; i < y.size(); i++){
    	  file << std::fixed << std::setprecision(14) << y[i] ;
    	  file << std::endl;
    	}
    	file.close();
      }
    }
  }
};

int main(int argc, char *argv[]){
  using app_t		= rompp::apps::UnsteadyNonLinAdvDiffReacFlame2dEigen;
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_residual_t	= typename app_t::residual_type;
  using app_jacob_t	= typename app_t::jacobian_type;
  constexpr auto zero = ::rompp::algebra::constants::zero<scalar_t>();

  constexpr int Nx = 12, Ny = 6;
  app_t appobj(Nx, Ny);
  appobj.setup();
  const auto y0n = appobj.getInitialState();
  const auto r0n = appobj.residual(y0n, zero);

  if (do_print){
    auto X = appobj.getX(); auto Y = appobj.getY();
    std::ofstream file; file.open( "xy.txt" );
    for(auto i=0; i < X.size(); i++){
      file << std::setprecision(14)
	   << X[i] << " " << Y[i];
      file << std::endl;
    }
    file.close();
  }

  using ode_state_t = rompp::algebra::Vector<app_state_t>;
  using ode_res_t   = rompp::algebra::Vector<app_residual_t>;
  using ode_jac_t   = rompp::algebra::Matrix<app_jacob_t>;

  ode_state_t y(y0n);
  constexpr auto ode_case = rompp::ode::ImplicitEnum::Euler;
  using stepper_t = rompp::ode::ImplicitStepper<
    ode_case, ode_state_t, ode_res_t, ode_jac_t, app_t>;
  stepper_t stepperObj(y, appobj);

  // define solver
  using lin_solver_t = rompp::solvers::iterative::EigenIterative<
    rompp::solvers::linear::iterative::Bicgstab, ode_jac_t>;
  rompp::solvers::NewtonRaphson<scalar_t, lin_solver_t> solverO;
  solverO.setTolerance(1e-6);
  solverO.setMaxIterations(200);

  // integrate in time
  constexpr scalar_t dt = 0.0001;
  constexpr auto Nsteps = 10;
  constexpr scalar_t fint = Nsteps*dt;
  Observer obs;
  rompp::ode::integrateNSteps(stepperObj, y, 0.0, dt, Nsteps, obs, solverO);
  std::cout << std::fixed << std::setprecision(14) << *y.data() << std::endl;
  {
    using namespace rompp::apps::test;
    checkSol(y,
  	     NonLinAdvDiffReacFlame2dImpGoldStates<ode_case>::get(Nx, Ny, dt, fint));
  }

  std::cout << checkStr << std::endl;
  return 0;
}
