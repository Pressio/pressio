
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "APPS_UNSTEADYNONLINADVDIFFREACTIONFLAME2D"
#include "../gold_states_explicit.hpp"

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

constexpr bool do_print = true;

template <typename T>
struct Observer{
  Observer() = default;

  void operator()(size_t step,
  		  double t,
  		  const T & y)
  {
    if (do_print){
      if (step % 50 == 0){
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
  using app_t		= pressio::apps::UnsteadyNonLinAdvDiffReacFlame2dEigen;
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;
  // constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();

  constexpr int Nx = 12, Ny = 6;
  app_t appobj(Nx, Ny);
  appobj.setup();
  auto y0n = appobj.getInitialState();

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

  using ode_state_t = pressio::containers::Vector<app_state_t>;
  using ode_res_t   = pressio::containers::Vector<app_rhs_t>;
  ode_state_t y(y0n);

  constexpr auto ode_case = pressio::ode::ExplicitEnum::RungeKutta4;
  using stepper_t = pressio::ode::ExplicitStepper<
    ode_case, ode_state_t, app_t, ode_res_t, scalar_t>;
  stepper_t stepperObj(y, appobj);

  // integrate in time
  constexpr scalar_t dt = 0.00001;
  constexpr auto Nsteps = 1100;
  constexpr scalar_t fint = Nsteps*dt;
  Observer<app_state_t> obs;
  pressio::ode::integrateNSteps(stepperObj, y, 0.0, dt, Nsteps, obs);
  std::cout << std::fixed << std::setprecision(15) << *y.data() << std::endl;
  {
    using namespace pressio::apps::test;
    checkSol(y,
  	     NonLinAdvDiffReacFlame2dExpGoldStates<ode_case>::get(Nx,
								  Ny,
								  dt,
								  fint));
  }

  std::cout << checkStr << std::endl;
  return 0;
}
