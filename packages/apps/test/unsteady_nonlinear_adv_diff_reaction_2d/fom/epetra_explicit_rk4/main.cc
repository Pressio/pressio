
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "APPS_UNSTEADYNONLINADVDIFFREACTION2D"
#include "../gold_states_explicit.hpp"

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
    const auto err = std::abs(y[i] - trueS[i]);
    std::cout << std::fixed << std::setprecision(15)
	      << " true = " << trueS[i]
	      << " y = " << y[i]
	      << " err = " << err
	      << std::endl;
    if ( err > eps or std::isnan(y[i])) checkStr = "FAILED";
  }
}

int main(int argc, char *argv[]){
  using app_t		= rompp::apps::UnsteadyNonLinAdvDiffReac2dEpetra;
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_residual_t	= typename app_t::residual_type;
  constexpr auto zero = ::rompp::utils::constants::zero<scalar_t>();

  MPI_Init(&argc,&argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  if(Comm.NumProc() != 1){
    checkStr = "FAIL";
    return 0;
  }

  constexpr int Nx = 11, Ny = Nx*2-1;
  app_t appobj(Comm, Nx, Ny);
  appobj.setup();
  const auto y0n = appobj.getInitialState();

  using ode_state_t = rompp::containers::Vector<app_state_t>;
  using ode_res_t   = rompp::containers::Vector<app_residual_t>;
  ode_state_t y(y0n);

  constexpr auto ode_case = rompp::ode::ExplicitEnum::RungeKutta4;
  using stepper_t = rompp::ode::ExplicitStepper<
    ode_case, ode_state_t, app_t, ode_res_t, scalar_t>;
  stepper_t stepperObj(y, appobj);

  // integrate in time
  constexpr scalar_t dt = 0.001;
  constexpr auto Nsteps = static_cast<unsigned int>(500);
  constexpr scalar_t fint = Nsteps*dt;
  rompp::ode::integrateNSteps(stepperObj, y, 0.0, dt, Nsteps);
  y.data()->Print(std::cout << std::setprecision(14));
  {
    using namespace rompp::apps::test;
    checkSol
      (y, NonLinAdvDiffReac2dExpGoldStates<ode_case>::get(Nx, Ny, dt, fint));
  }

  MPI_Finalize();
  std::cout << checkStr << std::endl;
  return 0;
}
