
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "APPS_BURGERS1D"
#include "../gold_states_explicit.hpp"

constexpr double eps = 1e-12;
std::string checkStr {"PASSED"};

template <typename T>
void checkSol(int rank, const T & y,
	      const std::vector<double> & trueS){

  int shift = 0;
  if (rank==1) shift = 5;
  else if (rank==2) shift = 10;
  else if (rank==3) shift = 15;

  for (auto i=0; i<5; i++){
    if (std::abs(y[i] - trueS[i+shift]) > eps) checkStr = "FAILED";
  }
}

int main(int argc, char *argv[]){
  using app_t		= rompp::apps::Burgers1dEpetra;
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_residual_t	= typename app_t::residual_type;

  MPI_Init(&argc,&argv);
  int rank; // My process ID
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  assert( Comm.NumProc() == 4);

  //-------------------------------
  std::vector<double> mu({5.0, 0.02, 0.02});
  const int Ncells = 20;
  app_t appobj(mu, Ncells, &Comm);
  appobj.setup();
  auto & y0n = appobj.getInitialState();
  auto r0n = appobj.residual(y0n, static_cast<scalar_t>(0));

  using ode_state_t = rompp::containers::Vector<app_state_t>;
  using ode_res_t   = rompp::containers::Vector<app_residual_t>;
  ode_state_t y(y0n);

  constexpr auto ode_case = rompp::ode::ExplicitEnum::Euler;
  using stepper_t = rompp::ode::ExplicitStepper<
    ode_case, ode_state_t, app_t, ode_res_t, scalar_t>;
  stepper_t stepperObj(y, appobj);

  // integrate in time
  scalar_t fint = 35;
  scalar_t dt = 0.01;
  auto Nsteps = static_cast<unsigned int>(fint/dt);
  rompp::ode::integrateNSteps(stepperObj, y, 0.0, dt, Nsteps);
  y.data()->Print(std::cout << std::setprecision(14));
  {
    using namespace rompp::apps::test;
    checkSol(rank, y, 
             Burgers1dExpGoldStates<ode_case>::get(Ncells, dt, fint));
  }

  MPI_Finalize();
  std::cout << checkStr << std::endl;
  return 0;
}
