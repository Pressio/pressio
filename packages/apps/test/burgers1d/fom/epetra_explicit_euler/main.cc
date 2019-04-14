
#include "CORE_ALL"
#include "ODE_ALL"
#include "APPS_BURGERS1D"
#include "../gold_states_explicit.hpp"

constexpr double eps = 1e-12;

template <typename T>
void checkSol(int rank, const T & y,
	      const std::vector<double> & trueS){
  if (rank==0){
    assert(std::abs(y[0] - trueS[0]) < eps);
    assert(std::abs(y[1] - trueS[1]) < eps);
    assert(std::abs(y[2] - trueS[2]) < eps);
    assert(std::abs(y[3] - trueS[3]) < eps);
    assert(std::abs(y[4] - trueS[4]) < eps);}

  if (rank==1){
    assert(std::abs(y[0] - trueS[5]) < eps);
    assert(std::abs(y[1] - trueS[6]) < eps);
    assert(std::abs(y[2] - trueS[7]) < eps);
    assert(std::abs(y[3] - trueS[8]) < eps);
    assert(std::abs(y[4] - trueS[9]) < eps);}

  if (rank==2){
    assert(std::abs(y[0] - trueS[10]) < eps);
    assert(std::abs(y[1] - trueS[11]) < eps);
    assert(std::abs(y[2] - trueS[12]) < eps);
    assert(std::abs(y[3] - trueS[13]) < eps);
    assert(std::abs(y[4] - trueS[14]) < eps);}

  if (rank==3){
    assert(std::abs(y[0] - trueS[15]) < eps);
    assert(std::abs(y[1] - trueS[16]) < eps);
    assert(std::abs(y[2] - trueS[17]) < eps);
    assert(std::abs(y[3] - trueS[18]) < eps);
    assert(std::abs(y[4] - trueS[19]) < eps);}
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

  using ode_state_t = rompp::core::Vector<app_state_t>;
  using ode_res_t   = rompp::core::Vector<app_residual_t>;
  ode_state_t y(y0n);
  ode_res_t r(r0n);

  constexpr auto ode_case = rompp::ode::ExplicitEnum::Euler;
  using stepper_t = rompp::ode::ExplicitStepper<
    ode_case, ode_state_t, app_t, ode_res_t>;
  stepper_t stepperObj(y, appobj, r);

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
  std::cout << "PASSED" << std::endl;
  return 0;
}
