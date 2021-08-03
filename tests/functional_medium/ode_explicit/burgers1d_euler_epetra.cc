
#include "pressio_ode_explicit.hpp"
#include "../testing_apps/apps.hpp"

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

int main(int argc, char *argv[])
{
  using app_t		= pressio::apps::Burgers1dEpetra;
  using scalar_t	= typename app_t::scalar_type;
  using state_t	= typename app_t::state_type;

  MPI_Init(&argc,&argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  assert( Comm.NumProc() == 4);

  std::vector<double> mu({5.0, 0.02, 0.02});
  const int Ncells = 20;
  app_t appObj(mu, Ncells, &Comm);
  auto & y0n = appObj.getInitialState();
  state_t y(y0n);

  auto stepperObj = pressio::ode::create_forward_euler_stepper(y, appObj);

  scalar_t fint = 35;
  scalar_t dt = 0.01;
  auto Nsteps = static_cast<::pressio::ode::step_count_type>(fint/dt);
  pressio::ode::advance_n_steps(stepperObj, y, 0.0, dt, Nsteps);
  // y.Print(std::cout << std::setprecision(14));
  {
    using namespace pressio::apps::test;
    checkSol(rank, y,
      Burgers1dExpGoldStatesEuler::get(Ncells, dt, fint));
  }

  MPI_Finalize();
  std::cout << checkStr << std::endl;
  return 0;
}
