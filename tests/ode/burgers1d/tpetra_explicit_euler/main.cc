
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "APPS_UNSTEADYBURGERS1D"

constexpr double eps = 1e-12;
std::string checkStr {"PASSED"};

template <typename T>
void checkSol(int rank, const T & y,
	      const std::vector<double> & trueS){
  auto y_v = y.data()->getData();

  int shift = 0;
  if (rank==1) shift = 5;
  else if (rank==2) shift = 10;
  else if (rank==3) shift = 15;

  for (auto i=0; i<5; i++){
    if (std::abs(y_v[i] - trueS[i+shift]) > eps)
      checkStr = "FAILED";
  }
}

int main(int argc, char *argv[]){
  using app_t		= pressio::apps::Burgers1dTpetra;
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_velocity_t	= typename app_t::velocity_type;

  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;

  // scope guard needed (MPI init within trilinos)
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rcpcomm_t Comm = Teuchos::rcp (new tcomm_t(MPI_COMM_WORLD));
    if (Comm->getSize() != 4){
      checkStr = "FAILED";
      return 0;
    }

    std::vector<double> mu({5.0, 0.02, 0.02});
    const int Ncells = 20;
    app_t appobj(mu, Ncells, Comm);
    appobj.setup();
    auto & y0n = appobj.getInitialState();

    using ode_state_t = pressio::containers::Vector<app_state_t>;
    using ode_res_t   = pressio::containers::Vector<app_velocity_t>;
    ode_state_t y(y0n);

    constexpr auto ode_case = pressio::ode::ExplicitEnum::Euler;
    using stepper_t = pressio::ode::ExplicitStepper<
      ode_case, ode_state_t, app_t, ode_res_t, scalar_t>;
    stepper_t stepperObj(y, appobj);

    // integrate in time
    scalar_t fint = 35;
    scalar_t dt = 0.01;
    auto Nsteps = static_cast<::pressio::ode::types::step_t>(fint/dt);
    pressio::ode::integrateNSteps(stepperObj, y, 0.0, dt, Nsteps);
    {
      using namespace pressio::apps::test;
      checkSol(rank, y, Burgers1dExpGoldStatesEuler::get(Ncells, dt, fint));
    }
  }

  std::cout << checkStr << std::endl;
  return 0;
}
