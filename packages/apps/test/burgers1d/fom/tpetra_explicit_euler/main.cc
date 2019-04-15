
#include "CORE_ALL"
#include "ODE_ALL"
#include "APPS_BURGERS1D"
#include "../gold_states_explicit.hpp"

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
  using app_t		= rompp::apps::Burgers1dTpetra;
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_residual_t	= typename app_t::residual_type;

  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;

  // MPI init
  MPI_Init(&argc,&argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  rcpcomm_t Comm = Teuchos::rcp (new tcomm_t(MPI_COMM_WORLD));
  assert( Comm->getSize() == 4);

  // scope guard needed for tpetra
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    std::vector<double> mu({5.0, 0.02, 0.02});
    const int Ncells = 20;
    app_t appobj(mu, Ncells, Comm);
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
    {
      using namespace rompp::apps::test;
      checkSol(rank, y, Burgers1dExpGoldStates<ode_case>::get(Ncells, dt, fint));
    }
  }

  MPI_Finalize();
  std::cout << checkStr << std::endl;
  return 0;
}
