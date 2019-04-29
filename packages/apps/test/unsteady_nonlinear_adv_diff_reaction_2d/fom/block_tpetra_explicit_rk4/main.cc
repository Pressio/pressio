
#include "CORE_ALL"
#include "ODE_ALL"
#include "APPS_UNSTEADYNONLINADVDIFFREACTION2D"
#include "../gold_states_explicit.hpp"
#include <Tpetra_Core.hpp>

constexpr double eps = 1e-12;
std::string checkStr {"PASSED"};

template <typename T>
void checkSol(T & y, //non Const because we need getVectorView
	      const std::vector<double> & trueS){
  auto y_tpv = y.data()->getVectorView();
  const auto arrrcp = y_tpv.getData();
  if (trueS.empty()) {
    std::cout << " true solution not found, empty " << std::endl;
    checkStr = "FAILED";
  }

  for (size_t i=0; i<trueS.size(); i++){
    std::cout << i << " "
	      << std::setprecision(14)
	      << arrrcp[i] << " " << trueS[i]
	      << std::endl;

    if (std::abs( arrrcp[i] - trueS[i]) > eps or
	std::isnan(arrrcp[i]) )
      checkStr = "FAILED";
  }
}


int main(int argc, char *argv[]){
  using app_t		= rompp::apps::UnsteadyNonLinAdvDiffReac2dBlockTpetra;
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_residual_t	= typename app_t::residual_type;
  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;
  constexpr auto zero = ::rompp::core::constants::zero<scalar_t>();

  // scope guard needed (MPI init within trilinos)
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rcpcomm_t Comm = Teuchos::rcp (new tcomm_t(MPI_COMM_WORLD));

    // create app object
    constexpr int Nx = 11, Ny = Nx*2-1;
    app_t appobj(Comm, Nx, Ny);
    appobj.setup();
    const auto y0n = appobj.getInitialState();
    const auto r0n = appobj.residual(y0n, zero);

    using ode_state_t = rompp::core::Vector<app_state_t>;
    using ode_res_t   = rompp::core::Vector<app_residual_t>;
    ode_state_t y(y0n);
    ode_res_t r(r0n);

    constexpr auto ode_case = rompp::ode::ExplicitEnum::RungeKutta4;
    using stepper_t = rompp::ode::ExplicitStepper<
      ode_case, ode_state_t, app_t, ode_res_t>;
    stepper_t stepperObj(y, appobj, r);

    // integrate in time
    constexpr scalar_t dt = 0.001;
    constexpr scalar_t fint = 0.5;
    constexpr auto Nsteps = static_cast<unsigned int>(fint/dt);
    rompp::ode::integrateNSteps(stepperObj, y, 0.0, dt, Nsteps);
    {
      using namespace rompp::apps::test;
      checkSol
	(y, NonLinAdvDiffReac2dExpGoldStates<ode_case>::get(Nx,
							    Ny,
							    dt,
							    fint));
    }
  }

  std::cout << checkStr << std::endl;
  return 0;
}
