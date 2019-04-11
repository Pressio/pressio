
#include "CORE_ALL"
#include "ODE_ALL"
#include "APPS_UNSTEADYLINADVDIFF1D"
#include "../gold_states_explicit.hpp"

constexpr double eps = 1e-12;

template <typename T>
void checkSol(int rank, const T & y,
	      const std::vector<double> & trueS){
  if (rank == 0){
    assert(std::abs(y[0] - trueS[0]) < eps);
    assert(std::abs(y[1] - trueS[1]) < eps);
    assert(std::abs(y[2] - trueS[2]) < eps);
    assert(std::abs(y[3] - trueS[3]) < eps);
    assert(std::abs(y[4] - trueS[4]) < eps);
  }
  if (rank == 1){
    assert(std::abs(y[0] - trueS[5]) < eps);
    assert(std::abs(y[1] - trueS[6]) < eps);
    assert(std::abs(y[2] - trueS[7]) < eps);
    assert(std::abs(y[3] - trueS[8]) < eps);
  }
}

int main(int argc, char *argv[]){
  using app_t          =rompp::apps::UnsteadyLinAdvDiff1dEpetra;
  using scalar_t       =typename app_t::scalar_type;
  using app_state_t    =typename app_t::state_type;
  using app_residual_t =typename app_t::residual_type;  //Need residual
  using native_state   =typename app_t:: state_type;
  //----------------------------------------------------------------------
  // Initialize MPI
  //----------------------------------------------------------------------
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  assert(Comm.NumProc() == 2);
  //----------------------------------------------------------------------
  // Parameters
  //----------------------------------------------------------------------
  std::vector<scalar_t> mu{-1 , 1, 1};
  std::vector<scalar_t> domain{0, 2.0, 0.2};
  std::vector<scalar_t> bc1D{0, 0};

  //----------------------------------------------------------------------
  // Problem setup
  //----------------------------------------------------------------------
  app_t appObj(Comm, mu, domain, bc1D);
  appObj.unsteadySetup();
  //States and forcing term to calculate residual
  const native_state y0n(*appObj.getInitialState());
  //Spatial Residual
  auto r0n = appObj.residual(y0n, static_cast<scalar_t>(0));
  //----------------------------------------------------------------------
  // Rompp time integrator
  //----------------------------------------------------------------------
  using ode_state_t = rompp::core::Vector<app_state_t>;
  using ode_res_t  = rompp::core::Vector<app_residual_t>;
  ode_state_t y(y0n);
  ode_res_t r(r0n);
  y.data()->Print(std::cout <<std::setprecision(14));

  constexpr auto ode_case = rompp::ode::ExplicitEnum::Euler;
  using stepper_t = rompp::ode::ExplicitStepper<
    ode_case, ode_state_t, app_t, ode_res_t>;
  stepper_t stepperObj(y, appObj, r);

  //Integrate in time
  scalar_t fint = 1;
  scalar_t dt = 0.01;
  auto Nsteps = static_cast<unsigned int>(fint/dt);
  rompp::ode::integrateNSteps(stepperObj, y, 0.0, dt, Nsteps);

  //----------------------------------------------------------------------
  // Results
  //----------------------------------------------------------------------
  y.data()->Print(std::cout <<std::setprecision(14));
  {
    using namespace rompp::apps::test;
    checkSol(rank, y,
	     UnsteadyLinAdvDiff1dExpGoldStates<ode_case>
	     ::get(mu, domain, bc1D, dt, fint));
  }

  MPI_Finalize();
  return 0;
}
