
#include "CORE_VECTOR"
#include "CORE_MATRIX"
#include "ODE_ALL"
#include "SOLVERS_EXP"
// app class
#include "apps_burgers1d_epetra.hpp"
#include "ode_observer.hpp"

#include "../apps_helper_ode.hpp"


struct mysizer{
  using state_t = core::vector<apps::burgers1dEpetra::state_type>;
 static size_t getSize(const state_t & obj){
   return obj.localSize();
 };
 static void matchSize(const state_t & src,
		       state_t & obj){
   obj.replaceDataMap(src.getDataMap());
 };
};

int main(int argc, char *argv[])
{
  using native_state_t = apps::burgers1dEpetra::state_type;
  using native_residual_t = native_state_t;
  using native_jac_t = apps::burgers1dEpetra::jacobian_type;
  using scalar_t = apps::burgers1dEpetra::scalar_type;
  using target_app_t = apps::burgers1dEpetra;

  MPI_Init(&argc,&argv);
  int rank; // My process ID
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  std::vector<double> mu({5.0, 0.02, 0.02});
  target_app_t appObj(mu, 20, &Comm);
  appObj.setup();
  
  // integrate in time starting from y0
  scalar_t final_t = 35;
  scalar_t dt = 0.001;

  // wrap with core structures
  using state_t = core::vector<native_state_t>;
  using jac_t = core::matrix<native_jac_t>;
  using residual_t = state_t;
  native_state_t y0n = appObj.getInitialState();
  snapshot_collector collObj;

  // data map from epetra
  auto const & dmap = y0n.Map();
    
  state_t y0(y0n);

  // using stepper_t = ode::explicitEulerStepper<
  //   state_t, residual_t, scalar_t,
  //   target_app_t, scalar_t, mysizer>;
  // stepper_t stepperObj(appObj, dmap);

  // ode::integrateNSteps(stepperObj, y, ti, dt, te/dt, collectorObj);
  
  using exstd_t = apps::test::appsExpOdeHelper<
    state_t, residual_t, scalar_t, target_app_t,
    scalar_t, mysizer, void, void, snapshot_collector>;

  exstd_t::explicitStandardEuler(y0, 0.0, final_t, dt, appObj, collObj, dmap);
  y0.data()->Print(std::cout << std::setprecision(14));
  // //printSol("Exp Euler", y0);
    
  MPI_Finalize();
  return 0;
}
