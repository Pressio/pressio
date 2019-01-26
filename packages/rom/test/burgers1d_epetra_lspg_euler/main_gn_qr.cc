
#include "CORE_ALL"
#include "ODE_ALL"
#include "SVD_BASIC"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "QR_BASIC"

#include "oneD_burgers1d_epetra.hpp"
#include "state_collector.hpp"
#include "../../burgers1d_gold_states.hpp"
#include "../../burgers1d_epetra_observer.hpp"
#include "../lspg_utils.hpp"
#include "../fom_gold.hpp"


int main(int argc, char *argv[]){

  //-------------------------------
  // declare native app types
  using app_t		= rompp_demo_apps::Burgers1dEpetra;
  using app_state_t	= app_t::state_type;
  using app_space_res_t	= app_t::residual_type;
  using scalar_t	= app_t::scalar_type;

  // declare wrapper types
  using app_state_w_t	= rompp::core::Vector<app_state_t>;
  using app_res_w_t	= rompp::core::Vector<app_space_res_t>;

  // declare type of POD modes
  using phi_t	 = rompp::core::MultiVector<Epetra_MultiVector>;
  // declare type of operator wrapping the POD modes
  using phi_op_t = rompp::rom::MultiVectorOperator<phi_t>;

  // declare ROM types
  using eig_stat_vec = Eigen::Matrix<scalar_t, -1, 1>;
  using rom_state_t  = rompp::core::Vector<eig_stat_vec>;
  using rom_res_t    = app_res_w_t;
  using rom_jac_t    = rompp::core::MultiVector<Epetra_MultiVector>;

  // declare residual policy type
  using res_pol_t = rompp::rom::RomLSPGResidualPolicy<
			app_state_w_t, app_res_w_t, phi_op_t, 1, 1>;
  using jac_pol_t = rompp::rom::RomLSPGJacobianPolicy<
			app_state_w_t, rom_jac_t, phi_op_t, 1>;

  // stepper needs to know about ROM types for doing time integration,
  using rom_stepper_t = rompp::ode::ImplicitStepper<
			rompp::ode::ImplicitEnum::Euler,
			rom_state_t, rom_res_t, rom_jac_t,
			app_t, void, res_pol_t, jac_pol_t>;

  // declare QR solver types
  using qr_algo = rompp::qr::Householder;
  using qr_type = rompp::qr::QRSolver<rom_jac_t, qr_algo>;

  // declare Gauss-newton solver
  using converged_when_t = rompp::solvers::iterative::default_convergence;
  using gnsolver_t	 = rompp::solvers::iterative::GaussNewtonQR<
				scalar_t, qr_type,
				converged_when_t, rom_stepper_t>;

  // collector for storing states during time integration
  using collector_t = state_collector<rom_state_t>;


  //-------------------------------
  // MPI init
  MPI_Init(&argc,&argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  assert(Comm.NumProc() == 2);

  //-------------------------------
  // app object
  int numCell = 50; // # fv cells
  std::vector<double> mu({5.0, 0.02, 0.02});
  app_t appobj(mu, numCell, &Comm);
  appobj.setup();
  auto t0 = static_cast<scalar_t>(0);
  auto & y0n = appobj.getInitialState();
  auto r0n = appobj.residual(y0n, t0);
  scalar_t dt = 0.01;

  // POD modes
  constexpr int romSize = 20;
  phi_t phi = demo_apps::test::getBasis(romSize, numCell, Comm,
					appobj.getDataMap());
  int numBasis = phi.globalNumVectors();
  assert( numBasis = romSize );

  //----------------------------------------
  // ROM
  //----------------------------------------
  // define operator wrapping the modes
  phi_op_t phiOp(phi);

  // define y0 and R0 for FOM
  app_state_w_t y0FOM(y0n);
  app_res_w_t r0FOM(r0n);

  // define residual and jacobian policies
  res_pol_t resObj(y0FOM, r0FOM, phiOp);
  jac_pol_t jacObj(y0FOM, phiOp);

  // define ROM state
  rom_state_t yROM(romSize);
  // initialize to zero (this has to be done)
  yROM.putScalar(0.0);

  // define stepper
  rom_stepper_t romStepperObj(appobj, resObj, jacObj, yROM);

  // define solver
  gnsolver_t solver(romStepperObj, yROM);
  solver.setTolerance(1e-14);
  solver.setMaxIterations(500);

  // integrate in time
  collector_t stateStore;
  rompp::ode::integrateNSteps(romStepperObj, yROM,
			      0.0, dt, 200,
			      stateStore, solver);

  // done with time integration

  //if (rank==0)  stateStore.printAll();

  // compute the fom corresponding to our rom final state
  app_state_w_t yRf_ = y0FOM;
  yRf_ += phiOp.apply(yROM);
  yRf_.data()->Print(std::cout << std::setprecision(14));

  // check against gold solution
  int shift = 0;
  if (rank==1)  shift = 25;
  int myn = yRf_.getDataMap().NumMyElements();
  for (auto i=0; i<myn; i++){
    assert(std::abs(yRf_[i] -
	   demo_apps::test::burger1DGoldDt001step200().bdf1Sol[i+shift])
	   < 1e-12 );
   }

  // print summary from timers
  #ifdef HAVE_TEUCHOS_TIMERS
  rompp::core::TeuchosPerformanceMonitor::stackedTimersReportMPI();
  #endif

  MPI_Finalize();
  return 0;
}
//------------------------------
// end
//------------------------------
