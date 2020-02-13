
#include "ROM_BASIC"
#include "APPS_UNSTEADYBURGERS1D"
#include "utils_tpetra.hpp"


int main(int argc, char *argv[]){
  std::string checkStr {"PASSED"};

  using fom_t		   = pressio::apps::Burgers1dTpetra;
  using scalar_t	   = typename fom_t::scalar_type;
  using fom_native_state_t = typename fom_t::state_type;
  using fom_dmat_t         = typename fom_t::dense_matrix_type;
  using fom_state_t        = ::pressio::containers::Vector<fom_native_state_t>;

  using eig_dyn_vec	   = Eigen::Matrix<scalar_t, -1, 1>;
  using eig_dyn_mat	   = Eigen::Matrix<scalar_t, -1, -1>;
  using wls_state_t	   = pressio::containers::Vector<eig_dyn_vec>;
  using hessian_t          = pressio::containers::Matrix<eig_dyn_mat>;

  using decoder_jac_t      = pressio::containers::MultiVector<fom_dmat_t>;
  using decoder_t	   = pressio::rom::LinearDecoder<decoder_jac_t, wls_state_t, fom_state_t>;

  using tcomm_t		   = Teuchos::MpiComm<int>;
  using rcpcomm_t	   = Teuchos::RCP<const tcomm_t>;

  constexpr auto zero = pressio::utils::constants::zero<scalar_t>();

  // scope guard needed for tpetra
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rcpcomm_t Comm = Teuchos::rcp (new tcomm_t(MPI_COMM_WORLD));

    // -------------------
    // fom object
    // -------------------
    constexpr int fomSize = 20;
    fom_t appObj( {5.0, 0.02, 0.02}, fomSize, Comm);
    // get initial condition
    auto & yFOM_IC_native = appObj.getInitialState();
    // wrap into pressio container
    fom_state_t yFOM_IC(yFOM_IC_native);
    //reference state is equal to the IC
    fom_state_t & yRef = yFOM_IC;

    // -------------------
    // decoder
    // -------------------
    int romSize = 11;
    const auto phiNative = pressio::rom::test::tpetra::readBasis("basis.txt", romSize, fomSize,
								 Comm, appObj.getDataMap());
    decoder_t decoderObj(phiNative);

    // -----------------
    // lin solver
    // -----------------
    using lin_solver_tag  = pressio::solvers::linear::direct::ColPivHouseholderQR;
    using linear_solver_t = pressio::solvers::direct::EigenDirect<lin_solver_tag, hessian_t>;
    linear_solver_t linear_solver;

    // -----------------
    // WLS problem
    // -----------------
    constexpr int numStepsInWindow = 5;
    using ode_tag	     = ::pressio::ode::implicitmethods::Euler;
    using wls_system_t = pressio::rom::wls::SystemHessianAndGradientApi<fom_t,wls_state_t,decoder_t,ode_tag,hessian_t>;
    // create the wls state
    wls_state_t  wlsState(romSize*numStepsInWindow);
    pressio::containers::ops::fill(wlsState, zero);
    // create the wls system
    wls_system_t wlsSystem(appObj, yFOM_IC, yRef, decoderObj, numStepsInWindow, romSize, linear_solver);

    // -----------------
    // solver
    // -----------------
    using gn_t		= pressio::solvers::iterative::GaussNewton<linear_solver_t, wls_system_t>;
    gn_t GNSolver(wlsSystem, wlsState, linear_solver);
    GNSolver.setTolerance(1e-13);
    GNSolver.setMaxIterations(50);

    // -----------------
    // solve wls problem
    // -----------------
    constexpr scalar_t finalTime = 0.1;
    constexpr scalar_t dt	       = 0.01;
    constexpr int numWindows     = static_cast<int>(finalTime/dt)/numStepsInWindow;

    auto startTime = std::chrono::high_resolution_clock::now();
    for (auto iWind = 0; iWind < numWindows; iWind++){
      wlsSystem.advanceOneWindow(wlsState, GNSolver, iWind, dt);
    }

    const auto finishTime = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> elapsed = finishTime - startTime;
    std::cout << "Walltime = " << elapsed.count() << '\n';

    // -----------------
    // process solution
    // -----------------
    const auto wlsCurrentState = pressio::containers::span(wlsState, (numStepsInWindow-1)*romSize, romSize);
    fom_state_t yFinal(yFOM_IC_native);
    using fom_state_reconstr_t = pressio::rom::FomStateReconstructor<scalar_t, fom_state_t, decoder_t>;
    fom_state_reconstr_t fomStateReconstructor(yRef, decoderObj);
    fomStateReconstructor(wlsCurrentState, yFinal);

    // this is a reproducing ROM test, so the final reconstructed state
    // has to match the FOM solution obtained with euler, same time-step, for 10 steps
    auto yFF_v = yFinal.data()->getData();
    int shift = (rank==0) ? 0 : 10;
    const int myn = yFinal.data()->getMap()->getNodeNumElements();
    const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(fomSize, dt, 0.10);
    for (auto i=0; i<myn; i++)
      if (std::abs(yFF_v[i] - trueY[i+shift]) > 1e-10) checkStr = "FAILED";

    std::cout << checkStr << std::endl;
  }

  return 0;
}
