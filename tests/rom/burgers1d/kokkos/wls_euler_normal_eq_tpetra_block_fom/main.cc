
#include "pressio_rom.hpp"
#include "pressio_apps.hpp"
#include "utils_kokkos.hpp"
#include "utils_tpetra.hpp"

namespace {

template <
  typename tcomm_t,
  typename lin_solver_tag,
  typename hessian_matrix_structure_tag,
  typename rcpcomm_t
  >
void doRun(rcpcomm_t & Comm, int rank)
{

  using fom_t		= pressio::apps::Burgers1dTpetraBlock;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t= typename fom_t::state_type;
  //using native_dmat_d_t = typename fom_t::dense_matrix_type;
  using fom_dmat_t         = typename fom_t::dense_matrix_type;
  using fom_state_t   = ::pressio::containers::Vector<native_state_t>;

  // wls state type
  using execution_space = Kokkos::DefaultExecutionSpace;
  using kll = Kokkos::LayoutLeft;
  using k1dLl_d = Kokkos::View<scalar_t*, kll, execution_space>;
  using k2dLl_d = Kokkos::View<scalar_t**, kll, execution_space>;

  using wls_state_d_t	= pressio::containers::Vector<k1dLl_d>;
  using hessian_d_t	= pressio::containers::Matrix<k2dLl_d>;

  // decoder jacobian type
  using decoder_jac_d_t	= pressio::containers::MultiVector<fom_dmat_t>;
  using decoder_d_t	= pressio::rom::LinearDecoder<decoder_jac_d_t, wls_state_d_t, fom_state_t>;

  std::string checkStr {"PASSED"};
  constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();

  // app object
  constexpr int numCell = 20;
  fom_t appObj({{5.0, 0.02, 0.02}}, numCell,Comm);
  constexpr scalar_t dt = 0.01;
  constexpr auto t0 = zero;

  int romSize = 11;

  // create/read jacobian of the decoder
  auto tpw_phi = pressio::rom::test::tpetra::readBasis("basis.txt", romSize,
                 numCell, Comm, appObj.getDataMap());
  fom_dmat_t tpb_phi(*tpw_phi.data(), *appObj.getDataMap(), 1);
  decoder_d_t decoderObj(tpb_phi);



  // for this problem, my reference state = initial state
  // get initial condition
  auto & yFOM_IC_native = appObj.getInitialState();

  // wrap into pressio container
  fom_state_t yFOM_IC(yFOM_IC_native);
  //reference state is equal to the IC
  fom_state_t & yRef = yFOM_IC;

  // -----------------
  // lin solver
  // -----------------
  using linear_solver_t = pressio::solvers::direct::KokkosDirect<lin_solver_tag, hessian_d_t>;
  linear_solver_t linear_solver;

  // -----------------
  // WLS problem
  // -----------------
  constexpr int numStepsInWindow = 5;
  using ode_tag      = ::pressio::ode::implicitmethods::Euler;
  using wls_system_t = pressio::rom::wls::SystemHessianAndGradientApi<fom_t, wls_state_d_t, decoder_d_t,
								      ode_tag, hessian_d_t,
								      hessian_matrix_structure_tag>;
  // create the wls system
  wls_system_t wlsSystem(appObj, yFOM_IC, yRef, decoderObj, numStepsInWindow,romSize,linear_solver);


  // create the wls state
  wls_state_d_t  wlsState(romSize*numStepsInWindow);
  pressio::ops::set_zero(wlsState);

  // -----------------
  // NL solver
  // -----------------
  using gn_t            = pressio::solvers::iterative::GaussNewton<linear_solver_t, wls_system_t>;
  gn_t GNSolver(wlsSystem, wlsState, linear_solver);
  GNSolver.setTolerance(1e-13);
  GNSolver.setMaxIterations(5);

  // -----------------
  // solve wls problem
  // -----------------
  constexpr scalar_t finalTime = 0.1;
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
  using fom_state_reconstr_t = pressio::rom::FomStateReconstructor<scalar_t, fom_state_t, decoder_d_t>;
  fom_state_reconstr_t fomStateReconstructor(yRef, decoderObj);
  fomStateReconstructor(wlsCurrentState, yFinal);

  // this is a reproducing ROM test, so the final reconstructed state
  // has to match the FOM solution obtained with euler, same time-step, for 10 steps
  auto yFF_v = yFinal.data()->getVectorView().getData();
  int shift = (rank==0) ? 0 : 10;
  const int myn = yFinal.data()->getMap()->getNodeNumElements();
  const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(numCell, dt, 0.10);
  for (auto i=0; i<myn; i++)
    if (std::abs(yFF_v[i] - trueY[i+shift]) > 1e-10) checkStr = "FAILED";

  std::cout << checkStr << std::endl;

}
}// end namespace


int main(int argc, char *argv[])
{
  using tcomm_t		   = Teuchos::MpiComm<int>;
  using rcpcomm_t	   = Teuchos::RCP<const tcomm_t>;


  // scope guard needed for tpetra
  Kokkos::initialize (argc, argv);
  {
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rcpcomm_t Comm = Teuchos::rcp (new tcomm_t(MPI_COMM_WORLD));
    //doRun< tcomm_t, pressio::solvers::linear::direct::potrsU, pressio::matrixUpperTriangular >(Comm, rank);
    doRun< tcomm_t, pressio::solvers::linear::direct::potrsL, pressio::matrixLowerTriangular >(Comm, rank);
  }
  Kokkos::finalize();
  }
  return 0;
}
