
#include "UTILS_ALL"
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "SOLVERS_EXPERIMENTAL"
#include "ROM_WLS"
#include "APPS_UNSTEADYBURGERS1D"
#include "utils_kokkos.hpp"

int main(int argc, char *argv[]){
  using fom_t		= pressio::apps::Burgers1dKokkos;
  using scalar_t	= typename fom_t::scalar_type;

  // using exe_space = typename fom_t::execution_space;
  using native_state_t_d = typename fom_t::state_type_d;
  using native_state_t_h = typename fom_t::state_type_h;
  using native_mv_t_d = typename fom_t::mv_d;


  using fom_state_t      = ::pressio::containers::Vector<native_state_t_d>;

  // using native_mv_t_h = typename fom_t::mv_h;

  // device wls state type
  using wls_state_d_t	= pressio::containers::Vector<native_state_t_h>;

  // device decoder jacobian type
  using decoder_jac_d_t	= pressio::containers::MultiVector<native_mv_t_d>;
  // host decoder jacobian type
  // using decoder_jac_h_t	= pressio::containers::MultiVector<native_mv_t_h>;

  // device decoder type
  using decoder_d_t	= pressio::rom::LinearDecoder<decoder_jac_d_t>;
  using hessian_t  = pressio::containers::Matrix<typename fom_t::mv_d>;

  std::string checkStr {"PASSED"};
  constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();

  Kokkos::initialize (argc, argv);
  {
    // app object
    constexpr int numCell = 20;
    fom_t appObj({{5.0, 0.02, 0.02}}, numCell);
    constexpr scalar_t dt = 0.01;
    constexpr auto t0 = zero;

    constexpr int romSize = 11;

    // create/read jacobian of the decoder
    decoder_jac_d_t phi("phi", numCell, 11);
    pressio::rom::test::kokkos::readBasis("basis.txt", romSize, numCell, *phi.data());
    if( phi.numVectors() != romSize ) return 0;

    // create decoder obj
    decoder_d_t decoderObj(phi);

    // for this problem, my reference state = initial state
    // get initial condition
    auto & yFOM_IC_native = appObj.getInitialState();

    // wrap into pressio container
    fom_state_t yFOM_IC(yFOM_IC_native);
    //reference state is equal to the IC
    fom_state_t & yRef = yFOM_IC;

    // -----------------
    // L solver
    // -----------------
    using lin_solver_tag  = pressio::solvers::linear::direct::getrs;
    using linear_solver_t = pressio::solvers::direct::KokkosDirect<lin_solver_tag, hessian_t>;
    linear_solver_t linear_solver;

    // -----------------
    // WLS problem
    // -----------------
    constexpr int numStepsInWindow = 5;
    using ode_tag	     = ::pressio::ode::implicitmethods::BDF2;
    using wls_system_t = pressio::rom::wls::SystemHessianAndGradientApi<fom_t,wls_state_d_t,decoder_d_t,ode_tag,hessian_t,linear_solver_t>;

    // create the wls state
    wls_state_d_t  wlsState("yRom",romSize*numStepsInWindow); wlsState.setZero();
    // create the wls system
    wls_system_t wlsSystem(appObj, yFOM_IC, yRef, decoderObj, numStepsInWindow,romSize,linear_solver);

    // -----------------
    // NL solver
    // -----------------
    using lin_solver_tag  = pressio::solvers::linear::direct::getrs;
    using linear_solver_t = pressio::solvers::direct::KokkosDirect<lin_solver_tag, hessian_t>;
    using gn_t            = pressio::solvers::iterative::GaussNewton<linear_solver_t, wls_system_t>;
    linear_solver_t linear_solver;
    gn_t GNSolver(wlsSystem, wlsState, linear_solver);
    GNSolver.setTolerance(1e-13);
    GNSolver.setMaxIterations(50);

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
    fom_state_t yFinal("yFF_d",numCell); //may not build

    using fom_state_reconstr_t = pressio::rom::FomStateReconstructor<fom_state_t, decoder_d_t>;
    fom_state_reconstr_t fomStateReconstructor(yRef, decoderObj);
    fomStateReconstructor(wlsCurrentState, yFinal);

    // create a host mirror for yFinal
    native_state_t_h yFinal_h("yFF_h", numCell);
    Kokkos::deep_copy(yFinal_h, *yFinal.data());

    // get true solution
    const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF2::get(numCell, dt, finalTime);
    for (auto i=0; i<numCell; i++){
      std::cout << std::setprecision(15) << yFinal_h(i) << " " << trueY[i] << std::endl;
      if (std::abs(yFinal_h(i) - trueY[i]) > 1e-10) checkStr = "FAILED";
    }
    std::cout << checkStr << std::endl;
  }
  Kokkos::finalize();
  return 0;
}
