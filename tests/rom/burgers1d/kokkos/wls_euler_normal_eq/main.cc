
#include "pressio_rom.hpp"
#include "pressio_apps.hpp"
#include "utils_kokkos.hpp"

int main(int argc, char *argv[]){

  using fom_t		= pressio::apps::Burgers1dKokkos;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_d_t= typename fom_t::state_type;
  using native_dmat_d_t = typename fom_t::dense_matrix_type;

  using fom_state_d_t   = ::pressio::containers::Vector<native_state_d_t>;

  // wls state type
  using wls_state_d_t	= pressio::containers::Vector<native_state_d_t>;
  using hessian_d_t	= pressio::containers::Matrix<native_dmat_d_t>;

  // decoder jacobian type
  using decoder_jac_d_t	= pressio::containers::MultiVector<native_dmat_d_t>;
  using decoder_d_t	= pressio::rom::LinearDecoder<decoder_jac_d_t, wls_state_d_t, fom_state_d_t>;

  std::string checkStr {"PASSED"};
  constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();

  Kokkos::initialize (argc, argv);
  {
    // app object
    constexpr int numCell = 20;
    fom_t appObj({{5.0, 0.02, 0.02}}, numCell);
    constexpr scalar_t dt = 0.01;
    constexpr auto t0 = zero;

    int romSize = 11;

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
    fom_state_d_t yFOM_IC(yFOM_IC_native);
    //reference state is equal to the IC
    fom_state_d_t & yRef = yFOM_IC;

    // -----------------
    // lin solver
    // -----------------
    using lin_solver_tag  = pressio::solvers::linear::direct::potrsL;
    using linear_solver_t = pressio::solvers::direct::KokkosDirect<lin_solver_tag, hessian_d_t>;
    linear_solver_t linear_solver;

    // -----------------
    // WLS problem
    // -----------------
    constexpr int numStepsInWindow = 5;
    using ode_tag      = ::pressio::ode::implicitmethods::Euler;
    using hessian_matrix_structure_tag = pressio::matrixLowerTriangular;
    using wls_system_t = pressio::rom::wls::SystemHessianAndGradientApi<fom_t,wls_state_d_t,decoder_d_t,ode_tag,hessian_d_t, hessian_matrix_structure_tag>;
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
    fom_state_d_t yFinal("yFF_d",numCell); //may not build

    using fom_state_reconstr_t = pressio::rom::FomStateReconstructor<scalar_t, fom_state_d_t, decoder_d_t>;
    fom_state_reconstr_t fomStateReconstructor(yRef, decoderObj);
    fomStateReconstructor(wlsCurrentState, yFinal);

    // create a host mirror for yFinal
    using native_state_t_h = typename fom_t::state_type_h;
    native_state_t_h yFinal_h("yFF_h", numCell);
    Kokkos::deep_copy(yFinal_h, *yFinal.data());

    // get true solution
    const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(numCell, dt, finalTime);
    for (auto i=0; i<numCell; i++){
      std::cout << std::setprecision(15) << yFinal_h(i) << " " << trueY[i] << std::endl;
      if (std::abs(yFinal_h(i) - trueY[i]) > 1e-8) checkStr = "FAILED";
    }
    std::cout << checkStr << std::endl;
  }
  Kokkos::finalize();
  return 0;
}
