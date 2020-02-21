
#include "pressio_rom.hpp"
#include "pressio_apps.hpp"
#include "utils_eigen.hpp"

int main(int argc, char *argv[]){
  std::string checkStr {"PASSED"};

  using fom_t		   = pressio::apps::Burgers1dEigen;
  using scalar_t	   = typename fom_t::scalar_type;
  using fom_native_state_t = typename fom_t::state_type;
  using fom_dense_matrix_t = typename fom_t::dense_matrix_type;

  using fom_state_t        = ::pressio::containers::Vector<fom_native_state_t>;

  using eig_dyn_mat	   = Eigen::Matrix<scalar_t, -1, -1>;
  using eig_dyn_vec	   = Eigen::Matrix<scalar_t, -1, 1>;
  using wls_state_t	   = pressio::containers::Vector<eig_dyn_vec>;
  using hessian_t          = pressio::containers::Matrix<eig_dyn_mat>;
  using decoder_jac_t	   = pressio::containers::MultiVector<fom_dense_matrix_t>;
  using decoder_t	   = pressio::rom::LinearDecoder<decoder_jac_t, wls_state_t, fom_state_t>;

  constexpr auto zero = pressio::utils::constants::zero<scalar_t>();

  // -------------------
  // fom object
  // -------------------
  constexpr int fomSize = 20;
  fom_t appObj( Eigen::Vector3d{5.0, 0.02, 0.02}, fomSize);
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
  const auto phiNative = pressio::rom::test::eigen::readBasis("basis.txt", romSize, fomSize);
  decoder_t decoderObj(phiNative);

  // -----------------
  // linear solver
  // -----------------
  using lin_solver_tag	= pressio::solvers::linear::direct::ColPivHouseholderQR;
  using linear_solver_t = pressio::solvers::direct::EigenDirect<lin_solver_tag, hessian_t>;
  linear_solver_t linearSolver;

  // -----------------
  // WLS problem
  // -----------------
  constexpr int numStepsInWindow = 5;
  using ode_tag	     = ::pressio::ode::implicitmethods::Euler;

  // The initial condition is based on L^2 projection onto trial space
  wls_state_t wlsStateIc(romSize);
  const auto & decoderJac = decoderObj.getReferenceToJacobian();
  pressio::rom::utils::set_gen_coordinates_L2_projection<scalar_t>(linearSolver, decoderJac, yFOM_IC, yRef, wlsStateIc);

  // create the wls system
  using wls_system_t = pressio::rom::wls::SystemHessianAndGradientApi<fom_t,wls_state_t,decoder_t,ode_tag,hessian_t>;
  wls_system_t wlsSystem(appObj, yFOM_IC, yRef, decoderObj, numStepsInWindow, romSize, wlsStateIc);

  // create the wls state
  wls_state_t  wlsState(romSize*numStepsInWindow);
  ::pressio::ops::set_zero(wlsState);

  // -----------------
  // nonlinear solver
  // -----------------
  using gn_t = pressio::solvers::iterative::GaussNewton<linear_solver_t, wls_system_t>;
  gn_t GNSolver(wlsSystem, wlsState, linearSolver);
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
  fom_state_t yFinal(fomSize);
  using fom_state_reconstr_t = pressio::rom::FomStateReconstructor<scalar_t, fom_state_t, decoder_t>;
  fom_state_reconstr_t fomStateReconstructor(yRef, decoderObj);
  fomStateReconstructor(wlsCurrentState, yFinal);

  // get true solution
  const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(fomSize, dt, finalTime);

  for (int i=0;i<fomSize;i++){
    std::cout << std::setprecision(15) << yFinal[i] << " " << trueY[i] << "\n";
    if (std::abs(yFinal[i] - trueY[i]) > 1e-8) checkStr = "FAILED";
  }
  std::cout << checkStr << std::endl;

  return 0;
}
