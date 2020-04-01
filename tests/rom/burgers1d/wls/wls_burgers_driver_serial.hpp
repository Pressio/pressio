
#ifndef PRESSIO_TESTS_WLS_BURGERS1D_DRIVER_SERIAL_HPP_
#define PRESSIO_TESTS_WLS_BURGERS1D_DRIVER_SERIAL_HPP_

#include "pressio_rom.hpp"
#include "pressio_apps.hpp"
#include "rom_data_type_eigen.hpp"
#include "burgers_fom_functions_eigen.hpp"

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "burgers_fom_functions_kokkos.hpp"
#include "rom_data_type_kokkos.hpp"
#endif

namespace
{

template <typename scalar_t>
auto readSol(::pressio::ode::implicitmethods::Euler odeTag, const std::size_t fomSize, const scalar_t dt)
  -> decltype(pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(fomSize, dt, 0.10))
{
  auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(fomSize, dt, 0.10);
  return trueY;
}

template <typename scalar_t>
auto readSol(::pressio::ode::implicitmethods::BDF2 odeTag, const std::size_t fomSize, const scalar_t dt)
  -> decltype(pressio::apps::test::Burgers1dImpGoldStatesBDF2::get(fomSize, dt, 0.10))
{
  auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF2::get(fomSize, dt, 0.10);
  return trueY;
}


template <typename fom_t>
fom_t appConstructor(std::size_t fomSize){};

template <>
  pressio::apps::Burgers1dEigen appConstructor<pressio::apps::Burgers1dEigen>(std::size_t fomSize)
{
  pressio::apps::Burgers1dEigen appObj( Eigen::Vector3d{5.0, 0.02, 0.02}, fomSize);
  return appObj;
}

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <>
  pressio::apps::Burgers1dKokkos appConstructor<pressio::apps::Burgers1dKokkos>(std::size_t fomSize)
{
  pressio::apps::Burgers1dKokkos appObj({{5.0, 0.02, 0.02}}, fomSize);
  return appObj;
}
#endif

} //end anonymous namespace


namespace pressio{ namespace testing{ namespace wls{

template <
  typename fom_t,
  typename rom_data_t,
  typename hessian_matrix_structure_tag,
  typename ode_tag
  >
std::string doRun()
{
  // types
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;
  using fom_dmat_t      = typename fom_t::dense_matrix_type;
  using fom_state_t     = pressio::containers::Vector<native_state_t>;
  using decoder_jac_d_t	= pressio::containers::MultiVector<fom_dmat_t>;
  using wls_state_d_t	= typename rom_data_t::wls_state_d_t;
  using wls_hessian_d_t	= typename rom_data_t::wls_hessian_d_t;
  using decoder_d_t	= pressio::rom::LinearDecoder<decoder_jac_d_t, wls_state_d_t, fom_state_t>;

  // app object
  constexpr std::size_t fomSize = 20;
  auto appObj = appConstructor<fom_t>(fomSize);
  constexpr scalar_t dt = 0.01;
  constexpr auto t0 = ::pressio::utils::constants::zero<scalar_t>();
  // for this problem, my reference state = initial state
  auto & yFOM_IC_native = appObj.getInitialState();
  // wrap into pressio container
  fom_state_t yFOM_IC(yFOM_IC_native);
  //reference state is equal to the IC
  fom_state_t & yRef = yFOM_IC;


  constexpr pressio::rom::wls::rom_size_t romSize = 11;
  constexpr pressio::rom::wls::window_size_t numStepsInWindow = 5;
  constexpr pressio::rom::wls::rom_size_t wlsSize = romSize*numStepsInWindow;
  constexpr scalar_t finalTime = 0.1;
  constexpr pressio::rom::wls::window_size_t numWindows = (finalTime/dt)/numStepsInWindow;

  //  jacobian of the decoder
  auto decoderObj = readBasis<decoder_d_t>(appObj,ode_tag(),romSize,fomSize);

  //  lin solver
  using linear_solver_t = typename rom_data_t::linear_solver_t;
  linear_solver_t linear_solver;

  //  WLS problem
  using wls_system_t = pressio::rom::wls::SystemHessianAndGradientApi<fom_t, wls_state_d_t, decoder_d_t,
								      ode_tag, wls_hessian_d_t,
								      hessian_matrix_structure_tag>;
  // create the wls system
  wls_system_t wlsSystem(appObj, yFOM_IC, yRef, decoderObj, numStepsInWindow, romSize, linear_solver);

  // create the wls state
  wls_state_d_t  wlsState(wlsSize);
  pressio::ops::set_zero(wlsState);

  // - NL solver
  using gn_t = pressio::solvers::iterative::GaussNewton<linear_solver_t, wls_system_t>;
  gn_t GNSolver(wlsSystem, wlsState, linear_solver);
  GNSolver.setTolerance(1e-13);
  GNSolver.setMaxIterations(5);

  //  solve wls problem
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
  const auto trueY = readSol(ode_tag(),fomSize, dt);

  std::string checkStr = checkSol(appObj ,yFinal,trueY,fomSize);
  return checkStr;
}

}}} //end namespace pressio::testing::wls
#endif
