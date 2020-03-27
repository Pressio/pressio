
#include "pressio_rom.hpp"
#include "pressio_apps.hpp"
#include "utils_eigen.hpp"
#include "rom_data_type_eigen.hpp"


#ifdef PRESSIO_ENABLE_TPL_TRILINOS
  #include "utils_tpetra.hpp"
#endif

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
  #include "utils_kokkos.hpp"
  #include "rom_data_type_kokkos.hpp"
#endif


namespace {


template <typename scalar_t>
auto readSol(::pressio::ode::implicitmethods::Euler & odeTag, const int fomSize, const scalar_t dt) 
  -> decltype(pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(fomSize, dt, 0.10)) 
{
  auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(fomSize, dt, 0.10);
  return trueY;
}

template <typename scalar_t>
auto readSol(::pressio::ode::implicitmethods::BDF2 & odeTag, const int fomSize, const scalar_t dt) 
  -> decltype(pressio::apps::test::Burgers1dImpGoldStatesBDF2::get(fomSize, dt, 0.10)) 
{
  auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF2::get(fomSize, dt, 0.10);
  return trueY;
}

//======= For eigen =============
template <typename decoder_d_t>
decoder_d_t readBasis( pressio::apps::Burgers1dEigen & appObj, ::pressio::ode::implicitmethods::Euler & odeTag, int  romSize, int  fomSize){

  const auto phiNative = pressio::rom::test::eigen::readBasis("basis_euler.txt", romSize, fomSize);
  decoder_d_t decoderObj(phiNative);

  return decoderObj;
}

template <typename decoder_d_t>
decoder_d_t readBasis( pressio::apps::Burgers1dEigen & appObj, ::pressio::ode::implicitmethods::BDF2 & odeTag, int  romSize, int  fomSize){

  const auto phiNative = pressio::rom::test::eigen::readBasis("basis_bdf2.txt", romSize, fomSize);
  decoder_d_t decoderObj(phiNative);

  return decoderObj;
}

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
//======= For tpetra =============
template <typename decoder_d_t, typename fom_dmat_t, typename rcpcomm_t>
decoder_d_t readBasis( pressio::apps::Burgers1dTpetra & appObj, ::pressio::ode::implicitmethods::Euler & odeTag, int  romSize, int  fomSize, rcpcomm_t Comm){

  const auto phiNative = pressio::rom::test::tpetra::readBasis("basis_euler.txt", romSize, fomSize,Comm, appObj.getDataMap());
  decoder_d_t decoderObj(phiNative);

  return decoderObj;
}

template <typename decoder_d_t, typename fom_dmat_t, typename rcpcomm_t>
decoder_d_t readBasis( pressio::apps::Burgers1dTpetra & appObj, ::pressio::ode::implicitmethods::BDF2 & odeTag, int  romSize, int  fomSize, rcpcomm_t Comm){

  const auto phiNative = pressio::rom::test::tpetra::readBasis("basis_bdf2.txt", romSize, fomSize,Comm, appObj.getDataMap());
  decoder_d_t decoderObj(phiNative);

  return decoderObj;
}

template<typename y1_t, typename y2_t>
std::string checkSol(pressio::apps::Burgers1dTpetra & appObj ,y1_t yFinal,y2_t trueY,int rank){
  std::string checkStr {"PASSED"};

  auto yFF_v = yFinal.data()->getData();
  int shift = (rank==0) ? 0 : 10;
  const int myn = yFinal.data()->getMap()->getNodeNumElements();
  for (auto i=0; i<myn; i++)
    if (std::abs(yFF_v[i] - trueY[i+shift]) > 1e-10) checkStr = "FAILED";

  return checkStr; 
}
//===================================


//for tpetra_block=================================
template <typename decoder_d_t, typename fom_dmat_t, typename rcpcomm_t>
decoder_d_t readBasis( pressio::apps::Burgers1dTpetraBlock & appObj, ::pressio::ode::implicitmethods::Euler & odeTag, int  romSize, int  fomSize, rcpcomm_t Comm){
  auto tpw_phi = pressio::rom::test::tpetra::readBasis("basis_euler.txt", romSize,fomSize, Comm, appObj.getDataMap());
  fom_dmat_t tpb_phi(*tpw_phi.data(), *appObj.getDataMap(), 1);
  decoder_d_t decoderObj(tpb_phi);
  return decoderObj;
}

template <typename decoder_d_t, typename fom_dmat_t, typename rcpcomm_t>
decoder_d_t readBasis( pressio::apps::Burgers1dTpetraBlock & appObj, ::pressio::ode::implicitmethods::BDF2 & odeTag, int  romSize, int  fomSize, rcpcomm_t Comm){
  auto tpw_phi = pressio::rom::test::tpetra::readBasis("basis_bdf2.txt", romSize,fomSize, Comm, appObj.getDataMap());
  fom_dmat_t tpb_phi(*tpw_phi.data(), *appObj.getDataMap(), 1);
  decoder_d_t decoderObj(tpb_phi);
  return decoderObj;
}

template<typename y1_t, typename y2_t>
std::string checkSol(pressio::apps::Burgers1dTpetraBlock & appObj ,y1_t yFinal,y2_t trueY,int rank){

  std::string checkStr {"PASSED"};
  auto yFF_v = yFinal.data()->getVectorView().getData();
  int shift = (rank==0) ? 0 : 10;
  const int myn = yFinal.data()->getMap()->getNodeNumElements();
  for (auto i=0; i<myn; i++)
    if (std::abs(yFF_v[i] - trueY[i+shift]) > 1e-10) checkStr = "FAILED";

  return checkStr;
}
//==============================================
 
/// MPI Version
template <
  typename fom_t,
  typename rom_data_t,
  typename tcomm_t,
  typename hessian_matrix_structure_tag,
  typename ode_tag,
  typename rcpcomm_t
  >
std::string doRun(rcpcomm_t & Comm, int rank)
{

  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t= typename fom_t::state_type;
  using fom_dmat_t         = typename fom_t::dense_matrix_type;
  using fom_state_t   = ::pressio::containers::Vector<native_state_t>;
  using decoder_jac_d_t	= pressio::containers::MultiVector<fom_dmat_t>;

  // wls state type
  using wls_state_d_t	= typename rom_data_t::wls_state_d_t;
  using hessian_d_t	= typename rom_data_t::hessian_d_t;

  // decoder jacobian type
  using decoder_d_t	= pressio::rom::LinearDecoder<decoder_jac_d_t, wls_state_d_t, fom_state_t>;

  constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();

  // app object
  constexpr int fomSize = 20;
  fom_t appObj({{5.0, 0.02, 0.02}}, fomSize,Comm);
  constexpr scalar_t dt = 0.01;
  constexpr auto t0 = zero;

  int romSize = 11;

  // create/read jacobian of the decoder
  ode_tag odeTag;
  auto decoderObj = readBasis<decoder_d_t,fom_dmat_t>(appObj,odeTag,romSize,fomSize,Comm);


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
  using linear_solver_t = typename rom_data_t::linear_solver_t;
  linear_solver_t linear_solver;

  // -----------------
  // WLS problem
  // -----------------
  constexpr int numStepsInWindow = 5;
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
  const auto trueY = readSol(odeTag,fomSize, dt);
  //const   auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(fomSize, dt, 0.10);
  std::string checkStr = checkSol(appObj ,yFinal,trueY,rank);
  return checkStr;
}

#endif

//Scalar version
template <
  typename fom_t,
  typename rom_data_t,
  typename hessian_matrix_structure_tag,
  typename ode_tag
  >
std::string doRun()
{
  std::string checkStr {"PASSED"};
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t= typename fom_t::state_type;
  using fom_dmat_t         = typename fom_t::dense_matrix_type;
  using fom_state_t   = ::pressio::containers::Vector<native_state_t>;
  using decoder_jac_d_t	= pressio::containers::MultiVector<fom_dmat_t>;

  // wls state type
  using wls_state_d_t	= typename rom_data_t::wls_state_d_t;
  using hessian_d_t	= typename rom_data_t::hessian_d_t;

  // decoder jacobian type
  using decoder_d_t	= pressio::rom::LinearDecoder<decoder_jac_d_t, wls_state_d_t, fom_state_t>;

  constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();

  // app object
  constexpr int fomSize = 20;
  fom_t appObj( Eigen::Vector3d{5.0, 0.02, 0.02}, fomSize);
  constexpr scalar_t dt = 0.01;
  constexpr auto t0 = zero;

  int romSize = 11;

  // create/read jacobian of the decoder
  ode_tag odeTag;
  auto decoderObj = readBasis<decoder_d_t>(appObj,odeTag,romSize,fomSize);

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
  using linear_solver_t = typename rom_data_t::linear_solver_t;
  linear_solver_t linear_solver;

  // -----------------
  // WLS problem
  // -----------------
  constexpr int numStepsInWindow = 5;
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
  const auto trueY = readSol(odeTag,fomSize, dt);

  for (int i=0;i<fomSize;i++){
    std::cout << std::setprecision(15) << yFinal[i] << " " << trueY[i] << "\n";
    if (std::abs(yFinal[i] - trueY[i]) > 1e-8) checkStr = "FAILED";
  }
  return checkStr;
}



}// end namespace



