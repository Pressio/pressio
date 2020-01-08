#include <iostream>
#include <string>
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "SOLVERS_EXPERIMENTAL"
#include "ROM_BASIC"
#include "rom/src/meta/wls_velocity_api/rom_model_meets_velocity_api_for_unsteady_wls.hpp"
#include "APPS_UNSTEADYBURGERS1D"
#include "utils_eigen.hpp"
#include "/Users/ejparis/pressio_repos/pressio/packages/rom/src/wls/apis/wls_apis.hpp"
//#include "/Users/ejparis/pressio_repos/pressio/packages/rom/src/wls/policies/wls_default.hpp"

//#include "../wls_files/my_gauss_newton.hpp"
//#include "../wls_files/wls_specializers.hpp"

using namespace std;

int main(int argc, char *argv[]){
  using fom_t   = pressio::apps::Burgers1dEigen;
  using scalar_t  = typename fom_t::scalar_type;
  using fom_native_state_t      = typename fom_t::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using wls_state_t	= pressio::containers::Vector<eig_dyn_vec>;
  using eig_dyn_mat	= Eigen::Matrix<scalar_t, -1, -1>;
  using decoder_jac_t	= pressio::containers::MultiVector<eig_dyn_mat>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t>;
  using fom_state_t     = ::pressio::containers::Vector<fom_native_state_t>;
  using hessian_t       = pressio::containers::Matrix<eig_dyn_mat>;
  std::string checkStr {"PASSED"};

  //---- Information for WLS
  int fomSize = 20;
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  scalar_t dt = 0.01;
  int romSize = 11;
  constexpr int numStepsInWindow = 5;
  // ODE information

  using ode_tag = ::pressio::ode::implicitmethods::Euler;

  const int t_stencil_width = 2;
  using aux_states_container_t = ::pressio::ode::AuxStatesContainer<false,fom_state_t,t_stencil_width>;
  //
  //-------------------------------
  // app object
  fom_t appObj( mu, fomSize);
  auto & yRef_native = appObj.getInitialState();
  fom_state_t yRef(yRef_native);
  auto t0 = static_cast<scalar_t>(0);
  // read from file the jacobian of the decoder
  // store modes computed before from fil
  decoder_jac_t phi =
    pressio::rom::test::eigen::readBasis("basis.txt", romSize, fomSize);
  int numBasis = phi.numVectors();
  if( numBasis != romSize ) return 0;
  // create decoder Obj
  decoder_t decoderObj(phi);

  
  // Create WLS state
  wls_state_t  wlsState(romSize*numStepsInWindow);
  wlsState.setZero();





  using hessian_gradient_policy_t = pressio::rom::wls::hessian_gradient_policy<fom_t,decoder_t,ode_tag>;
  using wls_system_t   = pressio::rom::wls::WlsSystemHessianAndGradientApi<fom_t,wls_state_t,decoder_t,ode_tag,hessian_gradient_policy_t,aux_states_container_t,hessian_t>;
  wls_system_t wlsSystem(appObj,decoderObj,yRef,numStepsInWindow,t_stencil_width);



  // Use pressio interface for solvers 
  using solver_tag   = pressio::solvers::linear::direct::ColPivHouseholderQR;
  using linear_solver_t = pressio::solvers::direct::EigenDirect<solver_tag, hessian_t>;
  linear_solver_t linear_solver; 
  using gn_t = pressio::solvers::iterative::GaussNewton<linear_solver_t, wls_system_t>;
  gn_t GNSolver(wlsSystem, wlsState, linear_solver);

  int numSteps = 10/numStepsInWindow;
  for (int step = 0; step < numSteps; step++)
  {
    wlsSystem.advanceOneWindow(wlsState,GNSolver,step,dt);
  }

  const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(fomSize, dt, 0.10);
  const auto wlsCurrentState = pressio::containers::span(wlsState,(numStepsInWindow-1)*romSize,romSize);
  fom_state_t yFinal(fomSize);

  using fom_state_reconstr_t    = pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  fom_state_reconstr_t  fomStateReconstructor(yRef,decoderObj);
  fomStateReconstructor(wlsCurrentState,yFinal);

  for (int i=0;i<fomSize;i++){
    cout << yFinal[i] << " " << trueY[i] << endl;
    if (std::abs(yFinal[i] - trueY[i]) > 1e-9) checkStr = "FAILED";
  }
  cout << checkStr << endl;

  return 0;
}
