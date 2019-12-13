#include <iostream>
#include <string>
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG_UNSTEADY"
#include "APPS_UNSTEADYBURGERS1D"
#include "utils_eigen.hpp"
#include "../wls_files/wls_apis.hpp"
#include "../wls_files/my_gauss_newton.hpp"
#include "../wls_files/wls_specializers.hpp"

// n is the window size
using namespace std;


int main(int argc, char *argv[]){
  using fom_t   = pressio::apps::Burgers1dEigen;
  using scalar_t  = typename fom_t::scalar_type;
  using fom_native_state_t      = typename fom_t::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using rom_state_t	= pressio::containers::Vector<eig_dyn_vec>;
  using eig_dyn_mat	= Eigen::Matrix<scalar_t, -1, -1>;
  using decoder_jac_t	= pressio::containers::MultiVector<eig_dyn_mat>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t>;
  using fom_native_state_t = pressio::apps::Burgers1dEigenResidualApi::state_type;

  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using wls_state_t    = ::pressio::containers::MultiVector<Eigen::Matrix<scalar_t, -1, -1,Eigen::ColMajor>>;

  //using jac_t = decoder_t::jacobian_t;



  std::string checkStr {"PASSED"};

  //---- Information for WLS
  int fomSize = 20;
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  scalar_t dt = 0.01;
  int romSize = 11;
  constexpr int numStepsInWindow = 10;
  const int t_stencil_width = 3;

  using ode_tag = ::pressio::ode::implicitmethods::BDF2;
  using aux_states_container_t = ::pressio::ode::AuxStatesContainer<false,fom_state_t,t_stencil_width>;// statesContainer(yFOM); 

  //-------------------------------
  // app object
  fom_t appObj( mu, fomSize);
  auto t0 = static_cast<scalar_t>(0);
  // read from file the jacobian of the decoder
  // store modes computed before from file
  decoder_jac_t phi =
    pressio::rom::test::eigen::readBasis("basis.txt", romSize, fomSize);
  int numBasis = phi.numVectors();
  if( numBasis != romSize ) return 0;
  // create decoder Obj
  decoder_t decoderObj(phi);
  decoderObj.getReferenceToJacobian();
  wls_state_t  wlsState(romSize,numStepsInWindow);
  wlsState.setZero();


  fom_state_t yFOM(fomSize);
  fom_native_state_t yRef(fomSize);
  (*yFOM.data()).setConstant(1.);

  using fom_state_reconstr_t    = pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  fom_state_reconstr_t  fomStateReconstructor(yFOM,decoderObj);


  using gradient_t = rom_state_t;
  using hessian_t = pressio::containers::Matrix<eig_dyn_mat>;
  using local_residual_policy_t = pressio::rom::wls::impl::local_residual_policy_velocityAPI<fom_t>;
  using local_jacobian_policy_t = pressio::rom::wls::impl::local_jacobian_policy_velocityAPI<fom_t,decoder_t>;
  using hessian_gradient_policy_t = pressio::rom::wls::impl::hessian_gradient_policy<wls_state_t,fom_t,hessian_t,gradient_t,decoder_t,local_residual_policy_t,local_jacobian_policy_t,ode_tag>;
  using wls_api_t   = pressio::rom::wls::WlsSystemHessianAndGradientApi<fom_t,wls_state_t,decoder_t,hessian_gradient_policy_t,aux_states_container_t>;
  using wls_system_t  = pressio::rom::wls::WlsSystem<wls_api_t, fom_t, wls_state_t,decoder_t, hessian_gradient_policy_t,aux_states_container_t>;
  // construct objects and test
  hessian_gradient_policy_t hessian_gradient_policy(romSize,fomSize,numStepsInWindow,t_stencil_width,phi);
  wls_system_t wlsSystem(appObj,hessian_gradient_policy,decoderObj,yFOM,dt,numStepsInWindow,romSize,fomSize,t_stencil_width);

  double t = 0.;
  int numSteps = 10/numStepsInWindow;

  // Use pressio interface for linear solvers 
  using solver_tag   = pressio::solvers::linear::direct::ColPivHouseholderQR;
  using linear_solver_t = pressio::solvers::direct::EigenDirect<solver_tag, hessian_t>;
  linear_solver_t linear_solver; 
  using gn_type = my_gauss_newton_class<wls_system_t,wls_state_t,hessian_t,gradient_t,linear_solver_t>;
  gn_type gn_solver(wlsSystem,wlsState,linear_solver);

  using app_jacob_t = typename fom_t::jacobian_type;
  using ode_jac_t   = pressio::containers::Matrix<app_jacob_t>;
  for (int step = 0; step < numSteps; step++)
  {
    wlsSystem.advanceOneWindow(wlsState,gn_solver,step);
  }
  fom_state_t yTest(fomSize);
  aux_states_container_t  fomStateContainerObj(yTest);
  ::pressio::ode::impl::time_discrete_residual<ode_tag>(yTest,yTest,fomStateContainerObj,dt);


  const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF2::get(fomSize, dt, 0.10);
  const auto wlsCurrentState = wlsState.viewColumnVector(numStepsInWindow-1);
  fom_state_t yFinal(fomSize);

  fomStateReconstructor(wlsCurrentState,yFinal);

  for (int i=0;i<fomSize;i++){
    cout << yFinal[i] << " " << trueY[i] << endl;
    if (std::abs(yFinal[i] - trueY[i]) > 1e-9) checkStr = "FAILED";
  }
  cout << checkStr << endl;
//  ::pressio::ode::impl::time_discrete_residual<ode_case_t,fom_state_t,fom_state_t,fom_state_t,scalar_t>time_discrete_residual(yFOM,yFOM,yFOM,dt);
//  ::pressio::ode::impl::time_discrete_residual<ode_case_t2>(yFOM,yFOM,statesContainer,dt);

/*

*/
  return 0;
}
