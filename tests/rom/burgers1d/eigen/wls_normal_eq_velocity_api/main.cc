#include <iostream>
#include <string>
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "SOLVERS_EXPERIMENTAL"
#include "ROM_LSPG_UNSTEADY"
#include "APPS_UNSTEADYBURGERS1D"
#include "utils_eigen.hpp"
#include "/Users/ejparis/pressio_repos/pressio/packages/rom/src/wls/wls_apis.hpp"
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
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;


  std::string checkStr {"PASSED"};

  //---- Information for WLS
  int fomSize = 20;
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  scalar_t dt = 0.01;
  int romSize = 11;
  constexpr int numStepsInWindow = 5;
  // ODE information

  using ode_tag = ::pressio::ode::implicitmethods::BDF2;

//  using ode_traits_t = pressio::ode::details::traits<ode_tag>;
//  ode_traits_t ode_traits;
  // update this to be inside
  const int t_stencil_width = 3;
  using aux_states_container_t = ::pressio::ode::AuxStatesContainer<false,fom_state_t,t_stencil_width>;
  //
  //-------------------------------
  // app object
  fom_t appObj( mu, fomSize);
  fom_state_t yFOM(fomSize);
  (*yFOM.data()).setConstant(1.);
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

  // Create WLS state
  wls_state_t  wlsState(romSize*numStepsInWindow);
  wlsState.setZero();



 // using gradient_t = wls_state_t; ( [Jphi]^T Jphi )
  using hessian_t = pressio::containers::Matrix<eig_dyn_mat>;


  using local_residual_policy_t = pressio::rom::wls::impl::local_residual_policy_velocityAPI;
  using local_jacobian_policy_t = pressio::rom::wls::impl::local_jacobian_policy_velocityAPI;
  /* Create policy to evaluate the hessian  and  gradient. Templates:
    fom_t: this is used to infer fom_state_t and  residual_t
    jac_t: this is used to infer the jacbobian type
    local_residual_policy_t: this is used to create a residual policy
    local_jacobian_policy_t: this is used to create a Jacobian policy 
    ode_tag: this is up in the air, but I'd like this to be used to infer the time stencil size. 
  */
  using hessian_gradient_policy_t = pressio::rom::wls::impl::hessian_gradient_policy<fom_t,decoder_t,local_residual_policy_t,local_jacobian_policy_t,ode_tag>;
  using wls_system_t   = pressio::rom::wls::WlsSystemHessianAndGradientApi<fom_t,wls_state_t,decoder_t,ode_tag,hessian_gradient_policy_t,aux_states_container_t,hessian_t>;
  wls_system_t wlsSystem(appObj,decoderObj,yFOM,numStepsInWindow,t_stencil_width);
// appObj, yREf, decoderObj, t0 
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

  const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF2::get(fomSize, dt, 0.10);
  const auto wlsCurrentState = pressio::containers::span(wlsState,(numStepsInWindow-1)*romSize,romSize);
  fom_state_t yFinal(fomSize);

  using fom_state_reconstr_t    = pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  fom_state_reconstr_t  fomStateReconstructor(yFOM,decoderObj);
  fomStateReconstructor(wlsCurrentState,yFinal);

  for (int i=0;i<fomSize;i++){
    cout << yFinal[i] << " " << trueY[i] << endl;
    if (std::abs(yFinal[i] - trueY[i]) > 1e-9) checkStr = "FAILED";
  }
  cout << checkStr << endl;

  return 0;
}
