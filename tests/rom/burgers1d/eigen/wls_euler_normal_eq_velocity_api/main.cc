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

  std::string checkStr {"PASSED"};

  //---- Information for WLS
  int fomSize = 20;
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  scalar_t dt = 0.01;
  int romSize = 11;
  constexpr int numStepsInWindow = 5;
  int t_stencil_width = 2;

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
  wls_state_t  wlsStateIC( romSize,t_stencil_width-1);
  wlsState.setZero();
  wlsStateIC.setZero();

  fom_state_t yFOM(fomSize);
  fom_native_state_t yRef(fomSize);
  (*yFOM.data()).setConstant(1.);
  using gradient_t = rom_state_t;
  using hessian_t = pressio::containers::Matrix<eig_dyn_mat>;
  hessian_t hessian(romSize*numStepsInWindow,romSize*numStepsInWindow);
  gradient_t gradient(romSize*numStepsInWindow);
  gradient_t grad(romSize*numStepsInWindow);

  using local_residual_policy_t = pressio::rom::wls::impl::local_residual_policy_velocityAPI<fom_t>;
  using local_jacobian_policy_t = pressio::rom::wls::impl::local_jacobian_policy_velocityAPI<fom_t,decoder_t>;
  using hessian_gradient_policy_t = pressio::rom::wls::impl::hessian_gradient_policy<wls_state_t,fom_t,hessian_t,gradient_t,decoder_t,local_residual_policy_t,local_jacobian_policy_t>;
  using wls_api_t   = pressio::rom::wls::WlsSystemHessianAndGradientApi<fom_t,wls_state_t,decoder_t,hessian_gradient_policy_t>;
  using wls_system_t  = pressio::rom::wls::WlsSystem<wls_api_t, fom_t, wls_state_t,decoder_t, hessian_gradient_policy_t>;

  // construct objects and test
  hessian_gradient_policy_t hessian_gradient_policy(romSize,fomSize,numStepsInWindow,2,phi);
  wls_system_t wlsSystem(appObj,hessian_gradient_policy,decoderObj,yFOM,dt,numStepsInWindow,romSize,fomSize,t_stencil_width);

  double t = 0.;
  int numSteps = 10/numStepsInWindow;

  // Use pressio interface for linear solvers 
  using solver_tag   = pressio::solvers::linear::direct::ColPivHouseholderQR;
  using linear_solver_t = pressio::solvers::direct::EigenDirect<solver_tag, hessian_t>;
  linear_solver_t linear_solver; 

  using gn_type = my_gauss_newton_class<wls_system_t,wls_state_t,hessian_t,gradient_t,linear_solver_t>;
  gn_type gn_solver(wlsSystem,wlsState,linear_solver);

  constexpr auto ode_case = pressio::ode::ImplicitEnum::Euler;
  using app_jacob_t = typename fom_t::jacobian_type;
  using ode_jac_t   = pressio::containers::Matrix<app_jacob_t>;
  using stepper_t = pressio::ode::ImplicitStepper<ode_case, fom_state_t, fom_state_t, ode_jac_t, fom_t>;
  stepper_t stepperObj(yFOM, appObj);
  using nm1 = pressio::ode::nMinusOne;
  auto & odeState_nm1 = stepperObj.auxStates_.template get<nm1>();
  cout << *(odeState_nm1).data() << endl;
  auto test = stepperObj.residual(yFOM);

  for (int step = 0; step < numSteps; step++)
  {
    gn_solver.my_gauss_newton(wlsSystem,wlsState,romSize,numStepsInWindow); 
    (*wlsSystem.wlsStateIC.data()).block(0,0,romSize,1) = (*wlsState.data()).block(0,numStepsInWindow-1,romSize,1); 
    cout << " step " << step << " completed " << endl;
  }

  const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(fomSize, dt, 0.1);

  using fom_state_reconstr_t    = pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  fom_state_reconstr_t  fomStateReconstructor(yFOM,decoderObj);
  fom_state_t yFinal(fomSize);
  const auto wlsCurrentState = wlsState.viewColumnVector(numStepsInWindow-1);
  fomStateReconstructor(wlsCurrentState,yFinal);

  for (int i=0;i<fomSize;i++){
    cout << yFinal[i] << " " << trueY[i] << endl;
    if (std::abs(yFinal[i] - trueY[i]) > 1e-7) checkStr = "FAILED";
  }
  cout << checkStr << endl;


//  using line_search_t = pressio::solvers::iterative::gn::noLineSearch;
//  using converged_t = pressio::solvers::iterative::converged_when::absoluteNormCorrectionBelowTol;

  using nlSolver_t =  pressio::solvers::iterative::GaussNewton<linear_solver_t, wls_system_t>;
//  nlSolver_t nlSolver(wlsSystem,wlsState,linear_solver);
/*
*/
  return 0;
}
