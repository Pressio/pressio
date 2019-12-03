#include <iostream>
#include <string>
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG_UNSTEADY"
#include "APPS_UNSTEADYBURGERS1D"
#include "utils_eigen.hpp"
#include "wls_apis.hpp"
#include "wls_containers.hpp"
//#include "/Users/ejparis/pressio_repos/pressio/packages/solvers/src/meta/solvers_basic_meta.hpp"
//#include "/Users/ejparis/pressio_repos/pressio/packages/solvers/src/meta/solvers_system_has_all_needed_jacobian_methods.hpp"
//#include "/Users/ejparis/pressio_repos/pressio/packages/solvers/src/meta/solvers_system_has_all_needed_residual_methods.hpp"
//#include "/Users/ejparis/pressio_repos/pressio/packages/rom/src/experimental/wls/rom_wls_sequential_residual_policy_for_residual_api.hpp"
//#include "/Users/ejparis/pressio_repos/pressio/packages/rom/src/experimental/wls/rom_wls_draft.hpp"

// n is the window size
using namespace std;
namespace pressio{ namespace rom{ namespace wls{


template<typename api_tag, typename fom_type, typename wls_state_type,typename ... rest>
struct Specializer{
  using type = void;
};

//struct DefaultApi{};

/*
template<
  typename fom_type, typename wls_state_type,
  typename residual_pol_type, typename jac_pol_type, typename decoder_type
>
struct Specializer<pressio::rom::wls::WlsSystemDefaultApi<fom_type,wls_state_type,residual_pol_type,jac_pol_type,decoder_type>, fom_type, wls_state_type, residual_pol_type, jac_pol_type,decoder_type>
{ 
//  static_assert( ::pressio::rom::meta::is_legitimate_residual_policy_for_wls<residual_pol_t>,
//     "You are trying to use Wls with default api but the residual policy passed \
is not admissible for this: maybe you have the wrong api? blas blas");
  using type = pressio::rom::wls::WlsSystemDefaultApi<fom_type, wls_state_type, residual_pol_type, jac_pol_type,decoder_type>;
};
*/

template<
  typename fom_type, typename wls_state_type,
  typename decoder_type, typename fom_state_container_type, typename hessian_gradient_pol_type
>
struct Specializer<pressio::rom::wls::WlsSystemJTJApi<fom_type,wls_state_type,decoder_type,fom_state_container_type, hessian_gradient_pol_type>, fom_type, wls_state_type, decoder_type, fom_state_container_type,hessian_gradient_pol_type>
{ 
  using type = pressio::rom::wls::WlsSystemJTJApi<fom_type, wls_state_type,decoder_type,fom_state_container_type, hessian_gradient_pol_type>;
};

//struct DefaultApi{};
//struct JTJRApi{};
template<typename api_type,typename fom_type,typename wls_state_type,typename ... rest>
using WlsSystem = typename pressio::rom::wls::Specializer<api_type, fom_type, wls_state_type, rest...>::type;
}}} // end namespace rom::wls



template <typename system_type, typename state_type, typename hessian_type, typename gradient_type, typename linear_solver_type>
class my_gauss_newton_class{
public:
  using scalar_type = typename system_type::scalar_type;
  hessian_type hessian_ = {};
  gradient_type gradient_ = {};
  gradient_type dx_ = {};
  linear_solver_type & linear_solver_ = {};

  my_gauss_newton_class(const system_type & system,const state_type & stateIn,linear_solver_type & linear_solver) : hessian_(system.createHessianObject(stateIn)) , gradient_(system.createGradientObject(stateIn)),dx_(system.createGradientObject(stateIn)), linear_solver_(linear_solver){}; 

  void my_gauss_newton(const system_type & sys, state_type & state, state_type & stateIC, int romSize, int numStepsInWindow)
  {
    double gnorm = 1.;
    double tol = 1e-8;
    scalar_type rnorm = 0.; 
    int iteration = 0;
    while (gnorm >= tol){
      sys.computeHessianAndGradient(state,stateIC,hessian_,gradient_); 
      linear_solver_.solve(hessian_,gradient_,dx_);
      iteration += 1;
      gnorm = (*gradient_.data()).squaredNorm();
      cout << "grad norm = " << gnorm << endl; 
      cout << "Hess norm = " << (*hessian_.data()).squaredNorm() << endl; 
      cout << "x    norm = " << (*dx_.data()).squaredNorm() << endl; 
      cout << "iteration = " << iteration << endl; 
      //cout << (*hessian.data()) << endl;
      hessian_.setZero();
      for (int n = 0; n < numStepsInWindow; n++){
        *(state[n]).data() += (*dx_.data()).block(n*romSize,0,romSize,1); 
        cout << "===============" << endl;
      }
    } 
  } 
};


int main(int argc, char *argv[]){
  using fom_t		= pressio::apps::Burgers1dEigenResidualApi;
  using scalar_t	= typename fom_t::scalar_type;
  using fom_native_state_t      = typename fom_t::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using rom_state_t	= pressio::containers::Vector<eig_dyn_vec>;
  using eig_dyn_mat	= Eigen::Matrix<scalar_t, -1, -1>;
  using decoder_jac_t	= pressio::containers::MultiVector<eig_dyn_mat>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t>;
  using fom_native_state_t = pressio::apps::Burgers1dEigenResidualApi::state_type;
  using fom_native_jacobian_t = pressio::apps::Burgers1dEigenResidualApi::jacobian_type;

  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using fom_native_residual_t = pressio::apps::Burgers1dEigenResidualApi::residual_type; 
  using fom_residual_t             = ::pressio::containers::Vector<fom_native_residual_t>;

  // Windowed least-squares types. std::vector will be replaced in the future with something custom; e.g., multivector
  using wls_state_t     = std::vector<rom_state_t>; 
//  using wls_state_t     = pressio::MultiVector<Eigen::Matrix<scalar_t, -1, -1,Eigen::ColMajor>>;


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
  //using stepper_stencil = ::pressio::ode::types::StepperTotalNumberOfStates<2>;
  //using fom_state_container_t = typename pressio::ode::AuxStatesContainer<fom_state_t, 2>; 
  using fom_state_container_t = std::vector<fom_state_t>; 

  // define WLS state
  wls_state_t  wlsState(numStepsInWindow,rom_state_t(romSize) );
  wls_state_t  wlsStateIC(t_stencil_width-1,rom_state_t(romSize) );

  // create initial conditions for WLS
  //wls_state_t  wlsStateIC(t_stencil_width-1,rom_state_t(romSize) );
  for (int i=0;i<numStepsInWindow; i++){wlsState[i].putScalar(0.0);}
  for (int i=0;i<numStepsInWindow; i++){wlsStateIC[i].putScalar(0.0);}

  fom_state_t yFOM(fomSize);
  fom_native_state_t yRef(fomSize);
  (*yFOM.data()).setConstant(1.);
  //yRef.setConstant(1);
  //auto & yFOM = appObj.getInitialState();


  using gradient_t = rom_state_t;
  using hessian_t = pressio::containers::Matrix<eig_dyn_mat>;
  hessian_t hessian(romSize*numStepsInWindow,romSize*numStepsInWindow);
  gradient_t gradient(romSize*numStepsInWindow);
  gradient_t grad(romSize*numStepsInWindow);

  using hessian_gradient_policy_t = pressio::rom::wls::impl::JTJ_JTR_policy_smart<wls_state_t,fom_t,hessian_t,gradient_t,decoder_t>;
  using wls_api_t   = pressio::rom::wls::WlsSystemJTJApi<fom_t,wls_state_t,decoder_t,fom_state_container_t,hessian_gradient_policy_t>;
  using wls_system_t  = pressio::rom::wls::WlsSystem<wls_api_t, fom_t, wls_state_t,decoder_t, fom_state_container_t, hessian_gradient_policy_t>;



  // construct objects and test
  hessian_gradient_policy_t hessian_gradient_policy(2,phi,romSize,numStepsInWindow,fomSize);

  wls_system_t wlsSystem(appObj,hessian_gradient_policy,decoderObj,yFOM,dt,numStepsInWindow,romSize,fomSize);

  double t = 0.;
  int numSteps = 10/numStepsInWindow;

  // Use pressio interface for linear solvers 
  using solver_tag   = pressio::solvers::linear::direct::ColPivHouseholderQR;
  using linear_solver_t = pressio::solvers::direct::EigenDirect<solver_tag, hessian_t>;
  linear_solver_t linear_solver; 

  using gn_type = my_gauss_newton_class<wls_system_t,wls_state_t,hessian_t,gradient_t,linear_solver_t>;
  gn_type gn_solver(wlsSystem,wlsState,linear_solver);

  for (int step = 0; step < numSteps; step++)
  {
    //my_gauss_newton<wls_system_t,wls_state_t,hessian_t,gradient_t,linear_solver_t>(wls,wlsState,wlsStateIC,hessian,gradient,grad,romSize,numStepsInWindow,linear_solver); 
    gn_solver.my_gauss_newton(wlsSystem,wlsState,wlsStateIC,romSize,numStepsInWindow); 

    *wlsStateIC[0].data() = *wlsState[numStepsInWindow-1].data();
    cout << " step " << step << " completed " << endl;
//    cout << "n1 " << (*wlsStateIC[0].data()).squaredNorm() << endl;
//    *wlsState[numStepsInWindow-1].data() = *wlsState[numStepsInWindow-1].data()*2; 
//    cout << "n2 " << (*wlsStateIC[0].data()).squaredNorm() << endl;

  }

  const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(fomSize, dt, 0.1);

  using fom_state_reconstr_t    = pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  fom_state_reconstr_t  fomStateReconstructor(yFOM,decoderObj);
  fom_state_t yFinal(fomSize);
  fomStateReconstructor(wlsState[numStepsInWindow - 1],yFinal);

  for (int i=0;i<fomSize;i++){
    cout << yFinal[i] << " " << trueY[i] << endl;
    if (std::abs(yFinal[i] - trueY[i]) > 1e-7) checkStr = "FAILED";
  }
  cout << checkStr << endl;


  using line_search_t = pressio::solvers::iterative::gn::noLineSearch;
  using converged_t = pressio::solvers::iterative::converged_when::absoluteNormCorrectionBelowTol;

  using t1 = wls_system_t::hessian_type;
  using t2 = wls_system_t::gradient_type;
  using t3 = wls_system_t::scalar_type;
/*
*/
//template <
//  typename system_type,
//  typename linear_solver_type,
//  typename scalar_type,
//  typename line_search_type,
//  typename converged_when
//  
//  using nlSolver_t =  pressio::solvers::iterative::GaussNewton<linear_solver_t, wls_system_t>;
  //pressio::solvers::iterative::impl::experimental::GaussNewtonHessianGradientApi<wls_system_t,linear_solver_t,scalar_t,line_search_t,converged_t>;
  //nlSolver_t nlSolver(wlsSystem,wlsState,linear_solver);
  //wls.JTJ_JTR(wlsState,wlsState,hessian,gradient);
  //cout << (*hessian.data()).size() << endl;
  //cout << (*gradient.data()).size() << endl;
  //auto x = (*hessian.data()).colPivHouseholderQr().solve((*gradient.data())); 

  /*
  //(*hessian.data()).block(1,1,testInt,testInt) = (*hessian2.data()); 
//  using resid_t = wls_state_t;
//  using jacobian_t = pressio::containers::MultiVector<eig_dyn_mat>;
//  resid_t residual(romSize);
//  jacobian_t jacobian(romSize,romSize);
//  using residual_policy_t = pressio::rom::wls::impl::residual_policy_naive<wls_state_t,fom_t,resid_t,decoder_t>;
//  using jacobian_policy_t = pressio::rom::wls::impl::jacobian_policy_naive<wls_state_t,fom_t,jacobian_t>;
//  using wls_naive_api_t   = pressio::rom::wls::WlsSystemDefaultApi<fom_t,wls_state_t,residual_policy_t,jacobian_policy_t,decoder_t>;
//  using wls_naive_system_t  = pressio::rom::wls::WlsSystem<wls_naive_api_t, fom_t, wls_state_t,residual_policy_t,jacobian_policy_t,decoder_t>;
//  //construct objects and test
//  residual_policy_t residual_policy;
//  jacobian_policy_t jacobian_policy;
//  wls_naive_system_t wls_naive(appObj,residual_policy,jacobian_policy,decoderObj,yFOM,dt);
  //==================================


  //This works:
  using fom_state_reconstr_t    = pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  const fom_state_reconstr_t fom_state_reconstr(yFOM,decoderObj);




  int step = 4;
  scalar_t time = 1.0;
  appObj.timeDiscreteResidual(step,time,dt,*fom_resid.data(),*yFOM.data(),*yFOM2.data()); //nested dependence templated  
  
  decoder_jac_t phi2(fomSize,romSize);

  // This doesn't work: 
  //appObj.applyTimeDiscreteJacobian(step,time,dt,*decoderObj.getReferenceToJacobian(),0,*phi2.data(),*yFOM.data(),*yFOM2.data());

  // This does:
  // 
  // computes A = J B
  const auto & phiC = decoderObj.getReferenceToJacobian();
  appObj.applyTimeDiscreteJacobian(step,time,dt,*phiC.data(),0,*phi2.data(),*yFOM.data(),*yFOM2.data());

  //appObj.applyTimeDiscreteJacobian(step,time,dt,*phi.data(),0,*phi2.data(),*yFOM.data(),*yFOM2.data());

  cout << fom_resid(4) << endl;

  using wls_container_t = pressio::rom::wls::WlsFomStatesContainer<numStepsInWindow,fom_state_t,fom_state_reconstr_t>;
  ///using fom_container_t = pressio::rom::FomStatesStaticContainer<fom_state_t, 2, fom_state_reconstr_t>;
  using fom_container_t = pressio::ode::StatesContainer<fom_state_t, 2>; // just use this fow now

  fom_container_t fom_container(yFOM);
  cout << fom_container[1][0] << endl;
  fom_state_reconstr(wlsState[0],fom_container[1]);
  //fom_container.reconstructCurrentFomState(wlsState[0]); 
  for (int i=1; i< numStepsInWindow; i++){
    fom_state_reconstr(wlsState[i],fom_container[0]);
    fom_state_reconstr(wlsState[i-1],fom_container[1]);
    appObj.timeDiscreteResidual(i,dt*i,dt,*fom_resid.data(),*fom_container[0].data(),*fom_container[1].data());
  }
//  appObj.template timeDiscreteResidual(step,time,dt,fom_resid,y1.data,y2.data); //nested dependence templated  
  using wls_container_t = pressio::rom::wls::WlsFomStatesContainer<numStepsInWindow,fom_state_t,fom_state_reconstr_t>;

  */
  return 0;
}
