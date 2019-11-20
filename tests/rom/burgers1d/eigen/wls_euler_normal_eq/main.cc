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

template<
  typename fom_type, typename wls_state_type,
  typename residual_pol_type, typename jac_pol_type
>
struct Specializer<pressio::rom::wls::WlsSystemDefaultApi<fom_type,wls_state_type,residual_pol_type,jac_pol_type>, fom_type, wls_state_type, residual_pol_type, jac_pol_type>
{ 
//  static_assert( ::pressio::rom::meta::is_legitimate_residual_policy_for_wls<residual_pol_t>,
//     "You are trying to use Wls with default api but the residual policy passed \
is not admissible for this: maybe you have the wrong api? blas blas");
  using type = pressio::rom::wls::WlsSystemDefaultApi<fom_type, wls_state_type, residual_pol_type, jac_pol_type>;
};
template<
  typename fom_type, typename wls_state_type,
  typename jtj_jtr_pol_type
>
struct Specializer<pressio::rom::wls::WlsSystemJTJApi<fom_type,wls_state_type,jtj_jtr_pol_type>, fom_type, wls_state_type, jtj_jtr_pol_type>
{ 
//  static_assert( ::pressio::rom::meta::is_legitimate_residual_policy_for_wls<residual_pol_t>,
//     "You are trying to use Wls with default api but the residual policy passed \
is not admissible for this: maybe you have the wrong api? blas blas");
  using type = pressio::rom::wls::WlsSystemJTJApi<fom_type, wls_state_type,jtj_jtr_pol_type>;
};

//struct DefaultApi{};
//struct JTJRApi{};
template<typename api_type,typename fom_type,typename wls_state_type,typename ... rest>
using WlsSystem = typename pressio::rom::wls::Specializer<api_type, fom_type, wls_state_type, rest...>::type;
}}} // end namespace rom::wls


int main(int argc, char *argv[]){
  using fom_t		= pressio::apps::Burgers1dEigenResidualApi;
  using scalar_t	= typename fom_t::scalar_type;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using rom_state_t	= pressio::containers::Vector<eig_dyn_vec>;
  using eig_dyn_mat	= Eigen::Matrix<scalar_t, -1, -1>;
  using decoder_jac_t	= pressio::containers::MultiVector<eig_dyn_mat>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t>;
  using wls_state_t     = std::vector<rom_state_t>;

  std::string checkStr {"PASSED"};

  //-------------------------------
  // app object
  int fomSize = 20;
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  fom_t appObj( mu, fomSize);
  auto t0 = static_cast<scalar_t>(0);
  scalar_t dt = 0.01;

  // read from file the jacobian of the decoder

  int romSize = 11;
  // store modes computed before from file
  decoder_jac_t phi =
    pressio::rom::test::eigen::readBasis("basis.txt", romSize, fomSize);
  int numBasis = phi.numVectors();
  if( numBasis != romSize ) return 0;

  // create decoder Obj
  decoder_t decoderObj(phi);

  constexpr std::size_t numStepsInWindow = 5;
  using stepper_stencil = ::pressio::ode::types::StepperTotalNumberOfStates<2>; 
  //using wls_problem      = pressio::rom::experimental:::WlsProblemGeneratorResidualApi<
//    DefaultWlsTypeGeneratorResidualApi, numStepsInWindow,
//    fom_t, wls_state_t, decoder_t, stepper_stencil, scalar_t>;

  // define standard rom state
  rom_state_t yROM(romSize);
  yROM.putScalar(2.0);

  // defined WLS state
  wls_state_t  wlsState(numStepsInWindow,rom_state_t(romSize) );

  using jtr_t = rom_state_t;
  using jtj_t = pressio::containers::MultiVector<eig_dyn_mat>;
  jtj_t jtj(romSize,romSize);
  jtr_t jtr(romSize);
  //using wls_t = pressio::rom::experimental::WlsSystem<fom_t,rom_state_t,JTJ_type,JTR_type,pressio::rom::wls::impl::JTJ_JTR_policy_smart>;
  using jtjr_policy_t = pressio::rom::wls::impl::JTJ_JTR_policy_smart<wls_state_t,fom_t,jtj_t,jtr_t>;
  using wls_api_t   = pressio::rom::wls::WlsSystemJTJApi<fom_t,wls_state_t,jtjr_policy_t>;
  using wls_system_t  = pressio::rom::wls::WlsSystem<wls_api_t, fom_t, wls_state_t,jtjr_policy_t>;
  // construct objects and test
  jtjr_policy_t jtjr_policy;
  wls_system_t wls(jtjr_policy);
  wls.JTJ_JTR(wlsState,jtj,jtr); 


  using resid_t = wls_state_t;
  using jacobian_t = pressio::containers::MultiVector<eig_dyn_mat>;
  resid_t residual(romSize);
  jacobian_t jacobian(romSize,romSize);
  using residual_policy_t = pressio::rom::wls::impl::residual_policy_naive<wls_state_t,fom_t,resid_t>;
  using jacobian_policy_t = pressio::rom::wls::impl::jacobian_policy_naive<wls_state_t,fom_t,jacobian_t>;
  using wls_naive_api_t   = pressio::rom::wls::WlsSystemDefaultApi<fom_t,wls_state_t,residual_policy_t,jacobian_policy_t>;
  using wls_naive_system_t  = pressio::rom::wls::WlsSystem<wls_naive_api_t, fom_t, wls_state_t,residual_policy_t,jacobian_policy_t>;
  //construct objects and test
  residual_policy_t residual_policy;
  jacobian_policy_t jacobian_policy;
  wls_naive_system_t wls_naive(residual_policy,jacobian_policy);
  wls_naive.residual(wlsState,residual); 




  //==================================
  using fom_native_state_t = pressio::apps::Burgers1dEigenResidualApi::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using fom_native_residual_t = pressio::apps::Burgers1dEigenResidualApi::residual_type; 
  using fom_residual_t             = ::pressio::containers::Vector<fom_native_residual_t>;

  fom_state_t yFOM(fomSize);

  fom_state_t yFOM2(fomSize);
  fom_residual_t fom_resid(fomSize);


  //This works:
  using fom_state_reconstr_t    = pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  const fom_state_reconstr_t fom_state_reconstr(yFOM,decoderObj);




  int step = 4;
  scalar_t time = 1.0;
  fom_state_reconstr(yROM,yFOM);

  appObj.timeDiscreteResidual(step,time,dt,*fom_resid.data(),*yFOM.data(),*yFOM2.data()); //nested dependence templated  
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
//  fom_state_t y1_const(fomSize);
//  y1_const.putScalar(1.0);
//  const fom_state_t *y1_const_ptr;
//  y1_const_ptr = &y1_const;
//  const rom_state_t yROM_const(romSize);
//  wls_container.reconstructFomStateAt(yROM,1);
  //std::array<double,4> test{4.};
  //cout << test[0] << endl;
//  std::cout << wls_container.fomStates_[0][19]  << std::endl;
  //.reconstructFomStateAt(yROM_const,1); 
//  wls_state_t wlsState(romSize,numStepsInWindow);
//  std::cout << wlsState(0,0) << std::endl;  
  //using fom_container_t = pressio::rom::FomStatesContainer<fom_state_t,numStepsInWindow,fom_state_reconstr_t>;
  //fom_container_t fom_container(fom_state_reconstr); 
  //using step_t        = typename fom_t::step_type;



  //fom_container.reconstructFomState(yROM);  
//  wls_system wls(residual_policy_t,jacobian_policy);

  /*
  //pressio::rom::experimental::WlsSystem<lspg_state_t,fom_t> wls(appobj); 
  // initialize to zero (this has to be done)
  yROM.putScalar(0.0);

  // define LSPG type
  //constexpr auto ode_case  = pressio::ode::ImplicitEnum::Euler;
  //using lspg_problem =  pressio::rom::LSPGUnsteadyProblem<pressio::rom::DefaultLSPGUnsteady, ode_case, fom_t, lspg_state_t, decoder_t>;
  using lspg_problem_t =  pressio::rom::LSPGUnsteadyProblem<pressio::rom::DefaultLSPGUnsteady,pressio::ode::ImplicitEnum::Euler, fom_t, lspg_state_t, decoder_t>;
  //using wls_problem_t = wlsProblemGenerator<fom_t, lspg_state_t, decoder_t>;
  using lspg_stepper_t = typename lspg_problem_t::lspg_stepper_t;
  //using wls_stepper_t = typename wls_problem_t::wls_stepper_t; 



  lspg_problem_t lspgProblem(appobj, yRef, decoderObj, yROM, t0);

  using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
  using hessian_t  = pressio::containers::Matrix<eig_dyn_mat>;
  using solver_tag   = pressio::solvers::linear::iterative::LSCG;
  using linear_solver_t = pressio::solvers::iterative::EigenIterative<solver_tag, hessian_t>;
  linear_solver_t linSolverObj;

  // GaussNewton solver
  // hessian comes up in GN solver, it is (J phi)^T (J phi)
  // rom is solved using eigen, hessian is wrapper of eigen matrix
  using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
  using gnsolver_t   =     pressio::solvers::iterative::GaussNewton<lspg_stepper_t, linear_solver_t>;

  gnsolver_t     solver   (lspgProblem.getStepperRef(), yROM, linSolverObj);

  //wlsProblem.residual(yROM);
  //wls_gnsolver_t solver2  (wlsProblem                  ,yROM, linSolverObj);
  solver.setTolerance(1e-13);
  // I know this should converge in few iters every step
  solver.setMaxIterations(2);

  // integrate in time
  pressio::ode::integrateNSteps(lspgProblem.getStepperRef(), yROM, 0.0, dt, 10, solver);

  // compute the fom corresponding to our rom final state
  auto yFomFinal = lspgProblem.getFomStateReconstructorCRef()(yROM);

  // this is a reproducing ROM test, so the final reconstructed state
  // has to match the FOM solution obtained with euler, same time-step, for 10 steps
  // const auto trueY = pressio::apps::test::Burg1DtrueImpEulerN20t010;
  const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(fomSize, dt, 0.10);
  for (auto i=0; i<yFomFinal.size(); i++){
    if (std::abs(yFomFinal[i] - trueY[i]) > 1e-10)
      checkStr = "FAILED";
  }

  std::cout << std::setprecision(14) << *yFomFinal.data() << std::endl;

  std::cout << checkStr <<  std::endl;
  */
  return 0;
}
