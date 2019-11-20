
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG_UNSTEADY"
#include "APPS_UNSTEADYBURGERS1D"
#include "utils_eigen.hpp"

#include "/Users/ejparis/pressio_repos/pressio/packages/solvers/src/meta/solvers_basic_meta.hpp"
#include "/Users/ejparis/pressio_repos/pressio/packages/solvers/src/meta/solvers_system_has_all_needed_jacobian_methods.hpp"
#include "/Users/ejparis/pressio_repos/pressio/packages/solvers/src/meta/solvers_system_has_all_needed_residual_methods.hpp"
#include "/Users/ejparis/pressio_repos/pressio/packages/rom/src/experimental/wls/rom_wls_sequential_residual_policy_for_residual_api.hpp"
//#include "/Users/ejparis/pressio_repos/pressio/packages/rom/src/experimental/wls/rom_wls_draft.hpp"
template<typename fom_t, typename lspg_state_t, typename decoder_t>
struct wlsProblemGenerator{
  //auto ode_case  = pressio::ode::ImplicitEnum::Euler;
  using lspg_problem_t = pressio::rom::LSPGUnsteadyProblem<pressio::rom::DefaultLSPGUnsteady,pressio::ode::ImplicitEnum::Euler, fom_t, lspg_state_t, decoder_t>;

  //using lspg_problem_t = pressio::rom::LSPGUnsteadyProblem<pressio::rom::DefaultLSPGUnsteady, ode_case                        , fom_t, lspg_state_t, decoder_t>;
  using scalar_t                = typename lspg_problem_t::scalar_t;
  using fom_native_state_t      = typename lspg_problem_t::fom_native_state_t;
  using lspg_stepper_t          = typename lspg_problem_t::lspg_stepper_t;
  using lspg_matrix_t           = typename lspg_problem_t::lspg_matrix_t;
  using lspg_residual_t   =  lspg_state_t;//fom_velocity_t;
  using wls_stepper_t = typename lspg_problem_t::lspg_stepper_t;


  using state_type      = lspg_state_t; 
  using residual_type   =  lspg_state_t; //this will change
  using jacobian_type   = lspg_matrix_t;
  using scalar_type = scalar_t;

  residual_type resid;
  jacobian_type Jphi;
  lspg_problem_t lspgProblem;
  wlsProblemGenerator(const fom_t & appobj,const fom_native_state_t &  yRef,decoder_t  & decoderObj,lspg_state_t  & yROM,scalar_t t0,int fomSize, int romSize) : lspgProblem(appobj, yRef, decoderObj, yROM, t0) , resid(fomSize) , Jphi(fomSize,romSize){};
  //lspg_stepper_t lspgStepper;
  void residual(const state_type & yROM, residual_type & resid) const{
    //lspgProblem.getStepperRef().residual(yROM,resid);
    //static_cast<const lspg_stepper_t>(lspgProblem.getStepperRef()).residual(yROM,resid);    
    //static_cast<const lspg_problem_t>(lspgProblem).getStepperRef().residual(yROM,resid);    
  };

  residual_type residual(const state_type & yROM) const{
     //return lspgProblem.getStepperRef().residual(yROM,resid);
     return static_cast<const lspg_stepper_t>(lspgProblem.getStepperRef()).residual(yROM,resid);  
  };

  void jacobian(const state_type & yROM, jacobian_type & Jphi) {
     //lspgProblem.getStepperRef().jacobian(yROM,Jphi);
  };
  jacobian_type jacobian(const state_type & yROM) const{
    return Jphi;//lspgProblem.getStepperRef().jacobian(yROM);
  };


};

int main(int argc, char *argv[]){
  using fom_t		= pressio::apps::Burgers1dEigen;
  //using fom_resid_t	= pressio::apps::Burgers1dEigenResidualApi::residual_type;
  //using fom_resid_t2	= fom_t::velocity_type;

  using scalar_t	= typename fom_t::scalar_type;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using eig_dyn_mat	= Eigen::Matrix<scalar_t, -1, -1>;
  using decoder_jac_t	= pressio::containers::MultiVector<eig_dyn_mat>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t>;

  std::string checkStr {"PASSED"};

  //-------------------------------
  // app object
  int fomSize = 20;
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  fom_t appobj( mu, fomSize);
  appobj.setup();
  auto t0 = static_cast<scalar_t>(0);
  scalar_t dt = 0.01;

  // read from file the jacobian of the decoder
  int romSize = 11;
  // store modes computed before from file
  decoder_jac_t phi =
    pressio::rom::test::eigen::readBasis("basis.txt", romSize, fomSize);
  int numBasis = phi.numVectors();
  if( numBasis != romSize ) return 0;

  // create decoder obj
  decoder_t decoderObj(phi);

  // for this problem, my reference state = initial state
  auto & yRef = appobj.getInitialState();

  // define ROM state
  lspg_state_t yROM(romSize);
  lspg_state_t resid(fomSize);

  // initialize to zero (this has to be done)
  yROM.putScalar(0.0);

  // define LSPG type
  //constexpr auto ode_case  = pressio::ode::ImplicitEnum::Euler;
  //using lspg_problem =  pressio::rom::LSPGUnsteadyProblem<pressio::rom::DefaultLSPGUnsteady, ode_case, fom_t, lspg_state_t, decoder_t>;
  using lspg_problem_t =  pressio::rom::LSPGUnsteadyProblem<pressio::rom::DefaultLSPGUnsteady,pressio::ode::ImplicitEnum::Euler, fom_t, lspg_state_t, decoder_t>;
  using wls_problem_t = wlsProblemGenerator<fom_t, lspg_state_t, decoder_t>;
  using lspg_stepper_t = typename lspg_problem_t::lspg_stepper_t;
  using wls_stepper_t = typename wls_problem_t::wls_stepper_t; 



  lspg_problem_t lspgProblem(appobj, yRef, decoderObj, yROM, t0);
  wls_problem_t wlsProblem(appobj, yRef, decoderObj, yROM, t0,fomSize, romSize);

  wlsProblem.residual(yROM,wlsProblem.resid);
//  wlsProblem.jacobian(yROM,wlsProblem.Jphi);

  //wls.lspgProblem;
  // linear solver
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
//  using wls_gnsolver_t   = pressio::solvers::iterative::GaussNewton< wls_problem_t, linear_solver_t>;

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
  return 0;
}
