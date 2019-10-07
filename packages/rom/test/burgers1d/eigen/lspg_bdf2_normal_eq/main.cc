
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "APPS_UNSTEADYBURGERS1D"
#include "utils_eigen.hpp"

int main(int argc, char *argv[]){
  using fom_t		= pressio::apps::Burgers1dEigen;
  using scalar_t	= typename fom_t::scalar_type;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using eig_dyn_mat	= Eigen::Matrix<scalar_t, -1, -1>;
  using decoder_jac_t	= pressio::containers::MultiVector<eig_dyn_mat>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t>;

  std::string checkStr {"PASSED"};

  //-------------------------------
  // app object
  constexpr int numCell = 20;
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  fom_t appobj( mu, numCell);
  appobj.setup();
  auto t0 = static_cast<scalar_t>(0);
  scalar_t dt = 0.01;

  // read from file the jacobian of the decoder
  constexpr int romSize = 11;
  // store modes computed before from file
  decoder_jac_t phi =
    pressio::rom::test::eigen::readBasis("basis.txt", romSize, numCell);
  const int numBasis = phi.numVectors();
  if( numBasis != romSize ) return 0;

  // create decoder obj
  decoder_t decoderObj(phi);

  // for this problem, my reference state = initial state
  auto & yRef = appobj.getInitialState();

  // define ROM state
  lspg_state_t yROM(romSize);
  // initialize to zero (this has to be done)
  yROM.putScalar(0.0);

  // define LSPG type
  constexpr auto ode_case  = pressio::ode::ImplicitEnum::BDF2;
  using lspg_problem_types = pressio::rom::DefaultLSPGTypeGenerator<
    fom_t, ode_case, decoder_t, lspg_state_t>;
  pressio::rom::LSPGUnsteadyProblemGenerator<lspg_problem_types> lspgProblem(
      appobj, yRef, decoderObj, yROM, t0);

  using lspg_stepper_t = typename lspg_problem_types::lspg_stepper_t;

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
  using gnsolver_t   = pressio::solvers::iterative::GaussNewton<
    lspg_stepper_t, linear_solver_t>;
  gnsolver_t solver(lspgProblem.getStepperRef(), yROM, linSolverObj);
  solver.setTolerance(1e-13);
  solver.setMaxIterations(200);

  // integrate in time
  pressio::ode::integrateNSteps(lspgProblem.getStepperRef(), yROM, 0.0, dt, 10, solver);

  // compute the fom corresponding to our rom final state
  auto yFomFinal = lspgProblem.getFomStateReconstructorCRef()(yROM);

  // this is a reproducing ROM test, so the final reconstructed state
  // has to match the FOM solution obtained with bdf2, same time-step, for 10 steps
  // const auto trueY = pressio::apps::test::Burg1DtrueImpBDF2N20t010;
  const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF2::get(numCell, dt, 0.10);  
  for (auto i=0; i<yFomFinal.size(); i++)
    if (std::abs(yFomFinal[i] - trueY[i]) > 1e-10) checkStr = "FAILED";

  std::cout << checkStr <<  std::endl;
  return 0;
}
