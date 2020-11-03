
#include "pressio_rom_lspg.hpp"
#include "pressio_apps.hpp"
#include "utils_eigen.hpp"
#include "custom_mapping.hpp"

int main(int argc, char *argv[])
{
  std::string checkStr {"PASSED"};

  using fom_t		= pressio::apps::Burgers1dEigenDiscreteTimeApi;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;
  using native_dmat_t   = typename fom_t::dense_matrix_type;
  using fom_state_t     = pressio::containers::Vector<native_state_t>;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using decoder_t	= MyCustomDecoder<native_dmat_t, fom_state_t>;

  // --- app object ---
  constexpr int numCell = 20;
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  fom_t appobj( mu, numCell);

  scalar_t dt = 0.01;
  constexpr int romSize = 11;

  // create decoder obj
  decoder_t decoderObj(romSize, numCell);

  // for this problem, my reference state = initial state
  native_state_t yRef(numCell);
  yRef.setConstant(1);

  // define ROM state
  lspg_state_t yROM_ = {};
  ::pressio::ops::resize(yROM_, romSize);
  ::pressio::ops::fill(yROM_, 0.0);

  // // define LSPG type
  // using ode_tag		 = pressio::ode::implicitmethods::Arbitrary;
  // using stepper_order    = ::pressio::ode::types::StepperOrder<1>;
  // using stepper_n_states = ::pressio::ode::types::StepperTotalNumberOfStates<2>;
  // using lspg_problem = typename pressio::rom::lspg::composeDefaultProblem<ode_tag, fom_t, 
  //       decoder_t, lspg_state_t, stepper_order, stepper_n_states>::type;
  // lspg_problem lspgProblem(appobj, decoderObj, yROM_, yRef);
  auto lspgProblem = pressio::rom::lspg::createDefaultProblemUnsteady<1,2>(
    appobj, decoderObj, yROM_, yRef);

  // linear solver
  using eig_dyn_mat	 = Eigen::Matrix<scalar_t, -1, -1>;
  using hessian_t	 = pressio::containers::DenseMatrix<eig_dyn_mat>;
  using solver_tag	 = pressio::solvers::linear::iterative::LSCG;
  using linear_solver_t  = pressio::solvers::linear::Solver<solver_tag, hessian_t>;
  linear_solver_t linSolverObj;

  // GaussNewton solver with normal equations
  auto solver = pressio::solvers::nonlinear::createGaussNewton(
    lspgProblem.stepperRef(), yROM_, linSolverObj);
  solver.setTolerance(1e-13);
  solver.setMaxIterations(4);

  // integrate in time
  pressio::rom::lspg::solveNSequentialMinimizations(lspgProblem, yROM_, 0.0, dt, 10, solver);

  // compute the fom corresponding to our rom final state
  auto yFomFinal = lspgProblem.fomStateReconstructorCRef()(yROM_);

  // this is a reproducing ROM test, so the final reconstructed state
  // has to match the FOM solution obtained with euler, same time-step, for 10 steps
  // const auto trueY = pressio::apps::test::Burg1DtrueImpEulerN20t010;
  const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(numCell, dt, 0.10);
  for (auto i=0; i<yFomFinal.extent(0); i++){
    if (std::abs(yFomFinal(i) - trueY[i]) > 1e-10)
      checkStr = "FAILED";
  }
  std::cout << checkStr <<  std::endl;

  return 0;
}
