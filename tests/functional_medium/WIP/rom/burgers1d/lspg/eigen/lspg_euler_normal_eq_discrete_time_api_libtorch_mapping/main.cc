
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
  constexpr int numCell = 100;
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  fom_t appobj( mu, numCell);

  auto t0 = static_cast<scalar_t>(0);
  scalar_t fint = 1.00;
  scalar_t dt = 0.01;
  auto Nsteps = static_cast<::pressio::ode::step_type>(fint/dt);
  constexpr int romSize = 8;

  // create decoder obj
  decoder_t decoderObj(romSize, numCell, "basis.pt");

  // for this problem, my reference state = initial state
  native_state_t yRef(numCell);
  yRef.setConstant(1);

  // define ROM state
  lspg_state_t yROM_ = {};
  ::pressio::ops::resize(yROM_, romSize);
  yROM_(0) = -0.0072;
  yROM_(1) = -0.7842;
  yROM_(2) = -0.6795;
  yROM_(3) = -0.6600;
  yROM_(4) = -0.7119;
  yROM_(5) = -0.6995;
  yROM_(6) = -0.7610;
  yROM_(7) = -0.3957;

  // // define LSPG type
  // using ode_tag		 = pressio::ode::implicitmethods::Arbitrary;
  // using stepper_order    = ::pressio::ode::StepperOrder<1>;
  // using stepper_n_states = ::pressio::ode::StepperTotalNumberOfStates<2>;
  //   using lspg_problem = typename pressio::rom::lspg::composeDefaultProblem<ode_tag, fom_t, 
  //     decoder_t, lspg_state_t, stepper_order, stepper_n_states>::type;
  // lspg_problem lspgProblem(appobj, decoderObj, yROM_, yRef);
  auto lspgProblem = pressio::rom::lspg::createDefaultProblemUnsteady<1,2>(
    appobj, decoderObj, yROM_, yRef);

  // linear solver
  using eig_dyn_mat	 = Eigen::Matrix<scalar_t, -1, -1>;
  using hessian_t	 = pressio::containers::DenseMatrix<eig_dyn_mat>;
  using solver_tag	 = pressio::solvers::linear::iterative::LSCG;
  using linear_solver_t  = pressio::solvers::linear::Solver<solver_tag, hessian_t>;
  linear_solver_t linSolverObj;

  // GaussNewton solver
  auto solver = pressio::rom::lspg::createGaussNewtonSolver(lspgProblem, yROM_, linSolverObj);
  solver.setTolerance(1e-13);
  solver.setMaxIterations(4);
  // solve
  pressio::rom::lspg::solveNSequentialMinimizations(lspgProblem, yROM_, 0.0, dt, Nsteps, solver);

  // compute the fom corresponding to our rom final state
  auto yFomFinal = lspgProblem.fomStateReconstructorCRef()(yROM_);
  
  // TODO: check accuracy of final state

  for (auto i=0; i<yFomFinal.extent(0); i++){
    std::cout << yFomFinal[i] << std::endl;
  }

  std::cout << checkStr <<  std::endl;

  return 0;
}
