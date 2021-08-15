
#include "pressio_rom_lspg.hpp"
#include "pressio_apps.hpp"
#include "utils_eigen.hpp"

int main(int argc, char *argv[])
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::trace});

  std::string checkStr {"PASSED"};

  using fom_t		= pressio::apps::Burgers1dEigen;
  using scalar_t	= typename fom_t::scalar_type;

  // -------------------------------------------------------
  // create FOM object
  // -------------------------------------------------------
  constexpr int numCell = 20;
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  fom_t appobj( mu, numCell);

  // for this problem, my reference state = initial state
  auto & yRef = appobj.getInitialState();

  // -------------------------------------------------------
  // read basis
  // -------------------------------------------------------
  using decoder_jac_t	= pressio::containers::MultiVector<Eigen::MatrixXd>;
  constexpr int romSize = 11;
  decoder_jac_t phi = pressio::rom::test::eigen::readBasis("basis.txt", romSize, numCell);
  const int numBasis = phi.numVectors();
  if( numBasis != romSize ) return 0;

  // -------------------------------------------------------
  // create decoder obj
  // -------------------------------------------------------
  using native_state_t  = typename fom_t::state_type;
  using fom_state_t  = pressio::containers::Vector<native_state_t>;
  using decoder_t = pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;
  decoder_t decoderObj(phi);

  // -------------------------------------------------------
  // create ROM problem
  // -------------------------------------------------------
  using lspg_state_t = pressio::containers::Vector<Eigen::Matrix<scalar_t,-1,1>>;

  // define ROM state
  lspg_state_t yROM(romSize);
  // initialize to zero (this has to be done)
  pressio::ops::fill(yROM, 0.0);

  // // define LSPG type
  using ode_tag  = pressio::ode::implicitmethods::BDF1;
  auto lspgProblem = pressio::rom::lspg::createDefaultProblemUnsteady<ode_tag>(
    appobj, decoderObj, yROM, yRef);

  // linear solver
  using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
  using hessian_t  = pressio::containers::DenseMatrix<eig_dyn_mat>;
  using solver_tag   = pressio::solvers::linear::iterative::LSCG;
  using linear_solver_t = pressio::solvers::linear::Solver<solver_tag, hessian_t>;
  linear_solver_t linSolverObj;

  // GaussNewton solver with normal equations
  auto solver = pressio::rom::lspg::create_gauss_newtonSolver(lspgProblem, yROM, linSolverObj);
  solver.setTolerance(1e-13);
  // I know this should converge in few iters every step
  solver.setMaxIterations(2);

  // solve
  scalar_t dt = 0.01;
  pressio::rom::lspg::solveNSequentialMinimizations(lspgProblem, yROM, 0.0, dt, 10, solver);

  // compute the fom corresponding to our rom final state
  auto yFomFinal = lspgProblem.fomStateReconstructorCRef()(yROM);

  // this is a reproducing ROM test, so the final reconstructed state
  // has to match the FOM solution obtained with euler, same time-step, for 10 steps
  // const auto trueY = pressio::apps::test::Burg1DtrueImpEulerN20t010;
  const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(numCell, dt, 0.10);
  for (auto i=0; i<yFomFinal.extent(0); i++){
    if (std::abs(yFomFinal(i) - trueY[i]) > 1e-10)
      checkStr = "FAILED";
  }

  std::cout << std::setprecision(14) << *yFomFinal.data() << std::endl;
  std::cout << checkStr <<  std::endl;
  pressio::log::finalize();
  
  return 0;
}
