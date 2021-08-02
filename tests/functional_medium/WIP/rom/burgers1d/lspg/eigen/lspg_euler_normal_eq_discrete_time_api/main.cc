
#include "pressio_rom_lspg.hpp"
#include "pressio_apps.hpp"
#include "utils_eigen.hpp"


struct EulerLSPGWithResidualApi
{
  using fom_t		= pressio::apps::Burgers1dEigenDiscreteTimeApi;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;
  using fom_state_t  = pressio::containers::Vector<native_state_t>;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using decoder_jac_t	= pressio::containers::MultiVector<Eigen::MatrixXd>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

  native_state_t fomSol_ = {};
  lspg_state_t yROM_ = {};

  EulerLSPGWithResidualApi()
  {
    std::string checkStr {"PASSED"};

    // app object
    constexpr int numCell = 20;
    Eigen::Vector3d mu(5.0, 0.02, 0.02);
    fom_t appobj( mu, numCell);
    scalar_t dt = 0.01;

    // read from file the jacobian of the decoder
    constexpr int romSize = 11;
    // store modes computed before from file
    decoder_jac_t phi =
      pressio::rom::test::eigen::readBasis("basis.txt", romSize, numCell);
    const int numBasis = phi.numVectors();
    if( numBasis != romSize )
      throw std::runtime_error("numBasis != romSize");

    // create decoder obj
    decoder_t decoderObj(phi);

    // for this problem, my reference state = initial state
    native_state_t yRef(numCell);
    yRef.setConstant(1);

    // define ROM state
    ::pressio::ops::resize(yROM_, romSize);
    ::pressio::ops::fill(yROM_, 0.0);

    // // define LSPG type
    // using ode_tag = pressio::ode::implicitmethods::Arbitrary;
    // using stepper_order    = ::pressio::ode::StepperOrder<1>;
    // using stepper_n_states = ::pressio::ode::StepperTotalNumberOfStates<2>;
    // using lspg_problem = typename pressio::rom::lspg::composeDefaultProblem<ode_tag, fom_t,
    //     decoder_t, lspg_state_t, stepper_order, stepper_n_states>::type;
    // lspg_problem lspgProblem(appobj, decoderObj, yROM_, yRef);
    auto lspgProblem = pressio::rom::lspg::createDefaultProblemUnsteady<1,2>(appobj, decoderObj, yROM_, yRef);

    // linear solver
    using eig_dyn_mat	 = Eigen::Matrix<scalar_t, -1, -1>;
    using hessian_t	 = pressio::containers::DenseMatrix<eig_dyn_mat>;
    using solver_tag	 = pressio::solvers::linear::iterative::LSCG;
    using linear_solver_t  = pressio::solvers::linear::Solver<solver_tag, hessian_t>;
    linear_solver_t linSolverObj;

    // GaussNewton solver with normal equations
    auto solver = pressio::rom::lspg::createGaussNewtonSolver(lspgProblem, yROM_, linSolverObj);
    solver.setTolerance(1e-13);
    solver.setMaxIterations(4);

    // solve
    pressio::rom::lspg::solveNSequentialMinimizations(lspgProblem, yROM_, 0.0, dt, 10, solver);

    // compute the fom corresponding to our rom final state
    auto yFomFinal = lspgProblem.fomStateReconstructorCRef()(yROM_);
    fomSol_ = *yFomFinal.data();
  }
};


struct EulerLSPGWithVelocityApi
{
  using fom_t		= pressio::apps::Burgers1dEigen;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;
  using fom_state_t  = pressio::containers::Vector<native_state_t>;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using decoder_jac_t	= pressio::containers::MultiVector<Eigen::MatrixXd>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

  native_state_t fomSol_ = {};
  lspg_state_t yROM_ = {};

  EulerLSPGWithVelocityApi()
  {
    std::string checkStr {"PASSED"};

    // app object
    constexpr int numCell = 20;
    Eigen::Vector3d mu(5.0, 0.02, 0.02);
    fom_t appobj( mu, numCell);
    scalar_t dt = 0.01;

    // read from file the jacobian of the decoder
    constexpr int romSize = 11;
    // store modes computed before from file
    decoder_jac_t phi =
      pressio::rom::test::eigen::readBasis("basis.txt", romSize, numCell);
    const int numBasis = phi.numVectors();
    if( numBasis != romSize ) throw std::runtime_error("numBasis != romSize");

    // create decoder obj
    decoder_t decoderObj(phi);

    // for this problem, my reference state = initial state
    native_state_t yRef(numCell);
    yRef.setConstant(1);

    // define ROM state
    ::pressio::ops::resize(yROM_, romSize);
    ::pressio::ops::fill(yROM_, 0.0);

    // define LSPG type
    using ode_tag = pressio::ode::implicitmethods::BDF1;
    auto lspgProblem = pressio::rom::lspg::createDefaultProblemUnsteady<ode_tag>(
      appobj, decoderObj, yROM_, yRef);

    // linear solver
    using eig_dyn_mat	 = Eigen::Matrix<scalar_t, -1, -1>;
    using hessian_t	 = pressio::containers::DenseMatrix<eig_dyn_mat>;
    using solver_tag	 = pressio::solvers::linear::iterative::LSCG;
    using linear_solver_t  = pressio::solvers::linear::Solver<solver_tag, hessian_t>;
    linear_solver_t linSolverObj;

    // GaussNewton solver with normal equations
    auto solver = pressio::rom::lspg::createGaussNewtonSolver(lspgProblem, yROM_, linSolverObj);
    solver.setTolerance(1e-13);
    solver.setMaxIterations(4);

    // integrate in time
    pressio::rom::lspg::solveNSequentialMinimizations(lspgProblem, yROM_, 0.0, dt, 10, solver);

    // compute the fom corresponding to our rom final state
    auto yFomFinal = lspgProblem.fomStateReconstructorCRef()(yROM_);
    fomSol_ = *yFomFinal.data();
  }
};


int main(int argc, char *argv[]){

  std::string checkStr {"PASSED"};

  EulerLSPGWithVelocityApi LSPGVeloApi;
  const auto veloFomSol = LSPGVeloApi.fomSol_;
  const auto veloRomSol = LSPGVeloApi.yROM_;

  EulerLSPGWithResidualApi LSPGResidApi;
  const auto residFomSol = LSPGResidApi.fomSol_;
  const auto residRomSol = LSPGResidApi.yROM_;

  std::cout << "check that gen coords match" << std::endl;
  // check the reconstructed rom state
  for (auto i=0; i<veloRomSol.extent(0); i++){
    std::cout << std::setprecision(14)
  	      << veloRomSol(i)
  	      << " "
  	      << residRomSol(i)
  	      << std::endl;

    if (std::abs(veloRomSol(i) - residRomSol(i)) > 1e-13)
      checkStr = "FAILED";
  }

  std::cout << "check that fom reconstructed state match" << std::endl;
  // check the reconstructed fom state
  for (auto i=0; i<veloFomSol.size(); i++){
    std::cout << std::setprecision(14)
  	      << veloFomSol(i)
  	      << " "
  	      << residFomSol[i]
  	      << std::endl;

    if (std::abs(veloFomSol(i) - residFomSol(i)) > 1e-13)
      checkStr = "FAILED";
  }

  std::cout << checkStr <<  std::endl;

  return 0;
}
