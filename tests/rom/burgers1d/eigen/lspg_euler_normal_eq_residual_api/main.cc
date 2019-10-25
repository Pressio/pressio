
#include "UTILS_ALL"
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG_UNSTEADY"
#include "APPS_UNSTEADYBURGERS1D"
#include "utils_eigen.hpp"


struct EulerLSPGWithResidualApi
{
  using fom_t		= pressio::apps::Burgers1dEigenResidualApi;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;
  using native_dmat_t   = typename fom_t::dense_matrix_type;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using decoder_jac_t	= pressio::containers::MultiVector<native_dmat_t>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t>;

  native_state_t fomSol_ = {};
  lspg_state_t yROM_ = {};

  EulerLSPGWithResidualApi()
  {
    std::string checkStr {"PASSED"};

    // app object
    constexpr int numCell = 20;
    Eigen::Vector3d mu(5.0, 0.02, 0.02);
    fom_t appobj( mu, numCell);
    auto t0 = static_cast<scalar_t>(0);
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
    yROM_.resize(romSize);
    yROM_.putScalar(0.0);

    // define LSPG type
    constexpr auto ode_case = pressio::ode::ImplicitEnum::Arbitrary;

    using stepper_order    = ::pressio::ode::types::StepperOrder<1>;
    using stepper_n_states = ::pressio::ode::types::StepperTotalNumberOfStates<2>;

    using lspg_problem	 = pressio::rom::LSPGUnsteadyProblem<
      pressio::rom::DefaultLSPGUnsteady, ode_case, fom_t,
      lspg_state_t, decoder_t, stepper_order, stepper_n_states, scalar_t>;
    using lspg_stepper_t	 = typename lspg_problem::lspg_stepper_t;
    lspg_problem lspgProblem(appobj, yRef, decoderObj, yROM_, t0);

    // linear solver
    using eig_dyn_mat	 = Eigen::Matrix<scalar_t, -1, -1>;
    using hessian_t	 = pressio::containers::Matrix<eig_dyn_mat>;
    using solver_tag	 = pressio::solvers::linear::iterative::LSCG;
    using linear_solver_t  = pressio::solvers::iterative::EigenIterative<solver_tag, hessian_t>;
    linear_solver_t linSolverObj;

    // GaussNewton solver
    // hessian comes up in GN solver, it is (J phi)^T (J phi)
    using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
    using gnsolver_t   = pressio::solvers::iterative::GaussNewton<lspg_stepper_t, linear_solver_t>;
    gnsolver_t solver(lspgProblem.getStepperRef(), yROM_, linSolverObj);
    solver.setTolerance(1e-13);
    solver.setMaxIterations(4);

    // integrate in time
    pressio::ode::integrateNSteps(lspgProblem.getStepperRef(), yROM_, 0.0, dt, 10, solver);

    // compute the fom corresponding to our rom final state
    auto yFomFinal = lspgProblem.getFomStateReconstructorCRef()(yROM_);
    fomSol_ = *yFomFinal.data();

    // this is a reproducing ROM test, so the final reconstructed state
    // has to match the FOM solution obtained with euler, same time-step, for 10 steps
    // const auto trueY = pressio::apps::test::Burg1DtrueImpEulerN20t010;
    const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(numCell, dt, 0.10);
    for (auto i=0; i<yFomFinal.size(); i++){
      if (std::abs(yFomFinal[i] - trueY[i]) > 1e-10)
        checkStr = "FAILED";
    }
    std::cout << checkStr <<  std::endl;
  }
};


struct EulerLSPGWithVelocityApi
{
  using fom_t		= pressio::apps::Burgers1dEigen;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;
  using native_dmat_t   = typename fom_t::dense_matrix_type;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using decoder_jac_t	= pressio::containers::MultiVector<native_dmat_t>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t>;

  native_state_t fomSol_ = {};
  lspg_state_t yROM_ = {};

  EulerLSPGWithVelocityApi()
  {
    std::string checkStr {"PASSED"};

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
    if( numBasis != romSize ) throw std::runtime_error("numBasis != romSize");

    // create decoder obj
    decoder_t decoderObj(phi);

    // for this problem, my reference state = initial state
    native_state_t yRef(numCell);
    yRef.setConstant(1);

    // define ROM state
    yROM_.resize(romSize);
    yROM_.putScalar(0.0);

    // define LSPG type
    constexpr auto ode_case = pressio::ode::ImplicitEnum::Euler;

    using lspg_problem	 = pressio::rom::LSPGUnsteadyProblem<
      pressio::rom::DefaultLSPGUnsteady, ode_case, fom_t, lspg_state_t, decoder_t>;
    using lspg_stepper_t	 = typename lspg_problem::lspg_stepper_t;
    lspg_problem lspgProblem(appobj, yRef, decoderObj, yROM_, t0);

    // linear solver
    using eig_dyn_mat	 = Eigen::Matrix<scalar_t, -1, -1>;
    using hessian_t	 = pressio::containers::Matrix<eig_dyn_mat>;
    using solver_tag	 = pressio::solvers::linear::iterative::LSCG;
    using linear_solver_t  = pressio::solvers::iterative::EigenIterative<solver_tag, hessian_t>;
    linear_solver_t linSolverObj;

    // GaussNewton solver
    // hessian comes up in GN solver, it is (J phi)^T (J phi)
    using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
    using gnsolver_t   = pressio::solvers::iterative::GaussNewton<lspg_stepper_t, linear_solver_t>;
    gnsolver_t solver(lspgProblem.getStepperRef(), yROM_, linSolverObj);
    solver.setTolerance(1e-13);
    solver.setMaxIterations(4);

    // integrate in time
    pressio::ode::integrateNSteps(lspgProblem.getStepperRef(), yROM_, 0.0, dt, 10, solver);

    // compute the fom corresponding to our rom final state
    auto yFomFinal = lspgProblem.getFomStateReconstructorCRef()(yROM_);
    fomSol_ = *yFomFinal.data();

    // this is a reproducing ROM test, so the final reconstructed state
    // has to match the FOM solution obtained with euler, same time-step, for 10 steps
    // const auto trueY = pressio::apps::test::Burg1DtrueImpEulerN20t010;
    const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(numCell, dt, 0.10);
    for (auto i=0; i<yFomFinal.size(); i++){
      if (std::abs(yFomFinal[i] - trueY[i]) > 1e-10)
        checkStr = "FAILED";
    }
    std::cout << checkStr <<  std::endl;
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
  for (auto i=0; i<veloRomSol.size(); i++){
    std::cout << std::setprecision(14)
  	      << veloRomSol[i]
  	      << " "
  	      << residRomSol[i]
  	      << std::endl;

    if (std::abs(veloRomSol[i] - residRomSol[i]) > 1e-13)
      checkStr = "FAILED";
  }

  std::cout << "check that fom reconstructed state match" << std::endl;
  // check the reconstructed fom state
  for (auto i=0; i<veloFomSol.size(); i++){
    std::cout << std::setprecision(14)
  	      << veloFomSol[i]
  	      << " "
  	      << residFomSol[i]
  	      << std::endl;

    if (std::abs(veloFomSol[i] - residFomSol[i]) > 1e-13)
      checkStr = "FAILED";
  }

  return 0;
}
