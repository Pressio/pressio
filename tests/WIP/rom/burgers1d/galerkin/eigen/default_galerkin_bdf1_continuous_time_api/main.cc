
#include "pressio_rom_galerkin.hpp"
#include "pressio_apps.hpp"
#include "utils_eigen.hpp"

struct GalerkinBDF1WithContinuousTimeApi
{
  using fom_t		= pressio::apps::Burgers1dEigen;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;
  using fom_state_t	= pressio::containers::Vector<native_state_t>;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using rom_state_t	= pressio::containers::Vector<eig_dyn_vec>;
  using decoder_jac_t	= pressio::containers::MultiVector<Eigen::MatrixXd>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

  native_state_t fomSol_ = {};
  rom_state_t yROM_ = {};

  GalerkinBDF1WithContinuousTimeApi()
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

    using ode_tag = pressio::ode::BDF1;
    using pressio::rom::galerkin::createDefaultProblem;
    auto Problem = createDefaultProblem<ode_tag>(appobj, decoderObj, yROM_, yRef);
    using problem_t = decltype(Problem);

    // linear solver
    using solver_tag	 = pressio::solvers::linear::iterative::LSCG;
    using matrix_t = typename problem_t::galerkin_jacobian_t;
      using linear_solver_t  = pressio::solvers::linear::Solver<solver_tag, matrix_t>;
    linear_solver_t linSolverObj;

    // nonlinear system
    auto solver = pressio::rom::galerkin::create_newton_raphsonSolver(Problem, yROM_, linSolverObj);
    solver.setMaxIterations(2);

    // Integrate in time
    pressio::rom::galerkin::solveNSteps(Problem, yROM_, 0.0, dt, 10, solver);

    // compute the fom corresponding to our rom final state
    auto yFomFinal = Problem.fomStateReconstructorCRef()(yROM_);
    fomSol_ = *yFomFinal.data();
  }
};

int main(int argc, char *argv[])
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  std::string checkStr {"PASSED"};

  GalerkinBDF1WithContinuousTimeApi GalerkinRun;
  const auto FomSol = GalerkinRun.fomSol_;
  const auto RomSol = GalerkinRun.yROM_;

  std::cout << "check that fom reconstructed state match" << std::endl;

  const std::vector<double> gold = {
    1.2392461852464,
    1.0051322414122,
    1.0025874405045,
    1.0028353027323,
    1.0031333375861,
    1.0034628720682,
    1.0038270644436,
    1.0042295591478,
    1.0046743843334,
    1.0051659919295,
    1.0057093019206,
    1.0063097518487,
    1.0069733511236,
    1.0077067410977,
    1.0085172613415,
    1.0094130238705,
    1.0104029932025,
    1.0114970770758,
    1.0127062247348,
    1.0140425372763};

  if ((std::size_t)FomSol.size() != gold.size())
    checkStr = "FAILED";

  // check the reconstructed fom state
  for (auto i=0; i<FomSol.size(); i++){
    std::cout << std::setprecision(14)
	      << gold[i]
	      << " "
	      << FomSol[i]
	      << std::endl;

    if (std::abs(gold[i] - FomSol[i]) > 1e-11){
      checkStr = "FAILED";
      //break;
    }
  }

  std::cout << checkStr <<  std::endl;
  pressio::log::finalize();
  return 0;
}
