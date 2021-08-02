
#include "pressio_rom_galerkin.hpp"
#include "pressio_apps.hpp"
#include "utils_eigen.hpp"

using fom_t		= pressio::apps::Burgers1dEigen;
using scalar_t		= typename fom_t::scalar_type;
using native_state_t	= typename fom_t::state_type;
using native_velo_t	= typename fom_t::velocity_type;
using fom_state_t	= pressio::containers::Vector<native_state_t>;
using fom_velo_t	= pressio::containers::Vector<native_velo_t>;
using decoder_jac_t	= pressio::containers::MultiVector<Eigen::MatrixXd>;
using decoder_t		= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;
using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
using rom_state_t	= pressio::containers::Vector<eig_dyn_vec>;

struct MyCollocator
{
  std::vector<int> rows_;
  MyCollocator(std::initializer_list<int> l) : rows_(l){}

  Eigen::MatrixXd sampleRows(const Eigen::MatrixXd & operand)
  {
    Eigen::MatrixXd result(rows_.size(), operand.cols());
    for (std::size_t i=0; i<(std::size_t)rows_.size(); ++i){
      for (std::size_t j=0; j<(std::size_t)operand.cols(); ++j){
	result(i,j) = operand(rows_[i], j);
      }
    }
    return result;
  }
};

struct Masker
{
  std::vector<int> rows_;
  Masker(std::initializer_list<int> l) : rows_(l){}

  native_velo_t createApplyMaskResult(const native_velo_t & src) const{
    return native_velo_t(rows_.size());
  }

  Eigen::MatrixXd createApplyMaskResult(const Eigen::MatrixXd & src) const{
    return Eigen::MatrixXd(rows_.size(), src.cols());
  }

  void applyMask(const native_velo_t & src, double /*t*/, native_velo_t & dest) const
  {
    assert((std::size_t)dest.rows() == (std::size_t)rows_.size());
    for (std::size_t i=0; i<rows_.size(); ++i){
      dest(i) = src(rows_[i]);
    }
  }

  void applyMask(const Eigen::MatrixXd & src, double /*t*/, Eigen::MatrixXd & dest) const
  {
    assert((std::size_t)dest.rows() == (std::size_t)rows_.size());
    for (std::size_t i=0; i<(std::size_t)rows_.size(); ++i){
      for (std::size_t j=0; j<(std::size_t)src.cols(); ++j){
	dest(i,j) = src(rows_[i], j);
      }
    }
  }
};

struct GalerkinRunnerContinuousTimeApi
{
  native_state_t fomSol_ = {};
  rom_state_t yROM_ = {};

  GalerkinRunnerContinuousTimeApi()
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

    // for this problem, my reference state = initial state which is all ones
    native_state_t yRef(numCell); yRef.setConstant(1);

    // define ROM state
    ::pressio::ops::resize(yROM_, romSize);
    ::pressio::ops::fill(yROM_, 0.0);

    // create decoder obj
    decoder_t decoderObj(phi);

    // set indices for masking
    // only pick subset of points (randomly, this is just a test)
    auto indices = {0,1,4,5,7,10,11,12,15,17,18,19};

    // since this is a masked Galerkin problem, we need the collocator to
    // create a projector suitable for the masked operator
    MyCollocator mapper(indices);
    auto projector = pressio::rom::galerkin::createCollocationProjector(decoderObj, std::move(mapper));

    // create masker
    Masker myMask(indices);

    using ode_tag = pressio::ode::implicitmethods::BDF1;
    using pressio::rom::galerkin::createMaskedVelocityProblem;
    auto Problem = createMaskedVelocityProblem<ode_tag>(appobj, decoderObj, yROM_, yRef, myMask, projector);

    // linear solver
    using solver_tag	   = pressio::solvers::linear::iterative::LSCG;
    using matrix_t	   = typename decltype(Problem)::galerkin_jacobian_t;
    using linear_solver_t  = pressio::solvers::linear::Solver<solver_tag, matrix_t>;
    linear_solver_t linSolverObj;

    // nonlinear system
    auto solver = pressio::rom::galerkin::createNewtonRaphsonSolver(Problem, yROM_, linSolverObj);
    solver.setMaxIterations(2);

    // Integrate in time
    pressio::rom::galerkin::solveNSteps(Problem, yROM_, 0.0, dt, 10, solver);

    // compute the fom corresponding to our rom final state
    // note that the fom is reconstructed using the decoder, so this is the fom state
    // on the FULL mesh
    auto yFomFinal = Problem.fomStateReconstructorCRef()(yROM_);
    fomSol_ = *yFomFinal.data();
  }
};

int main(int argc, char *argv[])
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  std::string checkStr {"PASSED"};

  GalerkinRunnerContinuousTimeApi GalerkinRun;
  const auto FomSol = GalerkinRun.fomSol_;
  const auto RomSol = GalerkinRun.yROM_;

  // note that the gold solution is one the full mesh
  const std::vector<double> gold = {
    1.2392461852463,
    1.0051322405528,
    1.0000008791913,
    0.99999102530353,
    1.0032108975228,
    1.0018038224142,
    0.99834268966878,
    1.004357342716,
    0.99913801977844,
    0.99761665316683,
    1.0026034173479,
    1.0081599176557,
    1.0053676941224,
    1.0029342607457,
    1.0039685412351,
    1.0076667218855,
    1.0040395527636,
    1.0113197844164,
    1.0108250707108,
    1.0131324724384};

  if ((std::size_t)FomSol.size() != gold.size())
    checkStr = "FAILED";

  // check the reconstructed fom state
  for (auto i=0; i<FomSol.size(); i++){
    std::cout << std::setprecision(14)
	      << gold[i]
	      << " "
	      << FomSol[i]
	      << std::endl;

    if (std::abs(gold[i] - FomSol[i]) > 1e-13){
      checkStr = "FAILED";
      break;
    }
  }

  std::cout << checkStr <<  std::endl;
  pressio::log::finalize();
  return 0;
}
