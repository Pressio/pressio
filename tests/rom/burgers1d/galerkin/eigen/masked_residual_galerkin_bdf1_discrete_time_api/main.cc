
#include "pressio_rom_galerkin.hpp"
#include "pressio_apps.hpp"
#include "utils_eigen.hpp"

using fom_t		= pressio::apps::Burgers1dEigenDiscreteTimeApi;
using scalar_t		= typename fom_t::scalar_type;
using native_state_t	= typename fom_t::state_type;
using native_residual_t = typename fom_t::discrete_time_residual_type;
using native_dmat_t	= typename fom_t::dense_matrix_type;

using fom_state_t	= pressio::containers::Vector<native_state_t>;
using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
using rom_state_t	= pressio::containers::Vector<eig_dyn_vec>;
using decoder_jac_t	= pressio::containers::MultiVector<native_dmat_t>;
using decoder_t		= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

struct MyCollocator
{
  std::vector<int> rows_;
  MyCollocator(std::initializer_list<int> l) : rows_(l){}

  native_dmat_t sampleRows(const native_dmat_t & operand)
  {
    native_dmat_t result(rows_.size(), operand.cols());
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

  native_residual_t createApplyMaskResult(const native_residual_t & src) const{
    return native_residual_t(rows_.size());
  }

  native_dmat_t createApplyMaskResult(const native_dmat_t & src) const{
    return native_dmat_t(rows_.size(), src.cols());
  }

  void applyMask(const native_residual_t & src, double /*t*/, native_residual_t & dest) const
  {
    assert((std::size_t)dest.rows() == (std::size_t)rows_.size());
    for (std::size_t i=0; i<rows_.size(); ++i){
      dest(i) = src(rows_[i]);
    }
  }

  void applyMask(const native_dmat_t & src, double /*t*/, native_dmat_t & dest) const
  {
    assert((std::size_t)dest.rows() == (std::size_t)rows_.size());
    for (std::size_t i=0; i<(std::size_t)rows_.size(); ++i){
      for (std::size_t j=0; j<(std::size_t)src.cols(); ++j){
	dest(i,j) = src(rows_[i], j);
      }
    }
  }
};

struct GalerkinBDF1WithResidualApi
{
  native_state_t fomSol_ = {};
  rom_state_t yROM_ = {};

  GalerkinBDF1WithResidualApi()
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

    // for this problem, my reference state = initial state
    native_state_t yRef(numCell); yRef.setConstant(1);

    // define ROM state
    ::pressio::ops::resize(yROM_, romSize);
    ::pressio::ops::fill(yROM_, 0.0);

    // create decoder obj
    decoder_t decoderObj(phi);

    // set indices for masking
    // only pick subset of points (randomly, this is just a test)
    auto indices = {0,1,7,8,9,14,15,18,19};

    // since this is a masked Galerkin problem, we need the collocator to
    // create a projector suitable for the masked operator
    MyCollocator mapper(indices);
    auto projector = pressio::rom::galerkin::createCollocationProjector(decoderObj,
									std::move(mapper));

    // create masker
    Masker myMask(indices);

    // create problem
    auto Problem = pressio::rom::galerkin::createMaskedResidualProblem
      <1, 2>(appobj, decoderObj, yROM_, yRef, myMask, projector);

    // linear solver
    using solver_tag	 = pressio::solvers::linear::iterative::LSCG;
    using rom_jacobian_t = decltype(Problem)::galerkin_jacobian_t;
    using linear_solver_t  = pressio::solvers::linear::Solver<solver_tag, rom_jacobian_t>;
    linear_solver_t linSolverObj;

    // nonlinear system
    auto solver = pressio::rom::galerkin::createNewtonRaphsonSolver(Problem, yROM_, linSolverObj);
    solver.setTolerance(1e-12);
    solver.setMaxIterations(3);

    // integrate in time
    pressio::rom::galerkin::solveNSteps(Problem, yROM_, 0.0, dt, 15, solver);

    // compute the fom corresponding to our rom final state
    auto yFomFinal = Problem.fomStateReconstructorCRef()(yROM_);
    fomSol_ = *yFomFinal.data();
  }
};

int main(int argc, char *argv[]){

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  std::string checkStr {"PASSED"};

  GalerkinBDF1WithResidualApi GalerkinResidApi;
  const auto residFomSol = GalerkinResidApi.fomSol_;
  const auto residRomSol = GalerkinResidApi.yROM_;

  std::cout << "check that fom reconstructed state match" << std::endl;

  // note that the gold solution is one the full mesh
  const std::vector<double> gold = {
    1.3567078458979,
    1.0098051142819,
    0.91146994459873,
    0.98959884299791,
    1.0041213990382,
    1.0051848587175,
    1.0057376621067,
    1.0063412972571,
    1.0070082159698,
    1.0077452708179,
    1.008559878282,
    1.0094601213094,
    1.0104549187023,
    1.0115545697744,
    1.0127697291376,
    1.0141127138905,
    1.0155971592115,
    1.0172371371053,
    1.0190500753216,
    1.0210535492599};

  if ((std::size_t)residFomSol.size() != gold.size())
    checkStr = "FAILED";

  // check the reconstructed fom state
  for (auto i=0; i<residFomSol.size(); i++){
    std::cout << std::setprecision(14)
	      << gold[i]
	      << " "
  	      << residFomSol[i]
  	      << std::endl;

    if (std::abs(gold[i] - residFomSol[i]) > 1e-13){
      checkStr = "FAILED";
      break;
    }
  }

  std::cout << checkStr <<  std::endl;
  return 0;
}
