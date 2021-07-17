
#include "pressio_rom_galerkin.hpp"
#include "pressio_apps.hpp"
#include "utils_eigen.hpp"

using fom_t		= pressio::apps::Burgers1dEigen;
using scalar_t		= typename fom_t::scalar_type;
using native_state_t	= typename fom_t::state_type;
using native_velo_t	= typename fom_t::velocity_type;
using fom_state_t	= pressio::containers::Vector<native_state_t>;
using decoder_jac_t	= pressio::containers::MultiVector<Eigen::MatrixXd>;
using decoder_t		= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;
using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
using rom_state_t	= pressio::containers::Vector<eig_dyn_vec>;

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

    // create decoder obj
    decoder_t decoderObj(phi);

    // for this problem, my reference state = initial state
    native_state_t yRef(numCell);
    yRef.setConstant(1);

    // define ROM state
    ::pressio::ops::resize(yROM_, romSize);
    ::pressio::ops::fill(yROM_, 0.0);

    using ode_tag = pressio::ode::explicitmethods::Euler;
    using pressio::rom::galerkin::createDefaultProblem;
    auto Problem = createDefaultProblem<ode_tag>(appobj, decoderObj, yROM_, yRef);

    // integrate in time
    pressio::rom::galerkin::solveNSteps(Problem, yROM_, 0.0, dt, 10);

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

  GalerkinRunnerContinuousTimeApi GalerkinRun;
  const auto FomSol = GalerkinRun.fomSol_;
  const auto RomSol = GalerkinRun.yROM_;

  const std::vector<double> gold = {
    1.2397709171427,
    1.0046222847819,
    1.0025778982014,
    1.0028357431389,
    1.0031339344242,
    1.003463531189,
    1.0038277943538,
    1.0042303664469,
    1.0046752770794,
    1.0051669786797,
    1.005710393128,
    1.0063109585301,
    1.0069746851771,
    1.0077082158072,
    1.0085188933591,
    1.0094148281951,
    1.0104049892671,
    1.0114992843306,
    1.0127086671614,
    1.0140452400396};

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
