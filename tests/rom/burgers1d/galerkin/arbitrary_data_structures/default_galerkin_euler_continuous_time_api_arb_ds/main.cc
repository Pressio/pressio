
#include "pressio_apps.hpp"
#include "pressio_rom_galerkin.hpp"
#include "utils_eigen.hpp"

namespace
{

const std::vector<double> goldFom = {5.0081542681376, 5.016629490569,
				     5.025433912557, 5.0345792953115,
				     5.0440827179355, 5.0539568295087,
				     5.0642107801363, 5.074857742734,
				     5.0859146000515,5.0974001265619,
				     5.1093302710968,5.1217197481536,
				     5.1345846667293, 5.1479436063682,
				     5.1618137609004, 5.1762071980595,
				     5.1911395190849, 5.2066322357211,
				     5.222706587389,5.2393822195142,
				     5.2566784890019, 5.274617970535,
				     5.2932246323729, 5.3125186218141,
				     5.3325236627322,5.3532729201416,
				     5.3747971779128, 5.3971189932731,
				     5.4202577535351, 5.4442348269811,
				     5.469078757402, 5.4948202159561,
				     5.5214859714822,5.5491009348394,
				     5.5776911098501, 5.6072849195866,
				     5.6379131952825,5.6696069037791,
				     5.7023980878343, 5.7363239274031,
				     5.7714263431002,5.807744410524,
				     5.8453128737429,5.884168702448,
				     5.9243510856362,5.9658923478856,
				     6.0088164545724, 6.0531503069487,
				     6.0989210765093, 6.1461565470309};

template <typename result_t>
void readBasis(std::string filename, result_t & phi)
{
  const auto nRows = phi.extent(0);
  const auto nCols = phi.extent(1);

  std::vector<std::vector<double>> A0;
  ::pressio::utils::readAsciiMatrixStdVecVec(filename, A0, nCols);
  for (std::size_t i=0; i<nRows; i++){
    for (std::size_t j=0; j<nCols; j++)
      phi(i,j) = A0[i][j];
  }
}

template <typename sc_t, typename dec_jac_t>
struct myOps
{
  //  y = beta * y + alpha*A*x
  template <typename x_t>
  void product(::pressio::nontranspose mode,
	       const sc_t alpha,
	       const dec_jac_t & A,
	       const x_t & x,
	       const sc_t beta,
	       pressio::apps::arbds::Vector<sc_t> & y) const
  {
    // x is subscriptable like a regular array, e.g. you can do x[i] or x(i)
    const auto nArows = A.extent(0);
    const auto nAcols = A.extent(1);
    for (std::size_t i=0; i<nArows; ++i)
    {
      y(i) = beta*y(i);
      for (std::size_t j=0; j<nAcols; ++j)
	y(i) += alpha * A(i,j) * x(j);
    }
  }

  //  y = beta * y + alpha*A^T*x
  template <typename y_t>
  void product(::pressio::transpose mode,
	       const sc_t alpha,
	       const dec_jac_t & A,
	       const pressio::apps::arbds::Vector<sc_t> & x,
	       const sc_t beta,
	       y_t & y) const
  {
    // y is subscriptable like a regular array
    const auto nArows = A.extent(0);
    const auto nAcols = A.extent(1);
    for (std::size_t j=0; j<nAcols; ++j){
      y(j) = beta*y(j);
      for (std::size_t i=0; i<nArows; ++i)
    	y(j) += alpha * A(i,j) * x(i);
    }
  }

  void deep_copy(pressio::apps::arbds::Vector<sc_t> & to,
  		 const pressio::apps::arbds::Vector<sc_t> & from) const
  {
    // here you need do a deep copy from -> to
    to = from;
  }

  void set_zero(pressio::apps::arbds::Vector<sc_t> & vec) const
  {
    for (std::size_t i=0; i<vec.extent(0); ++i)
     vec(i) = static_cast<sc_t>(0);
  }

  void axpy(sc_t alpha,
	    const pressio::apps::arbds::Vector<sc_t> & x,
	    pressio::apps::arbds::Vector<sc_t> & y) const
  {
    // compute y = y + alfa * x
    for (std::size_t i=0; i<y.extent(0); ++i)
      y(i) += alpha * x(i);
  }
};

struct EulerGalerkinWithVelocityApi
{
  using fom_t		= pressio::apps::Burgers1dArbDsContinuousTimeApiAdapter;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;
  using native_dmat_t   = typename fom_t::dense_matrix_type;

  using ops_t		= myOps<scalar_t, native_dmat_t>;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using rom_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using fom_state_t	= pressio::containers::Vector<native_state_t>;
  using decoder_jac_t	= pressio::containers::MultiVector<native_dmat_t>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t, ops_t>;

  static_assert(pressio::rom::concepts::continuous_time_system_without_user_provided_apply_jacobian<fom_t>::value, "");
  static_assert(pressio::rom::concepts::continuous_time_system<fom_t>::value, "");

  static_assert(pressio::rom::concepts::custom_ops_for_fom_state_reconstructor<
		ops_t, fom_state_t>::value, "");
  static_assert(pressio::rom::galerkin::concepts::custom_ops_continuous_time<
		ops_t, decoder_jac_t, rom_state_t, fom_state_t>::value, "");

  native_state_t fomSol_ = {};
  rom_state_t yROM_ = {};

  EulerGalerkinWithVelocityApi()
  {
    std::string checkStr {"PASSED"};

    // app object
    constexpr int numCell = 50;
    pressio::apps::Burgers1dArbDs appObj(numCell);
    // adapter
    fom_t fomObj(appObj);
    scalar_t dt = 0.01;

    ops_t myOps;

    // read from file the jacobian of the decoder
    constexpr int romSize = 20;
    // store modes from file
    decoder_jac_t phi(numCell, romSize);
    readBasis("basis.txt", phi);

    // create decoder obj
    decoder_t decoderObj(phi, myOps);

    // this is my reference state
    auto y0n = appObj.getInitialState();

    // for this problem, my reference state = initial state
    native_state_t yRef(numCell);
    for (std::size_t i=0; i<yRef.extent(0); ++i)
      yRef(i) = pressio::utils::constants<scalar_t>::one();

    // define ROM state
    pressio::ops::resize(yROM_, romSize);
    pressio::ops::fill(yROM_, pressio::utils::constants<scalar_t>::zero());

    static_assert
    (pressio::rom::concepts::continuous_time_system_without_user_provided_apply_jacobian<fom_t>::value,"");

    using ode_tag = pressio::ode::explicitmethods::Euler;
    // using problem_t  = pressio::rom::galerkin::composeDefaultProblem<
    //   ode_tag, fom_t, decoder_t, rom_state_t, ops_t>::type;
    // problem_t galerkinProb(fomObj, decoderObj, yROM_, y0n, myOps);
    auto galerkinProb =
      pressio::rom::galerkin::createDefaultProblem<ode_tag>
      (fomObj, decoderObj, yROM_, y0n, myOps);

    scalar_t fint = 35;
    auto nSteps = static_cast<::pressio::ode::types::step_t>(fint/dt);
    pressio::rom::galerkin::solveNSteps(galerkinProb, yROM_, 0.0, dt, nSteps);

    // compute the fom corresponding to our rom final state
    auto yFomFinal = galerkinProb.fomStateReconstructorCRef()(yROM_);
    fomSol_ = *yFomFinal.data();
  }
};
}

int main(int argc, char *argv[]){

  std::string checkStr {"PASSED"};
  EulerGalerkinWithVelocityApi probObj;

  std::cout << "check that fom reconstructed state match" << std::endl;
  // check the reconstructed fom state
  for (std::size_t i=0; i<probObj.fomSol_.extent(0); i++){
    std::cout << std::setprecision(14)
  	      << goldFom[i]
  	      << " "
  	      << probObj.fomSol_(i)
  	      << std::endl;

    if (std::abs(goldFom[i] - probObj.fomSol_(i)) > 1e-13)
      checkStr = "FAILED";
  }
  std::cout << checkStr <<  std::endl;

  return 0;
}
