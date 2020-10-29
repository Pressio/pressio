
#include "pressio_rom_galerkin.hpp"
#include "pressio_apps.hpp"
#include "utils_epetra.hpp"

const std::vector<double> bdf1Sol
{5.0081542681376, 5.016629490569, 5.025433912557, 5.0345792953115,
 5.0440827179355, 5.0539568295087, 5.0642107801363, 5.074857742734,
 5.0859146000515, 5.0974001265619, 5.1093302710968, 5.1217197481536,
 5.1345846667293, 5.1479436063682, 5.1618137609004, 5.1762071980595,
 5.1911395190849, 5.2066322357211, 5.222706587389, 5.2393822195142,
 5.2566784890019, 5.274617970535, 5.2932246323729, 5.3125186218141,
 5.3325236627322, 5.3532729201416, 5.3747971779128, 5.3971189932731,
 5.4202577535351, 5.4442348269811, 5.469078757402, 5.4948202159562,
 5.5214859714822, 5.5491009348394, 5.5776911098501, 5.6072849195866,
 5.6379131952825, 5.6696069037791, 5.7023980878343, 5.7363239274031,
 5.7714263431002, 5.807744410524, 5.8453128737429, 5.884168702448,
 5.9243510856362, 5.9658923478856, 6.0088164545724, 6.0531503069487,
 6.0989210765093, 6.1461565470309};

 int main(int argc, char *argv[]){

  using fom_t		= pressio::apps::Burgers1dEpetra;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;
  using fom_state_t  = pressio::containers::Vector<native_state_t>;

  using eig_dyn_vec = Eigen::Matrix<scalar_t, -1, 1>;
  using rom_state_t = pressio::containers::Vector<eig_dyn_vec>;

  using decoder_jac_t	= pressio::containers::MultiVector<Epetra_MultiVector>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

  std::string checkStr {"PASSED"};

  //-------------------------------
  // MPI init
  MPI_Init(&argc,&argv);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  if(Comm.NumProc() != 2) return 0;

  //-------------------------------
  // app object
  int numCell = 50;
  fom_t appobj({5.0, 0.02, 0.02}, numCell, &Comm);
  scalar_t dt = 0.01;

  // store (whichever way you want) the jacobian of the decoder
  constexpr int romSize = 20;
  decoder_jac_t phi =
    pressio::rom::test::epetra::readBasis("basis.txt", romSize, numCell,
					 Comm, appobj.getDataMap());
  decoder_t decoderObj(phi);

  // this is my reference state
  auto & y0n = appobj.getInitialState();

  // define ROM state
  rom_state_t yROM(romSize);
  // initialize to zero (this has to be done)
  pressio::ops::fill(yROM, 0.0);

  using ode_tag = pressio::ode::explicitmethods::Euler;
  // using problem_t =
  // pressio::rom::galerkin::composeDefaultProblem<ode_tag, fom_t, decoder_t, rom_state_t>::type;
  // problem_t galerkinProb(appobj, decoderObj, yROM, y0n);
  auto galerkinProb =
    pressio::rom::galerkin::createDefaultProblem<ode_tag>(appobj, decoderObj, yROM, y0n);

  scalar_t fint = 35;
  auto nSteps = static_cast<::pressio::ode::types::step_t>(fint/dt);
  pressio::ode::advanceNSteps(galerkinProb.stepperRef(), yROM, 0.0, dt, nSteps);

  std::cout << *yROM.data() << std::endl;

  // compute the fom corresponding to our rom final state
  auto yFomFinal = galerkinProb.fomStateReconstructorCRef()(yROM);
  yFomFinal.data()->Print(std::cout << std::setprecision(14));

  // check against gold solution
  int shift = 0;
  if (rank==1)  shift = 25;
  int myn = yFomFinal.data()->Map().NumMyElements();
  for (auto i=0; i<myn; i++){
    if(std::abs(yFomFinal[i] - bdf1Sol[i+shift]) > 1e-12 ){
      checkStr = "FAILED";
      break;
    }
  }

  MPI_Finalize();
  std::cout << checkStr <<  std::endl;
  return 0;
}
