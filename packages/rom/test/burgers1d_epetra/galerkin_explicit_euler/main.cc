
#include "CORE_ALL"
#include "ODE_ALL"
#include "SVD_BASIC"
#include "SOLVERS_NONLINEAR"
#include "ROM_GALERKIN"
#include "QR_BASIC"
#include "../lspg_utils.hpp"
#include "../burgers1dEpetra.hpp"

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

  using fom_t		= Burgers1dEpetra;
  using scalar_t	= typename fom_t::scalar_type;
  using decoder_jac_t	= rompp::core::MultiVector<Epetra_MultiVector>;
  using decoder_t	= rompp::rom::LinearDecoder<decoder_jac_t>;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using rom_state_t	= rompp::core::Vector<eig_dyn_vec>;

  //-------------------------------
  // MPI init
  MPI_Init(&argc,&argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  assert(Comm.NumProc() == 2);

  //-------------------------------
  // app object
  int numCell = 50; // # fv cells
  std::vector<double> mu({5.0, 0.02, 0.02});
  fom_t appobj(mu, numCell, &Comm);
  appobj.setup();
  auto t0 = static_cast<scalar_t>(0);
  scalar_t dt = 0.01;

  // store (whichever way you want) the jacobian of the decoder
  constexpr int romSize = 20;
  // store modes computed before from file
  decoder_jac_t phi = rompp::rom::test::getBasis(romSize, numCell, Comm,
  						 appobj.getDataMap());

  // this is my reference state
  auto & y0n = appobj.getInitialState();
  decoder_t decoderObj(phi);

  // define ROM state
  rom_state_t yROM(romSize);
  // initialize to zero (this has to be done)
  yROM.putScalar(0.0);

  using galerkin_types = rompp::rom::DefaultGalerkinTypeGenerator<
    fom_t, rompp::ode::ExplicitEnum::Euler, decoder_t, rom_state_t>;
  rompp::rom::GalerkinExplicitStepperObjectGenerator<galerkin_types> stGen(
      appobj, y0n, decoderObj, yROM, t0);

  scalar_t fint = 35;
  auto nSteps = static_cast<unsigned int>(fint/dt);
  rompp::ode::integrateNSteps(stGen.stepperObj_, yROM, 0.0, dt, nSteps);

  // compute the fom corresponding to our rom final state
  using fom_state_w_t = typename galerkin_types::fom_state_w_t;
  fom_state_w_t yRf(y0n);
  decoderObj.applyMapping(yROM, yRf);
  yRf += stGen.y0Fom_;
  yRf.data()->Print(std::cout << std::setprecision(14));

  // check against gold solution
  int shift = 0;
  if (rank==1)  shift = 25;
  int myn = yRf.getDataMap().NumMyElements();
  for (auto i=0; i<myn; i++){
    assert(std::abs(yRf[i] - bdf1Sol[i+shift]) < 1e-12 );
   }

  // // print summary from timers
  // #ifdef HAVE_TEUCHOS_TIMERS
  // rompp::core::TeuchosPerformanceMonitor::stackedTimersReportMPI();
  // #endif

  MPI_Finalize();
  return 0;
}
//------------------------------
// end
//------------------------------
