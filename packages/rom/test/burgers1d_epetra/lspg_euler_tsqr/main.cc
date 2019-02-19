
#include "CORE_ALL"
#include "ODE_ALL"
#include "SVD_BASIC"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "QR_BASIC"
#include "../lspg_utils.hpp"
#include "../burgers1dEpetra.hpp"

const std::vector<double> bdf1Sol =
  { 4.9683703619639, 4.6782198130332, 3.7809092143183, 2.3410370171735,
    1.3665227552633, 1.0998248856967, 1.0588744517773, 1.0552786290916,
    1.0565755841731, 1.05832566338, 1.0598499804131, 1.0612950695346,
    1.0634140889979, 1.0664875741963, 1.0701215280868, 1.0735734095262,
    1.0764754781593, 1.0790691949921, 1.0818228048193, 1.0850683068699,
    1.0887478932688, 1.092642408022, 1.0965446003949, 1.1004297242319,
    1.1043922522415, 1.1085490607723, 1.1129611388185, 1.117613658291,
    1.1224648232926, 1.1274778902247, 1.1326530965013, 1.13802239319,
    1.1436307647296, 1.1494987675359, 1.1556064137176, 1.161941826502,
    1.1685253861164, 1.1753816797655, 1.1825238154096, 1.1899594186277,
    1.1976959450073, 1.205743187805, 1.2141142750508, 1.2228251947121,
    1.2318925753159, 1.2413295315044, 1.2511467170661, 1.2613604601936,
    1.2719910028481, 1.2830509792253};

int main(int argc, char *argv[]){

  using fom_t		= Burgers1dEpetra;
  using scalar_t	= typename fom_t::scalar_type;
  using decoder_jac_t	= rompp::core::MultiVector<Epetra_MultiVector>;
  using decoder_t	= rompp::rom::LinearDecoder<decoder_jac_t>;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= rompp::core::Vector<eig_dyn_vec>;

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
  const int numBasis = phi.globalNumVectors();
  assert( numBasis == romSize );

  // this is my reference state
  auto & y0n = appobj.getInitialState();
  decoder_t decoderObj(phi);

  // define ROM state
  lspg_state_t yROM(romSize);
  // initialize to zero (this has to be done)
  yROM.putScalar(0.0);

  using lspg_problem_types = rompp::rom::DefaultLSPGTypeGenerator<
    fom_t, rompp::ode::ImplicitEnum::Euler, decoder_t, lspg_state_t>;
  rompp::rom::LSPGStepperObjectGenerator<lspg_problem_types> stGen(
      appobj, y0n, decoderObj, yROM, t0);

  using rom_stepper_t = typename lspg_problem_types::rom_stepper_t;
  using fom_state_w_t = typename lspg_problem_types::fom_state_w_t;
  using rom_jac_t     = typename lspg_problem_types::lspg_matrix_t;

  // GaussNewton solver
  using qr_algo = rompp::qr::TSQR;
  using qr_type = rompp::qr::QRSolver<rom_jac_t, qr_algo>;
  using converged_when_t = rompp::solvers::iterative::default_convergence;
  using gnsolver_t	 = rompp::solvers::iterative::GaussNewtonQR<
  				scalar_t, qr_type,
  				converged_when_t, rom_stepper_t>;
  gnsolver_t solver(stGen.stepperObj_, yROM);
  solver.setTolerance(1e-14);
  solver.setMaxIterations(100);

  // integrate in time
  rompp::ode::integrateNSteps(stGen.stepperObj_, yROM, 0.0, dt, 200, solver);

  // compute the fom corresponding to our rom final state
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

  MPI_Finalize();
  return 0;
}
//------------------------------
// end
//------------------------------
