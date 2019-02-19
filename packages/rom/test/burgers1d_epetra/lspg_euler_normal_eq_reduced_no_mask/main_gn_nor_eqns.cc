
#include "CORE_ALL"
#include "ODE_ALL"
#include "SVD_BASIC"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "QR_BASIC"
#include "../lspg_utils.hpp"
#include "../burgers1dEpetra_reduced_no_mask.hpp"

const std::vector<double> bdf1Sol =
  { 3.24554055320169, 1.89752376279213, 0.85479207280705,
    0.95870565961978, 1.13268225876757, 1.03283349158017,
    0.96598999130794, 0.98544746688408, 1.02979179623519,
    1.05105185279567, 1.04418409672547, 1.02642093543265,
    1.01506418485300, 1.01369384205116, 1.01660607396748,
    1.01934885770231, 1.02039498889101, 1.02041119534205,
    1.02041545687795, 1.02101238749315, 1.02209157087686,
    1.02333339700764, 1.02448690864655, 1.02552517453871,
    1.02653306205069, 1.0275820121412,  1.0286997136833,
    1.0298777957326,  1.0311068539427, 1.032383410013,
    1.0337051007736,  1.0350745588396,  1.0365026596887,
    1.0379958487248,  1.0395492054119,  1.0411610807793,
    1.0428391772815, 1.0445885471324, 1.0464078804562,
    1.0482990195458, 1.0502692141258, 1.0523226595578,
    1.0544584499389, 1.0566782919802, 1.058989288575,
    1.0613974462878, 1.0639040643275, 1.0665110949105,
    1.069225028814, 1.0720494917787};


int main(int argc, char *argv[]){

  using fom_t		= Burgers1dEpetraReducedNoMask;
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
  //assert(Comm.NumProc() == 2);

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

  // GaussNewton solver
  // hessian type: comes up in GN solver, it is (J phi)^T (J phi)
  // Since the rom is solved using eigen, hessian is wrapper of eigen matrix
  using eig_dyn_mat	 = Eigen::Matrix<scalar_t, -1, -1>;
  using hessian_t	 = rompp::core::Matrix<eig_dyn_mat>;
  using solver_tag	 = rompp::solvers::linear::LSCG;
  using converged_when_t = rompp::solvers::iterative::default_convergence;
  using gnsolver_t	 = rompp::solvers::iterative::GaussNewton<
    scalar_t, solver_tag, rompp::solvers::EigenIterative,
    converged_when_t, rom_stepper_t, hessian_t>;
  gnsolver_t solver(stGen.stepperObj_, yROM);
  solver.setTolerance(1e-14);
  solver.setMaxIterations(100);

  // integrate in time
  rompp::ode::integrateNSteps(stGen.stepperObj_, yROM, 0.0, dt, 50, solver);

  // compute the fom corresponding to our rom final state
  using fom_state_w_t = typename lspg_problem_types::fom_state_w_t;
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
