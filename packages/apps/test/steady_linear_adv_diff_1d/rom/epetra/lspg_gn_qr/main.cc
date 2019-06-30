
#include "CONTAINERS_ALL"
#include "SOLVERS_NONLINEAR"
#include "QR_BASIC"
#include "ROM_LSPG_STEADY"
#include "APPS_STEADYLINADVDIFF1D"
#include "utils_epetra.hpp"
#include "../../../fom/fom_gold_states.hpp"

int main(int argc, char *argv[]){
  using fom_t		= rompp::apps::SteadyLinAdvDiff1dEpetra;
  using scalar_t	= typename fom_t::scalar_type;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= rompp::containers::Vector<eig_dyn_vec>;

  using decoder_jac_t	= rompp::containers::MultiVector<Epetra_MultiVector>;
  using decoder_t	= rompp::rom::LinearDecoder<decoder_jac_t>;
  using native_state    = typename fom_t::state_type;

  std::string checkStr {"PASSED"};

  //-------------------------------
  // MPI init
  MPI_Init(&argc,&argv);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  if(Comm.NumProc() != 2) return 0;

  //-------------------------------
  //Parameters: diffusion, advection, expf
  std::vector<scalar_t> mu{-2.194214766762222, 0.984182281940700, 2.729726020165052};
  //1D spatial domain, xL, xR
  std::vector<scalar_t> domain{0, 2.0, 0.01};
  //Left and right boundary conditions
  std::vector<scalar_t> bc1D{0,0.25};
  //Create object
  fom_t  appObj(Comm, mu, domain, bc1D);
  appObj.setup();

  // run FOM model
  appObj.setup();
  appObj.calculateLinearSystem();
  appObj.calculateForcingTerm();
  appObj.solve();
  rompp::containers::Vector<native_state> yFom(*appObj.getState());

  // number of degrees of freedom
  const int numDof = appObj.getNumGlobalNodes();

  // read the jacobian of the decoder
  constexpr int romSize = 15;
  // store modes computed before from file
  decoder_jac_t phi =
    rompp::apps::test::epetra::readBasis("basis.txt", romSize, numDof,
  					 Comm, appObj.getDataMap());
  //print to terminal the basis
  phi.data()->Print(std::cout);
  // decoder object
  decoder_t decoderObj(phi);

  // my reference state
  auto yRef = appObj.getState();
  yRef->PutScalar( static_cast<scalar_t>(0) );

  // define ROM state
  lspg_state_t yROM(romSize);
  // initialize to zero (for now)
  yROM.putScalar(0.0);

  // define LSPG type
  using lspg_problem_type = rompp::rom::DefaultLSPGSteadyTypeGenerator<
    fom_t, decoder_t, lspg_state_t>;
  rompp::rom::LSPGSteadyProblemGenerator<lspg_problem_type> lspgProblem(
      appObj, *yRef, decoderObj, yROM);

  using rom_sys_t     = typename lspg_problem_type::lspg_system_t;
  using rom_jac_t     = typename lspg_problem_type::lspg_matrix_t;

  // GaussNewton solver
  using rom_jac_t     = typename lspg_problem_type::lspg_matrix_t;
  using qr_algo = rompp::qr::TSQR;
  using qr_type = rompp::qr::QRSolver<rom_jac_t, qr_algo>;
  using converged_when_t = rompp::solvers::iterative::default_convergence;
  using gnsolver_t  = rompp::solvers::iterative::GaussNewtonQR<
         rom_sys_t, qr_type, converged_when_t>;
  gnsolver_t solver(lspgProblem.systemObj_, yROM);
  solver.setTolerance(1e-11);
  solver.setMaxIterations(2);
  solver.solve(lspgProblem.systemObj_, yROM);

  // compute the fom corresponding to our rom final state
  auto yFomFinal = lspgProblem.yFomReconstructor_(yROM);
  yFomFinal.data()->Print(std::cout << std::setprecision(14));

  auto errorVec(yFom);
  errorVec = yFom-yFomFinal;
  const auto norm2err = rompp::containers::ops::norm2(errorVec);

  if (norm2err > 1e-10) checkStr = "FAILED";
  std::cout << std::setprecision(15) << norm2err << std::endl;

  MPI_Finalize();
  std::cout << checkStr << std::endl;
  return 0;
}
