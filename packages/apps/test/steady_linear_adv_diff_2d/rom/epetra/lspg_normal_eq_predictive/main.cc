
#include "CORE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG_STEADY"
#include "APPS_STEADYLINADVDIFF2D"
#include "utils_epetra.hpp"

int main(int argc, char *argv[]){
  using true_fom_t	= rompp::apps::SteadyLinAdvDiff2dEpetra;
  using fom_adapter_t	= rompp::apps::SteadyLinAdvDiff2dEpetraRomAdapter;
  using scalar_t	= typename fom_adapter_t::scalar_type;
  using native_state	= typename fom_adapter_t::state_type;

  // MPI
  MPI_Init(&argc,&argv);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  //assert(Comm.NumProc() == 2);

  // we run the FOM and LSPG for same values of parameters
  constexpr scalar_t Pr = 4.;
  constexpr scalar_t Re = 95.;

  // the discretization to use for solver
  const int Nx = 11, Ny = Nx*2-1;
  // we solve for inner grid along x
  const int numDof = (Nx-2)*Ny;

  // -------------------------
  // run FOM model

  true_fom_t  appObj(Comm, Nx, Ny, Pr, Re);
  appObj.setup();
  appObj.assembleMatrix();
  appObj.fillRhs();
  appObj.solve();
  appObj.printStateToFile("fom.txt");

  rompp::core::Vector<native_state> yFom(*appObj.getState());

  // -------------------------
  // LSPG ROM
  using native_state	= typename fom_adapter_t::state_type;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= rompp::core::Vector<eig_dyn_vec>;
  using decoder_jac_t	= rompp::core::MultiVector<Epetra_MultiVector>;
  using decoder_t	= rompp::rom::LinearDecoder<decoder_jac_t>;

  constexpr int romSize = 5;

  // app object for running rom
  fom_adapter_t  appObjROM(Comm, Nx, Ny, Pr, Re);

  // store modes computed before from file
  const decoder_jac_t phi =
    rompp::apps::test::epetra::readBasis("basis.txt", romSize, numDof,
  					 Comm, appObjROM.getDataMap());
  // decoder object
  decoder_t decoderObj(phi);

  // my reference state
  auto yRef = appObjROM.getState();
  yRef->PutScalar( static_cast<scalar_t>(0) );

  // define ROM state and set to zero
  lspg_state_t yROM(romSize);
  yROM.putScalar(0.0);

  // define LSPG type
  using lspg_problem_type = rompp::rom::DefaultLSPGSteadyTypeGenerator<
    fom_adapter_t, decoder_t, lspg_state_t>;
  rompp::rom::LSPGSteadyProblemGenerator<lspg_problem_type> lspgProblem(
      appObjROM, *yRef, decoderObj, yROM);

  // GaussNewton solver
  // hessian comes up in GN solver, it is (J phi)^T (J phi)
  // rom is solved using eigen, hessian is wrapper of eigen matrix
  using eig_dyn_mat	 = Eigen::Matrix<scalar_t, -1, -1>;
  using hessian_t	 = rompp::core::Matrix<eig_dyn_mat>;
  using solver_tag	 = rompp::solvers::linear::LSCG;
  using converged_when_t = rompp::solvers::iterative::default_convergence;
  using rom_system_t	 = typename lspg_problem_type::lspg_system_t;
  using gnsolver_t	 = rompp::solvers::iterative::GaussNewton<
    scalar_t, solver_tag, rompp::solvers::EigenIterative,
    converged_when_t, rom_system_t, hessian_t>;
  gnsolver_t solver(lspgProblem.systemObj_, yROM);
  solver.setTolerance(1e-14);
  solver.setMaxIterations(200);
  solver.solve(lspgProblem.systemObj_, yROM);

  // reconstruct the fom corresponding to our rom final state
  auto yFomApprox = lspgProblem.yFomReconstructor_(yROM);
  appObjROM.printStateToFile("rom.txt", *yFomApprox.data());

  /* this is a predictive run, so we should recover FOM
   * solution only approximately */
  auto normFomY = rompp::core::ops::norm2(yFom);
  auto errorVec(yFom); errorVec = yFom - yFomApprox;
  const auto norm2err = rompp::core::ops::norm2(errorVec);
  assert( (norm2err/normFomY)*100 < 0.1 ); // less than 0.1 %
  std::cout << std::setprecision(15) <<
    "% relative error (L2-norm) " << " " <<
    norm2err/normFomY*100. << std::endl;

  MPI_Finalize();
  return 0;
}
