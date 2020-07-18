
#include "pressio_rom.hpp"
#include "pressio_apps.hpp"
#include "utils_epetra.hpp"

int main(int argc, char *argv[]){
  using true_fom_t	= pressio::apps::SteadyLinAdvDiff2dEpetra;
  using fom_adapter_t	= pressio::apps::SteadyLinAdvDiff2dEpetraRomAdapter;
  using scalar_t	= typename fom_adapter_t::scalar_type;
  using native_state	= typename fom_adapter_t::state_type;

  std::string checkStr {"PASSED"};

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
  appObj.assembleMatrix();
  appObj.fillRhs();
  appObj.solve();
  appObj.printStateToFile("fom.txt");

  pressio::containers::Vector<native_state> yFom(*appObj.getState());

  // -------------------------
  // LSPG ROM
  using native_state	= typename fom_adapter_t::state_type;
  using fom_state_t = pressio::containers::Vector<native_state>;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;
  using decoder_jac_t	= pressio::containers::MultiVector<Epetra_MultiVector>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

  constexpr int romSize = 5;

  // app object for running rom
  fom_adapter_t  appObjROM(Comm, Nx, Ny, Pr, Re);

  // store modes computed before from file
  const decoder_jac_t phi =
    pressio::rom::test::epetra::readBasis("basis.txt", romSize, numDof,
  					 Comm, appObjROM.getDataMap());
  // decoder object
  decoder_t decoderObj(phi);

  // my reference state
  auto yRef = appObjROM.getState();
  yRef->PutScalar( static_cast<scalar_t>(0) );

  // define ROM state and set to zero
  lspg_state_t yROM(romSize);
  pressio::ops::fill(yROM, 0.0);

  // define LSPG type
  using lspg_problem_type = pressio::rom::lspg::composeDefaultProblem<
      fom_adapter_t, lspg_state_t, decoder_t>::type;
  lspg_problem_type lspgProblem(appObjROM, *yRef, decoderObj);  
  using rom_system_t = typename lspg_problem_type::lspg_system_t;

  // linear solver
  using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
  using hessian_t  = pressio::containers::Matrix<eig_dyn_mat>;
  using solver_tag   = pressio::solvers::linear::iterative::LSCG;
  using linear_solver_t = pressio::solvers::linear::Solver<solver_tag, hessian_t>;
  linear_solver_t linSolverObj;

  // GaussNewton solver
  using nls_t = pressio::solvers::nonlinear::composeGaussNewton_t<
    rom_system_t,
    pressio::solvers::nonlinear::DefaultUpdate,
    pressio::solvers::nonlinear::StopWhenCorrectionNormBelowTol,
    linear_solver_t>;
  nls_t solver(lspgProblem.getSystemRef(), yROM, linSolverObj);
  solver.setTolerance(1e-14);
  solver.setMaxIterations(200);
  solver.solve(lspgProblem.getSystemRef(), yROM);

  // reconstruct the fom corresponding to our rom final state
  auto yFomApprox = lspgProblem.getFomStateReconstructorCRef()(yROM);
  appObjROM.printStateToFile("rom.txt", *yFomApprox.data());

  /* this is a predictive run, so we should recover FOM
   * solution only approximately */
  auto normFomY = pressio::ops::norm2(yFom);
  auto errorVec(yFom); 
  pressio::ops::do_update(errorVec, yFom, 1., yFomApprox, -1.);
  const auto norm2err = pressio::ops::norm2(errorVec);
  if( (norm2err/normFomY)*100 > 0.1 ) checkStr = "FAILED";

  std::cout << std::setprecision(15) <<
    "% relative error (L2-norm) " << " " <<
    norm2err/normFomY*100. << std::endl;

  MPI_Finalize();
  std::cout << checkStr <<  std::endl;
  return 0;
}
