
#include "pressio_rom_lspg.hpp"
#include "pressio_apps.hpp"
#include "utils_epetra.hpp"
#include "utils_epetra_identity_preconditioner.hpp"

int main(int argc, char *argv[])
{
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
  constexpr scalar_t Pr = 2.48875378592597;
  constexpr scalar_t Re = 41.7029840887032;

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
  using prec_t = pressio::rom::test::EpetraIdentityPreconditioner;
  prec_t Prec;

  // using lspg_problem_type = 
  //   typename pressio::rom::lspg::composePreconditionedDefaultProblem<
  //   fom_adapter_t, decoder_t, lspg_state_t, prec_t>::type;
  // lspg_problem_type lspgProblem(appObjROM, *yRef, decoderObj, yROM, Prec);
  auto lspgProblem = pressio::rom::lspg::createPreconditionedDefaultProblemSteady(
    appObjROM, decoderObj, yROM, *yRef, Prec);

  // linear solver
  using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
  using hessian_t  = pressio::containers::DenseMatrix<eig_dyn_mat>;
  using solver_tag   = pressio::solvers::linear::iterative::LSCG;
  using linear_solver_t = pressio::solvers::linear::Solver<solver_tag, hessian_t>;
  linear_solver_t linSolverObj;

  // GaussNewton solver
  auto solver = pressio::solvers::nonlinear::createGaussNewton(
      lspgProblem.systemRef(), yROM, linSolverObj);
  solver.setTolerance(1e-14);
  solver.setMaxIterations(200);
  solver.solve(lspgProblem.systemRef(), yROM);

  /* the ROM is run for a parameter point that was used to generate
   * the basis, so we should recover the FOM solution exactly */
  // reconstruct the fom corresponding to our rom final state
  auto yFomApprox = lspgProblem.fomStateReconstructorCRef()(yROM);
  appObjROM.printStateToFile("rom.txt", *yFomApprox.data());
  auto errorVec(yFom); 
  pressio::ops::update(errorVec, yFom, 1., yFomApprox, -1.);
  const auto norm2err = pressio::ops::norm2(errorVec);
  if( norm2err > 1e-12 ) checkStr = "FAILED";

  std::cout << std::setprecision(15) << norm2err << std::endl;

  MPI_Finalize();
  std::cout << checkStr <<  std::endl;
  return 0;
}
