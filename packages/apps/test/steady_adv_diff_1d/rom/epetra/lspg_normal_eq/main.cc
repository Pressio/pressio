
#include "CORE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG_STEADY"
#include "APPS_STEADYADVDIFF1D"
#include "utils_epetra.hpp"
#include "../../../fom/fom_gold_states.hpp"

int main(int argc, char *argv[]){
  using fom_t		= rompp::apps::SteadyAdvDiff1dEpetra;
  using scalar_t	= typename fom_t::scalar_type;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= rompp::core::Vector<eig_dyn_vec>;

  using decoder_jac_t	= rompp::core::MultiVector<Epetra_MultiVector>;
  using decoder_t	= rompp::rom::LinearDecoder<decoder_jac_t>;

  //-------------------------------
  // MPI init
  MPI_Init(&argc,&argv);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  assert(Comm.NumProc() == 2);

  //-------------------------------
  //Parameters: diffusion, advection, expf
  std::vector<scalar_t> mu{-7, 1, 1};
  //1D spatial domain, xL, xR
  std::vector<scalar_t> domain{0, 2.0, 0.01};
  //Left and right boundary conditions
  std::vector<scalar_t> bc1D{0,0.25};
  //Create object
  fom_t  appObj(Comm, mu, domain, bc1D);
  appObj.setup();

  // number of degrees of freedom
  const int numDof = appObj.getNumGlobalNodes();
  // 21 because the app also accounts for left and right BC as unknown
  //assert(numDof == 21);

  // read the jacobian of the decoder
  constexpr int romSize = 5;
  // store modes computed before from file
  decoder_jac_t phi =
    rompp::apps::test::epetra::readBasis("basis.txt", romSize, numDof,
  					 Comm, appObj.getDataMap());
  //print to terminal the basis
  //phi.data()->Print(std::cout);
  // decoder object
  decoder_t decoderObj(phi);

  // my reference state
  auto yRef = appObj.getState();

  // define ROM state
  lspg_state_t yROM(romSize);
  // initialize to zero (for now)
  yROM.putScalar(0.0);

  // define LSPG type
  using lspg_problem_type = rompp::rom::DefaultLSPGSteadyTypeGenerator<
    fom_t, decoder_t, lspg_state_t>;
  rompp::rom::LSPGSteadyProblemGenerator<lspg_problem_type> lspgProblem(
      appObj, *yRef, decoderObj, yROM);

  // GaussNewton solver
  // hessian comes up in GN solver, it is (J phi)^T (J phi)
  // rom is solved using eigen, hessian is wrapper of eigen matrix
  using eig_dyn_mat	 = Eigen::Matrix<scalar_t, -1, -1>;
  using hessian_t	 = rompp::core::Matrix<eig_dyn_mat>;
  using solver_tag	 = rompp::solvers::linear::LSCG;
  using converged_when_t = rompp::solvers::iterative::default_convergence;
  using rom_system_t = typename lspg_problem_type::lspg_system_t;
  using gnsolver_t	 = rompp::solvers::iterative::GaussNewton<
    scalar_t, solver_tag, rompp::solvers::EigenIterative,
    converged_when_t, rom_system_t, hessian_t>;
  gnsolver_t solver(lspgProblem.systemObj_, yROM);
  solver.setTolerance(1e-6);
  solver.setMaxIterations(200);
  solver.solve(lspgProblem.systemObj_, yROM);

  // compute the fom corresponding to our rom final state
  auto yFomFinal = lspgProblem.yFomReconstructor_(yROM);
  yFomFinal.data()->Print(std::cout << std::setprecision(14));

  // /* if there is a reproducing test we can do, let's do it and we can
  //  * compare with the FOM solution.
  //  * Otherwise, we can compare ROM solution in some other way (?)
  //  * I leave the final check empty for now.
  //  */

  MPI_Finalize();
  return 0;
}
