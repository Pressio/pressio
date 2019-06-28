
#include "ALGEBRA_ALL"
#include "ODE_ALL"
#include "ROM_LSPG"
#include "SOLVERS_NONLINEAR"
#include "APPS_UNSTEADYLINADVDIFF1D"
#include "utils_epetra.hpp"
#include "../../../fom/gold_states_implicit.hpp"

constexpr double eps = 1e-7;
std::string checkStr = "PASSED";

template <typename T>
void checkSol(int rank, const T &y,
	      const std::vector<double> & trueS){
  int shift = (rank==0) ? 0 : 20;
  for (decltype(y.localSize()) i=0; i<y.localSize(); i++){
    if (std::abs(y[i] - trueS[i+shift]) > eps) checkStr = "FAILED";
  }
}

int main(int argc, char *argv[]){
  using fom_t          = rompp::apps::UnsteadyLinAdvDiff1dEpetra;
  using scalar_t       =typename fom_t::scalar_type;
  using  eig_dyn_vec   = Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t   = rompp::algebra::Vector<eig_dyn_vec>;
  using decoder_jac_t  = rompp::algebra::MultiVector<Epetra_MultiVector>;
  using decoder_t      = rompp::rom::LinearDecoder<decoder_jac_t>;
  using native_state    = typename fom_t::state_type;
  using app_state_t     = typename fom_t::state_type;
  using app_residual_t = typename fom_t::residual_type;
  using app_jacobian_t = typename fom_t::jacobian_type;

  //----------------------------------------------------------------------
  // MPI initialization
  //----------------------------------------------------------------------
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  assert(Comm.NumProc() == 2);

  //----------------------------------------------------------------------
  // Designate the paramerters and problem
  //----------------------------------------------------------------------
  const scalar_t dt = 0.1;
  const auto Nsteps = 50;//static_cast<unsigned int>(fint/dt);
  const scalar_t fint = dt*Nsteps;

  // Actual parameter inputs that were used to formulate the basis
  const std::vector<scalar_t> mu{-0.857241631161166, 0.104833925269630,
      -0.713183149274631};
  const std::vector<scalar_t> domain{0, 2.0, 0.05};
  const std::vector<scalar_t> bc1D{0, 0.25};

  //----------------------------------------------------------------------
  // App for UnsteadyLinAdvDiff1dEpetra Object
  //----------------------------------------------------------------------
  fom_t appObj(Comm, mu, domain, bc1D);
  appObj.unsteadySetup();

  //----------------------------------------------------------------------
  // ROM parameters and Basis functions
  //----------------------------------------------------------------------
  constexpr int romSize = 15;
  constexpr auto ode_case = rompp::ode::ImplicitEnum::Euler;
  auto t0 = static_cast<scalar_t>(0);
  const int numDof = appObj.getNumGlobalNodes();
  std::cout << numDof << "x" << romSize <<std::endl;

  //Retrieve basis found from an external source (MATLAB)
  const decoder_jac_t phi =
    rompp::apps::test::epetra::readBasis("basis.txt", romSize, numDof,
					 Comm, appObj.getDataMap());
  phi.data()->Print(std::cout);
  decoder_t decoderObj(phi);

  //Retrieve the initial states defined specifically for this problem
  auto yRef =  *appObj.getInitialState();
  lspg_state_t yROM(romSize);

  //initialize to zero since using y-y0 notation
  yROM.putScalar(0.0);

  //----------------------------------------------------------------------
  // Define LSPG Type
  //----------------------------------------------------------------------
  using lspg_problem_types = rompp::rom::DefaultLSPGTypeGenerator<
    fom_t, ode_case, decoder_t, lspg_state_t>;
  rompp::rom::LSPGUnsteadyProblemGenerator<lspg_problem_types>
    lspgProblem(appObj, yRef, decoderObj, yROM, t0);

  using lspg_stepper_t = typename lspg_problem_types::lspg_stepper_t;

  using eig_dyn_mat = Eigen::Matrix<scalar_t, -1, -1>;
  using hessian_t = rompp::algebra::Matrix<eig_dyn_mat>;
  using solver_tag = rompp::solvers::linear::iterative::LSCG;
  using lin_solver_t = rompp::solvers::iterative::EigenIterative<solver_tag, hessian_t>;
  lin_solver_t linSolverObj;

  // GN solver
  using converged_when_t = rompp::solvers::iterative::default_convergence;
  using gnsolver_t       = rompp::solvers::iterative::GaussNewton<
    lspg_stepper_t, converged_when_t, lin_solver_t>;
  gnsolver_t solver(lspgProblem.stepperObj_, yROM, linSolverObj);
  solver.setTolerance(1e-13);
  solver.setMaxIterations(200);

  // integrate in time
  rompp::ode::integrateNSteps(lspgProblem.stepperObj_, yROM, 0.0, dt, Nsteps, solver);

  // compute the fom corresponding to our rom final state
  auto yFomFinal = lspgProblem.yFomReconstructor_(yROM);
  yFomFinal.data()->Print(std::cout << std::setprecision(14));

  //----------------------------------------------------------------------
  // Compare the gold states from BDF1 to the reconstructed ROM
  //----------------------------------------------------------------------
  {
    using namespace rompp::apps::test;
    checkSol(rank, yFomFinal,
	     UnsteadyLinAdvDiff1dImplicitGoldStates<ode_case>
	     ::get(mu, domain, bc1D, dt, fint));
  }

  MPI_Finalize();
  std::cout << checkStr << std::endl;
  return 0;
}
